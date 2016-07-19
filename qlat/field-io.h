#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>
#include <qlat/mpi.h>
#include <qlat/geometry.h>

#include <omp.h>

#include <stdio.h>
#include <ctime>

QLAT_START_NAMESPACE       

template<class M>
void naive_serial_export(const qlat::Field<M> &origin, const std::string &export_addr){

	int MPI_rank_id;
	MPI_Comm_rank(getComm(), &MPI_rank_id);
	if(MPI_rank_id == 0) qlat::truncate(export_addr);

	std::ofstream output;
	
	Geometry geo_only_local;
	geo_only_local.init(origin.geo.geon, origin.geo.multiplicity, origin.geo.nodeSite);

	Field<M> field_only_local;
	field_only_local.init(geo_only_local);

	field_only_local = origin;

	for(int i = 0; i < geo_only_local.geon.numNode; i++){
		
		MPI_Comm_rank(getComm(), &MPI_rank_id);
		
		if(MPI_rank_id == i){
			
			std::cout << "Node ID: " << MPI_rank_id << std::endl;
		
			output.open(export_addr.c_str(), std::ios::app);	
			if(i == 0){
				time_t now = std::time(NULL);
				output << std::ctime(&now) << std::endl;
				output << "END HEADER" << std::endl;
			}
			M *ptr = getData(field_only_local).data();
			long size = sizeof(M) * geo_only_local.localVolume() * geo_only_local.multiplicity;
			std::cout << size << std::endl;
			output.write((char*)ptr, size);
			output.close();
		}

		syncNode();
	}
}

template<class M>
void naive_multiple_write(const qlat::Field<M> &origin, const std::string &export_addr){

	std::ofstream output;

	Geometry geo_only_local;
	geo_only_local.init(origin.geo.geon, origin.geo.multiplicity, origin.geo.nodeSite);

	Field<M> field_only_local;
	field_only_local.init(geo_only_local);

	field_only_local = origin;

	std::string real_export_addr = export_addr + "_node" + show((long)getIdNode());
	truncate(real_export_addr);

	output.open(real_export_addr.c_str());

	time_t now = std::time(NULL);
	output << std::ctime(&now) << std::endl;
	output << "END HEADER" << std::endl;

	M *ptr = getData(field_only_local).data();
	long size = sizeof(M) * geo_only_local.localVolume() * geo_only_local.multiplicity;
	output.write((char*)ptr, size);
	output.close();

}

// template<class M>
// void naive_multiple_read(const qlat::Field<M> &origin, const std::string &export_addr){
// 
// 	int MPI_rank_id;
// 	MPI_Comm_rank(getComm(), &MPI_rank_id);
// 
// 	std::ofstream output;
// 
// 	Geometry geo_only_local;
// 	geo_only_local.init(origin.geo.geon, origin.geo.multiplicity, origin.geo.nodeSite);
// 
// 	Field<M> field_only_local;
// 	field_only_local.init(geo_only_local);
// 
// 	field_only_local = origin;
// 
// 	std::string real_export_addr = export_addr + "_node" + show((long)getIdNode());
// 	truncate(real_export_addr);
// 
// 	output.open(real_export_addr.c_str());
// 
// 	time_t now = std::time(NULL);
// 	output << std::ctime(&now) << std::endl;
// 	output << "END HEADER" << std::endl;
// 
// 	M *ptr = getData(field_only_local).data();
// 	long size = sizeof(M) * geo_only_local.localVolume() * geo_only_local.multiplicity;
// 	output.write((char*)ptr, size);
// 	output.close();
// 
// }

template<class M>
void sophisticated_serial_write(const qlat::Field<M> &origin, 
		const std::string &write_addr, 
		const bool is_append = false,
		const bool does_skip_third = false){
	
	TIMER_FLOPS("sophisticated_serial_write");

	Geometry geo_only_local;
        geo_only_local.init(origin.geo.geon, origin.geo.multiplicity, origin.geo.nodeSite);

        Field<M> field_recv;
        field_recv.init(geo_only_local);

        Field<M> field_send;
        field_send.init(geo_only_local);
	field_send = origin;

	Field<M> field_rslt;
        field_rslt.init(geo_only_local);

	Coordinate totalSite(geo_only_local.totalSite(0),
			geo_only_local.totalSite(1), 
			geo_only_local.totalSite(2),
			geo_only_local.totalSite(3));

	long range_low = geo_only_local.localVolume() * getIdNode();
	long range_high = range_low + geo_only_local.localVolume();

	for(int i = 0; i < getNumNode(); i++){
	
		// std::cout << "Shuffle loop: " << i << std::endl;
		
		int id_send_node = (getIdNode() + i) % getNumNode();
		
		Coordinate coor_send_node; 
		qlat::coordinateFromIndex(coor_send_node, id_send_node, geo_only_local.geon.sizeNode);
#pragma omp parallel for
		for(int index = 0; index < geo_only_local.localVolume(); index++){
			Coordinate local_coor; geo_only_local.coordinateFromIndex(local_coor, index);
			Coordinate global_coor;
			for (int mu = 0; mu < 4; mu++) {
				global_coor[mu] = local_coor[mu] + coor_send_node[mu] * geo_only_local.nodeSite[mu];
			}
			long global_index = indexFromCoordinate(global_coor, totalSite);
			if(global_index >= range_low && global_index < range_high)
			{
				Coordinate local_coor_write; 
				geo_only_local.coordinateFromIndex(local_coor_write, global_index - range_low);
				assign(field_rslt.getElems(local_coor_write), field_send.getElemsConst(local_coor));
			}
		}

		getDataDir(getData(field_recv), getData(field_send), 0);
		swap(field_recv, field_send);	
	}
	
	field_send = field_rslt;

	FILE *outputFile = NULL;

	if(getIdNode() == 0){
		std::cout << "Node #0 open file!" << std::endl;
		if(is_append) outputFile = fopen(write_addr.c_str(), "a");
        	else outputFile = fopen(write_addr.c_str(), "w");
		assert(outputFile != NULL);

                std::ostringstream header_stream;
        
                header_stream << "BEGIN_HEADER" << std::endl;
                header_stream << "HDR_VERSION = 1.0" << std::endl;
                if(does_skip_third) header_stream << "DATATYPE = 4D_SU3_GAUGE" << std::endl;
                else header_stream << "DATATYPE = 4D_SU3_GAUGE_3x3" << std::endl;
                header_stream << "DIMENSION_1 = " << totalSite[0] << std::endl;
                header_stream << "DIMENSION_2 = " << totalSite[1] << std::endl;
                header_stream << "DIMENSION_3 = " << totalSite[2] << std::endl;
                header_stream << "DIMENSION_4 = " << totalSite[3] << std::endl;
                header_stream << "CHECKSUM = NOT yet implemented" << std::endl;
                header_stream << "LINK_TRACE =   NOT yet implemented" << std::endl;
                header_stream << "PLAQUETTE =   NOT yet implemented" << std::endl;
                header_stream << "CREATOR = RBC" << std::endl;
			time_t now = std::time(NULL);	
                header_stream << "ARCHIVE_DATE = " << std::ctime(&now);
                header_stream << "ENSEMBLE_LABEL = NOT yet implemented" << std::endl;
                header_stream << "FLOATING_POINT = IEEE64BIG" << std::endl;
                header_stream << "ENSEMBLE_ID = NOT yet implemented" << std::endl;
                header_stream << "SEQUENCE_NUMBER = NOT yet implemented" << std::endl;
                header_stream << "BETA = NOT yet implemented" << std::endl; 
                header_stream << "END_HEADER" << std::endl;
                
                fputs(header_stream.str().c_str(), outputFile);
	} 

	for(int i = 0; i < getNumNode(); i++){

		if(getIdNode() == 0){
			if(does_skip_third){
				char *ptr = (char*)(getData(field_send).data()); assert(ptr != NULL);
				char *cur = ptr;
				long size = geo_only_local.localVolume() * geo_only_local.multiplicity;
				int unit_size = sizeof(M) * 2 / 3;
				std::cout << "Writing Cycle: " << i << ", size = " << size << std::endl;
				for(long j = 0; j < size; j++){
					fwrite(cur, unit_size, 1, outputFile);
					cur += sizeof(M);
				}
			}else{
				M *ptr = getData(field_send).data(); assert(ptr != NULL);
				long size = sizeof(M) * geo_only_local.localVolume() * geo_only_local.multiplicity;
				std::cout << "Writing Cycle: " << i << ", size = " << size << std::endl;
				// std::cout << ((char*)ptr)[1] << std::endl;
				// output << "HAHHAHAHAHAHHA" << std::endl;
				// output.write((char*)ptr, 16);
				fwrite((char*)ptr, size, 1, outputFile); 
				fflush(outputFile);
			}
			// syncNode();
		}

		getDataDir(getData(field_recv), getData(field_send), 0);
		swap(field_recv, field_send);
	}

	if(getIdNode() == 0) fclose(outputFile);

	syncNode();
}

// std::string cps_Matrix_header_generator(const qlat::Field<cps::Matrix> &origin,
//                			const bool does_skip_third = false){
// 	NOT yet implemented :(
//	assert(false);	
//	return "NOT IMPLEMENTED.";
// }

void timer_fread(char* ptr, long size, FILE *inputFile){
	TIMER_VERBOSE("timer_free");
	fread(ptr, size, 1, inputFile);
}

template<class M>
void sophisticated_serial_read(qlat::Field<M> &destination, 
				const std::string &read_addr, 
				const int num_of_reading_threads = 0){

	// Tested to be working with 3x3 case.
	// Not yet able to read 3x2 file.
	// Not yet able to check checksum, plaquette, trace.
	// Just not able to do that. :(

	TIMER_VERBOSE("sophisticated_serial_read");

	Geometry geo_only_local;
        geo_only_local.init(destination.geo.geon, destination.geo.multiplicity, destination.geo.nodeSite);

        Field<M> field_recv;
        field_recv.init(geo_only_local);

        Field<M> field_send;
        field_send.init(geo_only_local);

	Field<M> field_rslt;
        field_rslt.init(geo_only_local);

	Coordinate totalSite(geo_only_local.totalSite(0), geo_only_local.totalSite(1), geo_only_local.totalSite(2), geo_only_local.totalSite(3));

	long range_low = geo_only_local.localVolume() * getIdNode();
	long range_high = range_low + geo_only_local.localVolume();

	
	// for every node:
	//
	FILE *inputFile = NULL;

// Well as you can see this is not really serial reading anymore. The sertial reading speed is unbearable slow.
// Anyway it is tested. And it seems to be right.
	inputFile = fopen(read_addr.c_str(), "rb"); assert(!ferror(inputFile));
	char line[1000];
	char indicator[] = "END_HEADER";

	int pos_ = -1;
	rewind(inputFile);
	while(fgets(line, 1000, inputFile) != NULL)
	{if(strstr(line, indicator) != NULL){
		pos_ = ftell(inputFile); break;
	}}
	assert(pos_ > -1); assert(!feof(inputFile));
		
	syncNode();
	int cycle_limit = 0;
	if(num_of_reading_threads > 0) 
		cycle_limit = (int)ceil((double)getNumNode() / num_of_reading_threads);
	else
		cycle_limit = 1;
	for(int cycle = 0; cycle < cycle_limit; cycle++){
		if(getIdNode() % cycle_limit == cycle){
			std::cout << "Reading Started: Node Number =\t" << getIdNode() << std::endl;
			M *ptr = getData(field_send).data();
			long size = sizeof(M) * geo_only_local.localVolume() * geo_only_local.multiplicity;
			assert(!fseek(inputFile, size * getIdNode(), SEEK_CUR));
			fread((char*)ptr, size, 1, inputFile);
			std::cout << "Reading Finished: Node Number =\t" << getIdNode() << std::endl;
			fclose(inputFile);
		}
		syncNode();
	}
	
// 	if(getIdNode() == 0){
//         	// input.open(read_addr.c_str());
// 		inputFile = fopen(read_addr.c_str(), "r");
// 		char line[1000];
// 		char indicator[] = "END_HEADER";
// 
// 		int pos_ = -1; fpos_t pos;
// 		rewind(inputFile);
// 		while(fgets(line, 1000, inputFile) != NULL)
// 		{if(strstr(line, indicator) != NULL){
// 			fgetpos(inputFile, &pos); pos_ = 1;  break;
// 		}}
// 		assert(pos_ > -1); assert(!feof(inputFile));
// 	} 
// 
// 	for(int i = 0; i < getNumNode(); i++){
// 		
// 		if(getIdNode() == 0){
// 			std::cout << "Reading Cycle: " << i << std::endl;
// 			M *ptr = getData(field_send).data();
//                         long size = sizeof(M) * geo_only_local.localVolume() * geo_only_local.multiplicity;
//                         timer_fread((char*)ptr, size, inputFile);
// 		// 	fflush(inputFile);
// 		}
// 		syncNode();	
// 		getDataDir(getData(field_recv), getData(field_send), 0);
//                 swap(field_recv, field_send);
// 	}
// 
// 	if(getIdNode() == 0) fclose(inputFile);

	field_rslt = field_send;

	for(int i = 0; i < getNumNode(); i++){
		
		int id_send_node = (getIdNode() + i) % getNumNode();
		
		Coordinate coor_send_node; 
		qlat::coordinateFromIndex(coor_send_node, id_send_node, geo_only_local.geon.sizeNode);
#pragma omp parallel for
		for(int index = 0; index < geo_only_local.localVolume(); index++){
			Coordinate local_coor; geo_only_local.coordinateFromIndex(local_coor, index);
			Coordinate global_coor;
			for (int mu = 0; mu < 4; mu++) {
				global_coor[mu] = local_coor[mu] + coor_send_node[mu] * geo_only_local.nodeSite[mu];
			}
			long global_index = indexFromCoordinate(global_coor, totalSite);
			if(global_index >= range_low && global_index < range_high)
			{
				Coordinate local_coor_read; 
				geo_only_local.coordinateFromIndex(local_coor_read, global_index - range_low);
				assign(field_send.getElems(local_coor), field_rslt.getElemsConst(local_coor_read));
			}
		}

		getDataDir(getData(field_recv), getData(field_send), 0);
		swap(field_recv, field_send);
		if(getIdNode() == 0) std::cout << "Shuffling Cycle:\t" << i << std::endl;
	}
	
	destination = field_send;

	syncNode();
}

QLAT_END_NAMESPACE

