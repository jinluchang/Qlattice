#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>
#include <qlat/mpi.h>
#include <qlat/geometry.h>

#include <omp.h>

#include <stdio.h>

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
void naive_multiple_export(const qlat::Field<M> &origin, const std::string &export_addr){

	int MPI_rank_id;
	MPI_Comm_rank(getComm(), &MPI_rank_id);

	std::ofstream output;

	Geometry geo_only_local;
	geo_only_local.init(origin.geo.geon, origin.geo.multiplicity, origin.geo.nodeSite);

	Field<M> field_only_local;
	field_only_local.init(geo_only_local);

	field_only_local = origin;

	std::string real_export_addr = export_addr + "_node" + show((long)MPI_rank_id);
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

	Coordinate totalSite(geo_only_local.totalSite(0), geo_only_local.totalSite(1), geo_only_local.totalSite(2), geo_only_local.totalSite(3));

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
                header_stream << "CHECKSUM = f1bdcf61" << std::endl;
                header_stream << "LINK_TRACE =   2.606798585E-05" << std::endl;
                header_stream << "PLAQUETTE =   5.880740538E-01" << std::endl;
                header_stream << "CREATOR = RBC" << std::endl;
                header_stream << "ARCHIVE_DATE = Thu Jul 14 15:44:49 2016" << std::endl;
                header_stream << "ENSEMBLE_LABEL = IWASAKI_Nf2p1_24c64s16_b2.13M1.80mu0.005ms0.04_rhm" << std::endl;
                header_stream << "FLOATING_POINT = IEEE64BIG" << std::endl;
                header_stream << "ENSEMBLE_ID = 147" << std::endl;
                header_stream << "SEQUENCE_NUMBER = 5000" << std::endl;
                header_stream << "BETA = 2.13" << std::endl; 
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
				std::cout << "Write cycle: " << i << ", size = " << size << std::endl;
				for(long j = 0; j < size; j++){
					fwrite(cur, unit_size, 1, outputFile);
					cur += sizeof(M);
				}
			}else{
				M *ptr = getData(field_send).data(); assert(ptr != NULL);
				long size = sizeof(M) * geo_only_local.localVolume() * geo_only_local.multiplicity;
				std::cout << "Write cycle: " << i << ", size = " << size << std::endl;
				// std::cout << ((char*)ptr)[1] << std::endl;
				// output << "HAHHAHAHAHAHHA" << std::endl;
				// output.write((char*)ptr, 16);
				fwrite((char*)ptr, size, 1, outputFile); 
				fflush(outputFile);
			}

		}

		getDataDir(getData(field_recv), getData(field_send), 0);
		swap(field_recv, field_send);
	}

	if(getIdNode() == 0) fclose(outputFile);

	syncNode();
}

template<class M>
void sophisticated_serial_read(qlat::Field<M> &destination, const std::string &read_addr){
	
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

	std::ifstream input;
	if(getIdNode() == 0){
        	input.open(read_addr.c_str());
	} 

	for(int i = 0; i < getNumNode(); i++){
		
		if(getIdNode() == 0){
			M *ptr = getData(field_send).data();
                        long size = sizeof(M) * geo_only_local.localVolume() * geo_only_local.multiplicity;
                        input.read((char*)ptr, size);
		}
	
		getDataDir(getData(field_recv), getData(field_send), 0);
                swap(field_recv, field_send);
	}

	if(getIdNode() == 0) input.close();

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
	}
	
	destination = field_send;

	syncNode();
}

QLAT_END_NAMESPACE

