#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>
#include <qlat/mpi.h>
#include <qlat/geometry.h>

#include <hash-cpp/crc32.h>

#include <omp.h>

#include <stdio.h>
#include <ctime>
#include <array>

#include <fstream>
#include <iostream>

QLAT_START_NAMESPACE       

typedef std::array<Complex, 6> MatrixTruncatedSU3;
typedef std::array<Complex, 9> MatrixSU3;

class rePort{
public:
	std::ostream *os;
	rePort(){
		os = &std::cout;
	}
};
template<class T>
const rePort& operator<<(const rePort &p, const T &data){
	if(getIdNode() == 0) *(p.os) << data;
	return p;
}
const rePort& operator<<(const rePort &p, std::ostream &(*func)(std::ostream&)){
	if(getIdNode() == 0) *(p.os) << func;
	return p;
}
static const rePort report;

template<class M, class N>
void castTruncated(M &x, const N &y)
{
	assert(sizeof(M) <= sizeof(N));
	memcpy(&x, &y, sizeof(M));
}

template<class M, class N>
void fieldCastTruncated(Field<M> &dest, const Field<N> &src)
{
	TIMER("fieldCastTruncated");
	const Geometry& geo = src.geo;
	dest.init(geo);
#pragma omp parallel for
	for (long index = 0; index < geo.localVolume(); ++index) {
		Coordinate xl; geo.coordinateFromIndex(xl, index);
		const Vector<N> s = src.getElemsConst(xl);
		Vector<M> d = dest.getElems(xl);
		for (int m = 0; m < geo.multiplicity; ++m) {
			castTruncated(d[m], s[m]);
		}
	}
}

template<class M>
uint32_t fieldChecksumSum32(const Field<M> &f)
{
	TIMER("fieldChecksumSum32");
	assert(f.geo.isOnlyLocal());
	assert(sizeof(M) % sizeof(uint32_t) == 0);
	long sum = 0;
	const uint32_t *data = (const uint32_t *)f.field.data();
	const long size = f.field.size() * sizeof(M) / sizeof(uint32_t);
	for (long i = 0; i < size; ++i) {
		sum += data[i];
	}
	glbSum(sum);
	uint32_t cs = sum;
	DisplayInfo(cname, fname, "check sum = %x\n", cs);
	return cs;
}

template<class M>
std::string field_hash_crc32(const qlat::Field<M> &origin){
	// somehow this checksum function does not agree with CPS's one.
	// Do not know why. But I am not sure what algorithm CPS uses.
	
	TIMER("field_hash_crc32");

	Geometry geo_only_local;
        geo_only_local.init(origin.geo.geon, \
		origin.geo.multiplicity, origin.geo.nodeSite);
	CRC32 crc32;
	void *buffer = (void *)&crc32;
	for(int id_node = 0; id_node < getNumNode(); id_node++){
		if(getIdNode() == id_node){
			crc32.add((void *)getData(origin).data(), \
				getData(origin).size() * sizeof(M));
		}
		syncNode();
		MPI_Bcast(buffer, 4, MPI_BYTE, id_node, getComm());
	}
	return crc32.getHash();			
}

void timer_fwrite(char* ptr, long size, FILE *outputFile){
	TIMER("timer_fwrite");
	fwrite(ptr, size, 1, outputFile);
}

template<class M>
void sophisticated_make_to_order(Field<M> &result, const Field<M> &origin){
	TIMER("sophisticated_make_to_order");

	Geometry geo_only_local;
        geo_only_local.init(origin.geo.geon, \
		origin.geo.multiplicity, origin.geo.nodeSite);

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
		qlat::coordinateFromIndex(coor_send_node, \
			id_send_node, geo_only_local.geon.sizeNode);
#pragma omp parallel for
		for(int index = 0; index < geo_only_local.localVolume(); index++){
			Coordinate local_coor; 
			geo_only_local.coordinateFromIndex(local_coor, index);
			Coordinate global_coor;
			for (int mu = 0; mu < 4; mu++) {
				global_coor[mu] = local_coor[mu] + coor_send_node[mu] \
					* geo_only_local.nodeSite[mu];
			}
			long global_index = indexFromCoordinate(global_coor, totalSite);
			if(global_index >= range_low && global_index < range_high)
			{
				Coordinate local_coor_write; 
				geo_only_local.coordinateFromIndex(local_coor_write, \
					global_index - range_low);
				assign(field_rslt.getElems(local_coor_write), \
					field_send.getElemsConst(local_coor));
			}
		}

		getDataDir(getData(field_recv), getData(field_send), 0);
		swap(field_recv, field_send);	
	}
	result.init(geo_only_local);
	result = field_rslt;
}

template<class M>
void sophisticated_serial_write(const qlat::Field<M> &origin, 
		const std::string &write_addr, 
		const bool is_append = false){
	
	TIMER("sophisticated_serial_write");


	Geometry geo_only_local;
        geo_only_local.init(origin.geo.geon, \
		origin.geo.multiplicity, origin.geo.nodeSite);

        Field<M> field_recv;
        field_recv.init(geo_only_local);

        Field<M> field_send;
        field_send.init(geo_only_local);
	field_send = origin;

	FILE *outputFile = NULL;
	if(getIdNode() == 0){
		if(is_append) outputFile = fopen(write_addr.c_str(), "a");
		else outputFile = fopen(write_addr.c_str(), "w");
	}

	for(int i = 0; i < getNumNode(); i++){

		if(getIdNode() == 0){
			M *ptr = getData(field_send).data(); 
			assert(ptr != NULL);
			long size = sizeof(M) * geo_only_local.localVolume() * \
				    geo_only_local.multiplicity;
			std::cout << "Writing Cycle: " << i << "\tsize = " << size << std::endl;
			timer_fwrite((char *)ptr, size, outputFile); 
			fflush(outputFile);
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
	TIMER_VERBOSE("timer_fread");
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

// Well as you can see this is not really serial reading anymore. The sertial reading speed is unbearablly slow.
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
			std::cout << "Reading Started: Node Number =\t" 
				<< getIdNode() << std::endl;
			M *ptr = getData(field_send).data();
			long size = sizeof(M) * geo_only_local.localVolume() * \
				geo_only_local.multiplicity;
			assert(!fseek(inputFile, size * getIdNode(), SEEK_CUR));
			timer_fread((char*)ptr, size, inputFile);
			std::cout << "Reading Finished: Node Number =\t" 
				<< getIdNode() << std::endl;
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

