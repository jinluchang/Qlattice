#pragma once

#include <hash-cpp/crc32.h>

#include <qlat/config.h>
#include <qlat/utils.h>
#include <qlat/mpi.h>
#include <qlat/geometry.h>

#include <timer.h>

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
	if(get_id_node() == 0) *(p.os) << data;
	return p;
}
const rePort& operator<<(const rePort &p, std::ostream &(*func)(std::ostream&)){
	if(get_id_node() == 0) *(p.os) << func;
	return p;
}
static const rePort report;

inline std::string str_printf(const char *format, ...){
	char cstr[512];
	va_list args; va_start(args, format);
	vsnprintf(cstr, sizeof(cstr), format, args);
	return std::string(cstr);
}

inline int Printf(const char *format, ...){
	if(!get_id_node()){
		va_list args; va_start(args, format);
		return vprintf(format, args);
	}else{
		return 0;
	}
}

inline FILE* Fopen(const char* filename, const char* mode){
	if(!get_id_node()){
		return fopen(filename, mode);
	}else{
		return NULL;
	}
}

inline int Fprintf(FILE *pFile, const char *format, ...){
	if(!get_id_node()){
		va_list args; va_start(args, format);
		return vfprintf(pFile, format, args);
	}else{
		return 0;
	}
}

inline int Fflush(FILE *pFile){
	if(!get_id_node()){
		return fflush(pFile);
	}else{
		return 0;
	}
}

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
	for (long index = 0; index < geo.local_volume(); ++index) {
		Coordinate xl = geo.coordinate_from_index(index);
		const Vector<N> s = src.get_elems_const(xl);
		Vector<M> d = dest.get_elems(xl);
		for (int m = 0; m < geo.multiplicity; ++m) {
			castTruncated(d[m], s[m]);
		}
	}
}

template<class M>
uint32_t fieldChecksumSum32(const Field<M> &f)
{
	TIMER("fieldChecksumSum32");
	assert(f.geo.is_only_local());
	assert(sizeof(M) % sizeof(uint32_t) == 0);
	long sum = 0;
	const uint32_t *data = (const uint32_t *)f.field.data();
	const long size = f.field.size() * sizeof(M) / sizeof(uint32_t);
	for (long i = 0; i < size; ++i) {
		sum += data[i];
	}
	glb_sum(sum);
	uint32_t cs = sum;
	DisplayInfo(cname, fname, "check sum = %x\n", cs);
	return cs;
}

template<class M>
std::string field_hash_crc32(const qlat::Field<M> &origin){
	// somehow this checksum function does not agree with CPS's one.
	// Do not know why. But I am not sure what algorithm CPS uses.
	
	TIMER("field_hash_crc32");

	Geometry geo_only_local = geo_resize(origin.geo, 0);
	CRC32 crc32;
	void *buffer = (void *)&crc32;
	for(int id_node = 0; id_node < get_num_node(); id_node++){
		if(get_id_node() == id_node){
			crc32.add((void *)get_data(origin).data(), 
							get_data(origin).size() * sizeof(M));
		}
		sync_node();
		MPI_Bcast(buffer, 4, MPI_BYTE, id_node, get_comm());
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

	Geometry geo_only_local = geo_resize(origin.geo, 0);;

        Field<M> field_recv;
        field_recv.init(geo_only_local);

        Field<M> field_send;
        field_send.init(geo_only_local);
	field_send = origin;

	Field<M> field_rslt;
        field_rslt.init(geo_only_local);

	Coordinate total_site(geo_only_local.total_site(0),
			geo_only_local.total_site(1), 
			geo_only_local.total_site(2),
			geo_only_local.total_site(3));

	long range_low = geo_only_local.local_volume() * get_id_node();
	long range_high = range_low + geo_only_local.local_volume();

	for(int i = 0; i < get_num_node(); i++){
	
		// std::cout << "Shuffle loop: " << i << std::endl;
		
		int id_send_node = (get_id_node() + i) % get_num_node();
		
		Coordinate coor_send_node = qlat::coordinate_from_index(
									id_send_node, geo_only_local.geon.size_node);
#pragma omp parallel for
		for(int index = 0; index < geo_only_local.local_volume(); index++){
			Coordinate local_coor = geo_only_local.coordinate_from_index(index);
			Coordinate global_coor;
			for (int mu = 0; mu < 4; mu++) {
				global_coor[mu] = local_coor[mu] + coor_send_node[mu] 
													* geo_only_local.node_site[mu];
			}
			long global_index = index_from_coordinate(global_coor, total_site);
			if(global_index >= range_low && global_index < range_high)
			{
				Coordinate local_coor_write = geo_only_local.coordinate_from_index(
													global_index - range_low);
				assign(field_rslt.get_elems(local_coor_write), 
											field_send.get_elems_const(local_coor));
			}
		}

		get_data_dir(get_data(field_recv), get_data(field_send), 0);
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


	Geometry geo_only_local = geo_resize(origin.geo, 0);

        Field<M> field_recv;
        field_recv.init(geo_only_local);

        Field<M> field_send;
        field_send.init(geo_only_local);
	field_send = origin;

	FILE *outputFile = NULL;
	if(get_id_node() == 0){
		if(is_append) outputFile = fopen(write_addr.c_str(), "a");
		else outputFile = fopen(write_addr.c_str(), "w");
	}

	for(int i = 0; i < get_num_node(); i++){

		if(get_id_node() == 0){
			M *ptr = get_data(field_send).data(); 
			assert(ptr != NULL);
			long size = sizeof(M) * geo_only_local.local_volume() 
											* geo_only_local.multiplicity;
			std::cout << "Writing CYCLE: " << i << "\tSIZE = " << size << std::endl;
			timer_fwrite((char *)ptr, size, outputFile); 
			fflush(outputFile);
		}

		get_data_dir(get_data(field_recv), get_data(field_send), 0);
		swap(field_recv, field_send);
	}

	if(get_id_node() == 0) fclose(outputFile);
    
	report << "Export file CLOSED" << std::endl;

	sync_node();
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

	Geometry geo_only_local = geo_resize(destination.geo, 0);;

        Field<M> field_recv;
        field_recv.init(geo_only_local);

        Field<M> field_send;
        field_send.init(geo_only_local);

	Field<M> field_rslt;
        field_rslt.init(geo_only_local);

	Coordinate total_site(geo_only_local.total_site(0), geo_only_local.total_site(1), geo_only_local.total_site(2), geo_only_local.total_site(3));

	long range_low = geo_only_local.local_volume() * get_id_node();
	long range_high = range_low + geo_only_local.local_volume();

	
	// for every node:
	//
	FILE *inputFile = NULL;

// Well as you can see this is not really serial reading anymore. The sertial reading speed is unbearablly slow.
// Anyway it is tested. And it seems to be right.
	inputFile = fopen(read_addr.c_str(), "rb"); 
	assert(inputFile != NULL);
	assert(!ferror(inputFile));
	char line[1000];
	char indicator[] = "END_HEADER";

	int pos_ = -1;
	rewind(inputFile);
	while(fgets(line, 1000, inputFile) != NULL)
	{if(strstr(line, indicator) != NULL){
		pos_ = ftell(inputFile); break;
	}}
	assert(pos_ > -1); assert(!feof(inputFile));
		
	sync_node();
	int cycle_limit = 0;
	if(num_of_reading_threads > 0) 
		cycle_limit = (int)ceil((double)get_num_node() / num_of_reading_threads);
	else
		cycle_limit = 1;
	for(int cycle = 0; cycle < cycle_limit; cycle++){
		if(get_id_node() % cycle_limit == cycle){
			std::cout << "Reading STARTED:  Node Number =\t" 
				<< get_id_node() << std::endl;
			M *ptr = get_data(field_send).data();
			long size = sizeof(M) * geo_only_local.local_volume() 
													* geo_only_local.multiplicity;
			assert(!fseek(inputFile, size * get_id_node(), SEEK_CUR));
			timer_fread((char*)ptr, size, inputFile);
			std::cout << "Reading FINISHED: Node Number =\t" 
				<< get_id_node() << std::endl;
			fclose(inputFile);
		}
		sync_node();
	}
	
// 	if(get_id_node() == 0){
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
// 	for(int i = 0; i < get_num_node(); i++){
// 		
// 		if(get_id_node() == 0){
// 			std::cout << "Reading Cycle: " << i << std::endl;
// 			M *ptr = get_data(field_send).data();
//                         long size = sizeof(M) * geo_only_local.local_volume() * geo_only_local.multiplicity;
//                         timer_fread((char*)ptr, size, inputFile);
// 		// 	fflush(inputFile);
// 		}
// 		sync_node();	
// 		get_data_dir(get_data(field_recv), get_data(field_send), 0);
//                 swap(field_recv, field_send);
// 	}
// 
// 	if(get_id_node() == 0) fclose(inputFile);

	field_rslt = field_send;

	for(int i = 0; i < get_num_node(); i++){
		
		int id_send_node = (get_id_node() + i) % get_num_node();
		
		Coordinate coor_send_node = qlat::coordinate_from_index(id_send_node, geo_only_local.geon.size_node);
#pragma omp parallel for
		for(int index = 0; index < geo_only_local.local_volume(); index++){
			Coordinate local_coor = geo_only_local.coordinate_from_index(index);
			Coordinate global_coor;
			for (int mu = 0; mu < 4; mu++) {
				global_coor[mu] = local_coor[mu] + coor_send_node[mu] * geo_only_local.node_site[mu];
			}
			long global_index = index_from_coordinate(global_coor, total_site);
			if(global_index >= range_low && global_index < range_high)
			{
				Coordinate local_coor_read = geo_only_local.coordinate_from_index(global_index - range_low);
				assign(field_send.get_elems(local_coor), field_rslt.get_elems_const(local_coor_read));
			}
		}

		get_data_dir(get_data(field_recv), get_data(field_send), 0);
		swap(field_recv, field_send);
		if(get_id_node() == 0) std::cout << "Shuffling CYCLE:\t" << i << std::endl;
	}
	
	destination = field_send;

	sync_node();
}

QLAT_END_NAMESPACE

