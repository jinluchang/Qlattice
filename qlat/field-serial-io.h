#pragma once

#include <qlat/field-io.h>

#include <omp.h>

#include <stdio.h>
#include <ctime>
#include <array>

#include <fstream>
#include <iostream>

QLAT_START_NAMESPACE       

template<class M>
void field_import_serial(qlat::Field<M>& destination, 
    const std::string& read_addr, 
    const long offset = 0,
    const int whence = SEEK_SET,
    const int num_of_reading_threads = 0)
{
  TIMER_VERBOSE("field_import_serial");

  Geometry geo_only_local = geo_resize(destination.geo, 0);;

  Field<M> field_recv;
  field_recv.init(geo_only_local);

  Field<M> field_send;
  field_send.init(geo_only_local);

  Field<M> field_rslt;
  field_rslt.init(geo_only_local);

  Coordinate total_site = geo_only_local.total_site();

  long range_low = geo_only_local.local_volume() * get_id_node();
  long range_high = range_low + geo_only_local.local_volume();

  // for every node:
  //
  FILE *inputFile = NULL;

  // Well as you can see this is not really serial reading anymore. The sertial reading speed is unbearablly slow.
  // Anyway it is tested. And it seems to be right.
  inputFile = fopen(read_addr.c_str(), "r"); 
  qassert(inputFile != NULL);
  qassert(!ferror(inputFile));
  fseek(inputFile, offset, whence);
  qassert(!feof(inputFile));

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
  // 		qassert(pos_ > -1); qassert(!feof(inputFile));
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
