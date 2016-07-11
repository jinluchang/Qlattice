#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>
#include <qlat/mpi.h>
#include <qlat/geometry.h>

QLAT_START_NAMESPACE       

template<class M>
void naive_serial_export(const Field<M> &origin, const std::string &export_addr){

	int MPI_rank_id;
	MPI_Comm_rank(getComm(), &MPI_rank_id);
	if(MPI_rank_id == 0) qlat::truncate(export_addr);

	std::ofstream output;
	
	geometry geo_only_local;
	geo_only_local.init(origin.geon, origin.geo.multiplicity, origin.geo.nodeSite);

	Field<M> field_only_local;
	field_only_local.init(geo_only_local);

	field_only_local = origin;

	for(int i = 0; i < geo.geon.numNode; i++){
		
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
			output.write((char*)ptr, size);
			output.close();
		}

		syncNode();
	}
}


QLAT_END_NAMESPACE

