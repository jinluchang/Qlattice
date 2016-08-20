#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>

#include <mpi.h>
#include <timer.h>

#include <array>
#include <map>
#include <vector>

QLAT_START_NAMESPACE

template <class M>
void fetch_expanded(Field<M> &field_comm){

	// tested for expansion = 2 case.

	TIMER_VERBOSE("fetch_expanded");

	std::map<Coordinate, std::vector<M> > send_map;
	std::map<Coordinate, int> send_map_consume;

	Coordinate pos; // coordinate position of a site relative to this node
	Coordinate local_pos; // coordinate position of a site relative to its home node
	Coordinate node_pos; // home node coordinate of a site in node space

// populate send_map with the data that we need to send to other nodes
	long record_size = field_comm.geo.localVolumeExpanded();
	for(long record = 0; record < record_size; record++){
		field_comm.geo.coordinateFromRecord(pos, record);
		if(field_comm.geo.isLocal(pos)) continue;
		for(int mu = 0; mu < DIM; mu++){
			local_pos[mu] = pos[mu] % field_comm.geo.nodeSite[mu];
			node_pos[mu] = pos[mu] / field_comm.geo.nodeSite[mu];
			if(local_pos[mu] < 0){
				local_pos[mu] += field_comm.geo.nodeSite[mu];
				node_pos[mu]--;
			}
		}
		std::vector<M> &vec = send_map[node_pos];
		for(int mu = 0; mu < field_comm.geo.multiplicity; mu++)
			vec.push_back(field_comm.getElemsConst(local_pos)[mu]);
	}

	std::vector<M> recv_vec;
	// will store data received from other nodes
	// Iterate over all the nodes to which we need to send data.
	// We ultimately copy the received data into the corresponding
	// value of sendmap.
	typename std::map<Coordinate, std::vector<M> >::iterator it;
	for(it = send_map.begin(); it != send_map.end(); it++){
		node_pos = it->first;
		std::vector<M> &send_vec = it->second;
		long size = send_vec.size();
		size_t size_bytes = size * sizeof(M);
		recv_vec.resize(std::max((long)2500, size));

		M *send = send_vec.data();
		M *recv = recv_vec.data();

		Coordinate coor_this, coort, coorf;
		int id_this, idt, idf;
		// assuming periodic boundary condition. maybe need some fixing?
		id_this = getIdNode();
		qlat::coordinateFromIndex(coor_this, id_this, \
			field_comm.geo.geon.sizeNode);
		coort = coor_this - node_pos; 
		regularize(coort, field_comm.geo.geon.sizeNode);
		coorf = coor_this + node_pos;
		regularize(coorf, field_comm.geo.geon.sizeNode);
		
		idt = qlat::indexFromCoordinate(coort, field_comm.geo.geon.sizeNode);
		idf = qlat::indexFromCoordinate(coorf, field_comm.geo.geon.sizeNode);
			
		MPI_Request req;
		MPI_Isend((void*)send, size_bytes, MPI_BYTE, idt, 0, getComm(), &req);
		const int ret = MPI_Recv((void*)recv, size_bytes, MPI_BYTE, \
			idf, 0, getComm(), MPI_STATUS_IGNORE);
		MPI_Wait(&req, MPI_STATUS_IGNORE);
		assert(!ret);

		memcpy(send, recv, size_bytes);

		send_map_consume[node_pos] = 0;
	}
	// Now send_map[node_pos] is the vector of data recieved from the node
	// pointed to by key.
	for(long record = 0; record < record_size; record++){
		field_comm.geo.coordinateFromRecord(pos, record);
		if(field_comm.geo.isLocal(pos)) continue;
		for(int mu = 0; mu < DIM; mu++){
			local_pos[mu] = pos[mu] % field_comm.geo.nodeSite[mu];
			node_pos[mu] = pos[mu] / field_comm.geo.nodeSite[mu];
			if(local_pos[mu] < 0){
				local_pos[mu] += field_comm.geo.nodeSite[mu];
				node_pos[mu]--;
			}
		}
		// send_map_consume[key] keeps track of our progress in consuming the
		// received data in sendmap[key], so that we know which offset of
		// send_map[node_pos] corresponds to which site.
		int consume = send_map_consume[node_pos];
		std::vector<M> &vec = send_map[node_pos];
		for(int mu = 0; mu < field_comm.geo.multiplicity; mu++){
			field_comm.getElems(pos)[mu] = vec[consume];
			consume++;
		}
		send_map_consume[node_pos] = consume;
	}
}

QLAT_END_NAMESPACE
