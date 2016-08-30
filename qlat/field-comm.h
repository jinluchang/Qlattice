#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>

#include <mpi.h>
#include <timer.h>

#include <array>
#include <map>
#include <set>
#include <vector>

QLAT_START_NAMESPACE

template <class M>
void fetch_expanded(Field<M> &field_comm){

	// tested for expansion = 2 case.

	TIMER("fetch_expanded");

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
	
	// pure communication
	
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

template<class M>
class Chart: public std::map<Coordinate, std::vector<Coordinate> > {
public:
	Coordinate expansionLeft, expansionRight;
	Geometry geo;
	std::map<Coordinate, std::vector<M> > send_map;
};

enum gActionType {WILSON, IWASAKI};
class gAction{
public:
	gActionType type;
	double c1;
	gAction(){c1 = 0.; type = WILSON;}
};

template<class M>
void produce_chart_envelope(Chart<M> &chart, const Geometry geometry, const gAction &gA){
	TIMER("produce_chart_envelope()");
	
	chart.geo = geometry;
	std::set<Coordinate> target;

	int muP, nuP;
	switch(gA.type){
		case WILSON: 	muP = 1;
				nuP = 1;
				chart.expansionLeft = Coordinate(1, 1, 1, 1);
				chart.expansionRight = Coordinate(1, 1, 1, 1);
				break;
		case IWASAKI:	muP = 2;
				nuP = 1;
				chart.expansionLeft = Coordinate(2, 2, 2, 2);
				chart.expansionRight = Coordinate(2, 2, 2, 2);
				break;
		default:	assert(false);
	}

	Coordinate index_pos;
	Coordinate index_pos_m;
	for(long index = 0; index < geometry.localVolume(); index++){
		geometry.coordinateFromIndex(index_pos, index);
		for(int mu = 0; mu < DIM; mu++){
		for(int nu = 0; nu < DIM; nu++){
			if(mu == nu) continue;
			for(int muI = -muP; muI <= muP; muI++){
			for(int nuI = -nuP; nuI <= nuP; nuI++){
				index_pos_m = index_pos;
				index_pos_m[mu] += muI;
				index_pos_m[nu] += nuI;
				if(!chart.geo.isLocal(index_pos_m))
					target.insert(index_pos_m);
			}}
		}}
	}

	Coordinate pos; // coordinate position of a site relative to this node
	Coordinate local_pos; // coordinate position of a site relative to its home node
	Coordinate node_pos; // home node coordinate of a site in node space

	chart.clear();
	std::set<Coordinate>::const_iterator it;
	for(it = target.begin(); it != target.end(); it++){
		pos = *it;
		for(int mu = 0; mu < DIM; mu++){
			local_pos[mu] = pos[mu] % geometry.nodeSite[mu];
			node_pos[mu] = pos[mu] / geometry.nodeSite[mu];
			if(local_pos[mu] < 0){
				local_pos[mu] += geometry.nodeSite[mu];
				node_pos[mu]--;
			}
		}
		chart[node_pos].push_back(local_pos);
	}

	typename Chart<M>::iterator it_chart;
	long size;
	for(it_chart = chart.begin(); it_chart != chart.end(); it_chart++){
		size = geometry.multiplicity * it_chart->second.size();
		chart.send_map[it_chart->first].resize(size);
	}
}

template<class M>
void produce_chart_geo(Chart<M> &chart, const Geometry geometry){
	
	Coordinate pos; // coordinate position of a site relative to this node
	Coordinate local_pos; // coordinate position of a site relative to its home node
	Coordinate node_pos; // home node coordinate of a site in node space
	
	chart.geo = geometry;
	chart.expansionLeft = geometry.expansionLeft;
	chart.expansionRight = geometry.expansionRight;

	chart.clear();
	long record_size = geometry.localVolumeExpanded();
	for(long record = 0; record < record_size; record++){
		geometry.coordinateFromRecord(pos, record);
		if(geometry.isLocal(pos)) continue;
		for(int mu = 0; mu < DIM; mu++){
			local_pos[mu] = pos[mu] % geometry.nodeSite[mu];
			node_pos[mu] = pos[mu] / geometry.nodeSite[mu];
			if(local_pos[mu] < 0){
				local_pos[mu] += geometry.nodeSite[mu];
				node_pos[mu]--;
			}
		}
		chart[node_pos].push_back(local_pos);
	}
	
	typename Chart<M>::iterator it_chart;
	long size;
	for(it_chart = chart.begin(); it_chart != chart.end(); it_chart++){
		size = geometry.multiplicity * it_chart->second.size();
		chart.send_map[it_chart->first].resize(size);
	}
}

template <class M>
void fetch_expanded_chart(Field<M> &field_comm, Chart<M> &send_chart){
	TIMER("fetch_expanded_chart");

	assert(isMatchingGeo(send_chart.geo, field_comm.geo));

	Coordinate node_pos; // home node coordinate of a site in node space

// std::map<Coordinate, std::vector<Coordinate> >
//           node_pos          local_pos
//        it_chart->first  it_chart->second
//                                 ^
//                                 |
//                              it_coor

// populate send_map with the data that we need to send to other nodes
	typename Chart<M>::const_iterator it_chart;
	std::vector<Coordinate>::const_iterator it_coor;
	for(it_chart = send_chart.begin(); it_chart != send_chart.end(); it_chart++){
		node_pos = it_chart->first;
		long consume = 0;
		std::vector<M> &vec = send_chart.send_map[node_pos];
		for(it_coor = it_chart->second.begin(); 
					it_coor != it_chart->second.end(); it_coor++){
			for(int mu = 0; mu < field_comm.geo.multiplicity; mu++){
				vec[consume] = field_comm.getElemsConst(*it_coor)[mu];
				consume++;
			}
		}
	}

	static std::vector<M> recv_vec;
// will store data received from other nodes
// Iterate over all the nodes to which we need to send data.
// We ultimately copy the received data into the corresponding
// value of sendmap.
	typename std::map<Coordinate, std::vector<M> >::iterator it;
	
	{
	TIMER("pure_comm.");
	for(it = send_chart.send_map.begin(); it != send_chart.send_map.end(); it++){
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
	}
	
	}
	Coordinate pos;
// Now send_map[node_pos] is the vector of data recieved from the node
// pointed to by key.
	for(it_chart = send_chart.begin(); it_chart != send_chart.end(); it_chart++){
		node_pos = it_chart->first;
		long consume = 0;
		std::vector<M> &vec = send_chart.send_map[node_pos];
		for(it_coor = it_chart->second.begin(); 
					it_coor != it_chart->second.end(); it_coor++){
			for(int mu = 0; mu < field_comm.geo.multiplicity; mu++){
				pos = node_pos * field_comm.geo.nodeSite + *it_coor;
				field_comm.getElems(pos)[mu] = vec[consume];
				consume++;
			}
		}
	}

}

QLAT_END_NAMESPACE
