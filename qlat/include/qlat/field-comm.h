#pragma once

#include <qlat/setup.h>
#include <qlat/mpi.h>

#include <map>
#include <set>
#include <vector>

namespace qlat { //

template <class M>
void fetch_expanded(Field<M> &field_comm)
{
  // tested for expansion = 2 case.

  TIMER("fetch_expanded");

  std::map<Coordinate, std::vector<M> > send_map;
  std::map<Coordinate, int> send_map_consume;

  Coordinate pos;  // coordinate position of a site relative to this node
  Coordinate
      local_pos;  // coordinate position of a site relative to its home node
  Coordinate node_pos;  // home node coordinate of a site in node space

  // populate send_map with the data that we need to send to other nodes
  long record_size = field_comm.geo().local_volume_expanded();
  for (long record = 0; record < record_size; record++) {
    pos = field_comm.geo().coordinateFromRecord(record);
    if (field_comm.geo().is_local(pos)) continue;
    for (int mu = 0; mu < DIMN; mu++) {
      local_pos[mu] = pos[mu] % field_comm.geo().node_site[mu];
      node_pos[mu] = pos[mu] / field_comm.geo().node_site[mu];
      if (local_pos[mu] < 0) {
        local_pos[mu] += field_comm.geo().node_site[mu];
        node_pos[mu]--;
      }
    }
    std::vector<M> &vec = send_map[node_pos];
    for (int mu = 0; mu < field_comm.geo().multiplicity; mu++)
      vec.push_back(field_comm.get_elems_const(local_pos)[mu]);
  }

  std::vector<M> recv_vec;
  // will store data received from other nodes
  // Iterate over all the nodes to which we need to send data.
  // We ultimately copy the received data into the corresponding
  // value of sendmap.
  typename std::map<Coordinate, std::vector<M> >::iterator it;

  // pure communication

  for (it = send_map.begin(); it != send_map.end(); it++) {
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
    id_this = get_id_node();
    coor_this =
        qlat::coordinate_from_index(id_this, field_comm.geo().geon.size_node);
    coort = coor_this - node_pos;
    regularize_coordinate(coort, field_comm.geo().geon.size_node);
    coorf = coor_this + node_pos;
    regularize_coordinate(coorf, field_comm.geo().geon.size_node);

    idt = qlat::index_from_coordinate(coort, field_comm.geo().geon.size_node);
    idf = qlat::index_from_coordinate(coorf, field_comm.geo().geon.size_node);

    MPI_Request req;
    MPI_Isend((void *)send, size_bytes, MPI_BYTE, idt, 0, get_comm(), &req);
    const int ret = MPI_Recv((void *)recv, size_bytes, MPI_BYTE, idf, 0,
                             get_comm(), MPI_STATUS_IGNORE);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    qassert(!ret);

    memcpy(send, recv, size_bytes);

    send_map_consume[node_pos] = 0;
  }
  // Now send_map[node_pos] is the vector of data recieved from the node
  // pointed to by key.
  for (long record = 0; record < record_size; record++) {
    pos = field_comm.geo().coordinateFromRecord(record);
    if (field_comm.geo().is_local(pos)) continue;
    for (int mu = 0; mu < DIMN; mu++) {
      local_pos[mu] = pos[mu] % field_comm.geo().node_site[mu];
      node_pos[mu] = pos[mu] / field_comm.geo().node_site[mu];
      if (local_pos[mu] < 0) {
        local_pos[mu] += field_comm.geo().node_site[mu];
        node_pos[mu]--;
      }
    }
    // send_map_consume[key] keeps track of our progress in consuming the
    // received data in sendmap[key], so that we know which offset of
    // send_map[node_pos] corresponds to which site.
    int consume = send_map_consume[node_pos];
    std::vector<M> &vec = send_map[node_pos];
    for (int mu = 0; mu < field_comm.geo().multiplicity; mu++) {
      field_comm.get_elems(pos)[mu] = vec[consume];
      consume++;
    }
    send_map_consume[node_pos] = consume;
  }
}

template <class M>
class Chart : public std::map<Coordinate, std::vector<Coordinate> >
{
 public:
  Coordinate expansion_left, expansion_right;
  Geometry geo;
  std::map<Coordinate, std::vector<M> > send_map;
};

enum GAUGE_TYPE { WILSON, IWASAKI };
class Gauge
{
 public:
  GAUGE_TYPE type;
  double c1;
  Gauge()
  {
    c1 = 0.;
    type = WILSON;
  }
};

template <class M>
void produce_chart_envelope(Chart<M> &chart, const Geometry geometry,
                            const Gauge &gauge)
{
  TIMER("produce_chart_envelope()");

  chart.geo() = geometry;
  std::set<Coordinate> target;

  int muP, nuP;
  switch (gauge.type) {
    case WILSON:
      muP = 1;
      nuP = 1;
      chart.expansion_left = Coordinate(1, 1, 1, 1);
      chart.expansion_right = Coordinate(1, 1, 1, 1);
      break;
    case IWASAKI:
      muP = 2;
      nuP = 1;
      chart.expansion_left = Coordinate(2, 2, 2, 2);
      chart.expansion_right = Coordinate(2, 2, 2, 2);
      break;
    default:
      qassert(false);
  }

  Coordinate index_pos;
  Coordinate index_pos_m;
  for (long index = 0; index < geometry.local_volume(); index++) {
    index_pos = geometry.coordinate_from_index(index);
    for (int mu = 0; mu < DIMN; mu++) {
      for (int nu = 0; nu < DIMN; nu++) {
        if (mu == nu) continue;
        for (int muI = -muP; muI <= muP; muI++) {
          for (int nuI = -nuP; nuI <= nuP; nuI++) {
            index_pos_m = index_pos;
            index_pos_m[mu] += muI;
            index_pos_m[nu] += nuI;
            if (!chart.geo().is_local(index_pos_m)) target.insert(index_pos_m);
          }
        }
      }
    }
  }

  Coordinate pos;  // coordinate position of a site relative to this node
  Coordinate
      local_pos;  // coordinate position of a site relative to its home node
  Coordinate node_pos;  // home node coordinate of a site in node space

  chart.clear();
  std::set<Coordinate>::const_iterator it;
  for (it = target.begin(); it != target.end(); it++) {
    pos = *it;
    for (int mu = 0; mu < DIMN; mu++) {
      local_pos[mu] = pos[mu] % geometry.node_site[mu];
      node_pos[mu] = pos[mu] / geometry.node_site[mu];
      if (local_pos[mu] < 0) {
        local_pos[mu] += geometry.node_site[mu];
        node_pos[mu]--;
      }
    }
    chart[node_pos].push_back(local_pos);
  }

  typename Chart<M>::iterator it_chart;
  long size;
  for (it_chart = chart.begin(); it_chart != chart.end(); it_chart++) {
    size = geometry.multiplicity * it_chart->second.size();
    chart.send_map[it_chart->first].resize(size);
  }
}

// TODO: FIXME!!!
// template<class M>
// void produce_chart_envelope(Chart<M> &chart, const Geometry geometry,
//								array<int, DIMN - 1> &R,
//int &T){ 	TIMER("produce_chart_envelope()");
//
//	chart.geo() = geometry;
//	std::set<Coordinate> target;
//
//	int muP, nuP;
//	switch(gauge.type){
//		case WILSON: 	muP = 1;
//				nuP = 1;
//				chart.expansion_left = Coordinate(1, 1, 1, 1);
//				chart.expansion_right = Coordinate(1, 1, 1, 1);
//				break;
//		case IWASAKI:	muP = 2;
//				nuP = 1;
//				chart.expansion_left = Coordinate(2, 2, 2, 2);
//				chart.expansion_right = Coordinate(2, 2, 2, 2);
//				break;
//		default:	qassert(false);
//	}
//
//	Coordinate index_pos;
//	Coordinate index_pos_m;
//	for(long index = 0; index < geometry.local_volume(); index++){
//		geometry.coordinate_from_index(index_pos, index);
//		for(int mu = 0; mu < DIMN; mu++){
//		for(int nu = 0; nu < DIMN; nu++){
//			if(mu == nu) continue;
//			for(int muI = -muP; muI <= muP; muI++){
//			for(int nuI = -nuP; nuI <= nuP; nuI++){
//				index_pos_m = index_pos;
//				index_pos_m[mu] += muI;
//				index_pos_m[nu] += nuI;
//				if(!chart.geo().is_local(index_pos_m))
//					target.insert(index_pos_m);
//			}}
//		}}
//	}
//
//	Coordinate pos; // coordinate position of a site relative to this node
//	Coordinate local_pos; // coordinate position of a site relative to its
// home node 	Coordinate node_pos; // home node coordinate of a site in node
// space
//
//	chart.clear();
//	std::set<Coordinate>::const_iterator it;
//	for(it = target.begin(); it != target.end(); it++){
//		pos = *it;
//		for(int mu = 0; mu < DIMN; mu++){
//			local_pos[mu] = pos[mu] % geometry.node_site[mu];
//			node_pos[mu] = pos[mu] / geometry.node_site[mu];
//			if(local_pos[mu] < 0){
//				local_pos[mu] += geometry.node_site[mu];
//				node_pos[mu]--;
//			}
//		}
//		chart[node_pos].push_back(local_pos);
//	}
//
//	typename Chart<M>::iterator it_chart;
//	long size;
//	for(it_chart = chart.begin(); it_chart != chart.end(); it_chart++){
//		size = geometry.multiplicity * it_chart->second.size();
//		chart.send_map[it_chart->first].resize(size);
//	}
//}

template <class M>
void produce_chart_geo(Chart<M> &chart, const Geometry geometry)
{
  Coordinate pos;  // coordinate position of a site relative to this node
  Coordinate
      local_pos;  // coordinate position of a site relative to its home node
  Coordinate node_pos;  // home node coordinate of a site in node space

  chart.geo() = geometry;
  chart.expansion_left = geometry.expansion_left;
  chart.expansion_right = geometry.expansion_right;

  chart.clear();
  long record_size = geometry.local_volume_expanded();
  for (long record = 0; record < record_size; record++) {
    pos = geometry.coordinateFromRecord(record);
    if (geometry.is_local(pos)) continue;
    for (int mu = 0; mu < DIMN; mu++) {
      local_pos[mu] = pos[mu] % geometry.node_site[mu];
      node_pos[mu] = pos[mu] / geometry.node_site[mu];
      if (local_pos[mu] < 0) {
        local_pos[mu] += geometry.node_site[mu];
        node_pos[mu]--;
      }
    }
    chart[node_pos].push_back(local_pos);
  }

  typename Chart<M>::iterator it_chart;
  long size;
  for (it_chart = chart.begin(); it_chart != chart.end(); it_chart++) {
    size = geometry.multiplicity * it_chart->second.size();
    chart.send_map[it_chart->first].resize(size);
  }
}

template <class M>
void fetch_expanded_chart(Field<M> &field_comm, Chart<M> &send_chart)
{
  TIMER("fetch_expanded_chart");

  qassert(is_matching_geo_mult(send_chart.geo(), field_comm.geo()));

  Coordinate node_pos;  // home node coordinate of a site in node space

  // std::map<Coordinate, std::vector<Coordinate> >
  //           node_pos          local_pos
  //        it_chart->first  it_chart->second
  //                                 ^
  //                                 |
  //                              it_coor

  // populate send_map with the data that we need to send to other nodes
  typename Chart<M>::const_iterator it_chart;
  std::vector<Coordinate>::const_iterator it_coor;
  for (it_chart = send_chart.begin(); it_chart != send_chart.end();
       it_chart++) {
    node_pos = it_chart->first;
    long consume = 0;
    std::vector<M> &vec = send_chart.send_map[node_pos];
    for (it_coor = it_chart->second.begin(); it_coor != it_chart->second.end();
         it_coor++) {
      for (int mu = 0; mu < field_comm.geo().multiplicity; mu++) {
        vec[consume] = field_comm.get_elems_const(*it_coor)[mu];
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
    for (it = send_chart.send_map.begin(); it != send_chart.send_map.end();
         it++) {
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
      id_this = get_id_node();
      coor_this =
          qlat::coordinate_from_index(id_this, field_comm.geo().geon.size_node);
      coort = coor_this - node_pos;
      regularize_coordinate(coort, field_comm.geo().geon.size_node);
      coorf = coor_this + node_pos;
      regularize_coordinate(coorf, field_comm.geo().geon.size_node);

      idt = qlat::index_from_coordinate(coort, field_comm.geo().geon.size_node);
      idf = qlat::index_from_coordinate(coorf, field_comm.geo().geon.size_node);

      MPI_Request req;
      MPI_Isend((void *)send, size_bytes, MPI_BYTE, idt, 0, get_comm(), &req);
      const int ret = MPI_Recv((void *)recv, size_bytes, MPI_BYTE, idf, 0,
                               get_comm(), MPI_STATUS_IGNORE);
      MPI_Wait(&req, MPI_STATUS_IGNORE);
      qassert(!ret);

      memcpy(send, recv, size_bytes);
    }
  }
  Coordinate pos;
  // Now send_map[node_pos] is the vector of data recieved from the node
  // pointed to by key.
  for (it_chart = send_chart.begin(); it_chart != send_chart.end();
       it_chart++) {
    node_pos = it_chart->first;
    long consume = 0;
    std::vector<M> &vec = send_chart.send_map[node_pos];
    for (it_coor = it_chart->second.begin(); it_coor != it_chart->second.end();
         it_coor++) {
      for (int mu = 0; mu < field_comm.geo().multiplicity; mu++) {
        pos = node_pos * field_comm.geo().node_site + *it_coor;
        field_comm.get_elems(pos)[mu] = vec[consume];
        consume++;
      }
    }
  }
}

}
