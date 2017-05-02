// vim: set ts=2 sw=2 expandtab:

#pragma once

#include <qlat/config.h>
#include <qlat/utils.h>

#include <mpi.h>

#include <array>
#include <map>
#include <set>
#include <vector>

QLAT_START_NAMESPACE

template <class M>
void refresh_expanded(Field<M>& field_comm)
{
  TIMER("refresh_expanded");
  // tested for expansion = 2 case.
  std::map<Coordinate, std::vector<M> > send_map;
  std::map<Coordinate, int> send_map_consume;

  Coordinate pos; // coordinate position of a site relative to this node
  Coordinate local_pos; // coordinate position of a site relative to its home node
  Coordinate node_pos; // home node coordinate of a site in node space

  // populate send_map with the data that we need to send to other nodes
  long record_size = field_comm.geo.local_volume_expanded();
  for(long record = 0; record < record_size; record++){
    pos = field_comm.geo.coordinateFromRecord(record);
    if(field_comm.geo.is_local(pos)) continue;
    for(int mu = 0; mu < DIMN; mu++){
      local_pos[mu] = pos[mu] % field_comm.geo.node_site[mu];
      node_pos[mu] = pos[mu] / field_comm.geo.node_site[mu];
      if(local_pos[mu] < 0){
        local_pos[mu] += field_comm.geo.node_site[mu];
        node_pos[mu]--;
      }
    }
    std::vector<M> &vec = send_map[node_pos];
    for(int mu = 0; mu < field_comm.geo.multiplicity; mu++)
      vec.push_back(field_comm.get_elems_const(local_pos)[mu]);
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
    id_this = get_id_node();
    coor_this = qlat::coordinate_from_index(id_this, \
        field_comm.geo.geon.size_node);
    coort = coor_this - node_pos; 
    regularize_coordinate(coort, field_comm.geo.geon.size_node);
    coorf = coor_this + node_pos;
    regularize_coordinate(coorf, field_comm.geo.geon.size_node);

    idt = qlat::index_from_coordinate(coort, field_comm.geo.geon.size_node);
    idf = qlat::index_from_coordinate(coorf, field_comm.geo.geon.size_node);

    MPI_Request req;
    MPI_Isend((void*)send, size_bytes, MPI_BYTE, idt, 0, get_comm(), &req);
    const int ret = MPI_Recv((void*)recv, size_bytes, MPI_BYTE, \
        idf, 0, get_comm(), MPI_STATUS_IGNORE);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    qassert(!ret);

    memcpy(send, recv, size_bytes);

    send_map_consume[node_pos] = 0;
  }
  // Now send_map[node_pos] is the vector of data recieved from the node
  // pointed to by key.
  for(long record = 0; record < record_size; record++){
    pos = field_comm.geo.coordinateFromRecord(record);
    if(field_comm.geo.is_local(pos)) continue;
    for(int mu = 0; mu < DIMN; mu++){
      local_pos[mu] = pos[mu] % field_comm.geo.node_site[mu];
      node_pos[mu] = pos[mu] / field_comm.geo.node_site[mu];
      if(local_pos[mu] < 0){
        local_pos[mu] += field_comm.geo.node_site[mu];
        node_pos[mu]--;
      }
    }
    // send_map_consume[key] keeps track of our progress in consuming the
    // received data in sendmap[key], so that we know which offset of
    // send_map[node_pos] corresponds to which site.
    int consume = send_map_consume[node_pos];
    std::vector<M> &vec = send_map[node_pos];
    for(int mu = 0; mu < field_comm.geo.multiplicity; mu++){
      field_comm.get_elems(pos)[mu] = vec[consume];
      consume++;
    }
    send_map_consume[node_pos] = consume;
  }
}

QLAT_END_NAMESPACE
