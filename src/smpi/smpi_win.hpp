/* Copyright (c) 2010, 2013-2015. The SimGrid Team.
 * All rights reserved.                                                     */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SMPI_WIN_HPP_INCLUDED
#define SMPI_WIN_HPP_INCLUDED

#include "private.h"
#include <vector>

namespace simgrid{
namespace smpi{

class Win {
  private :
  void* base_;
  MPI_Aint size_;
  int disp_unit_;
  int assert_;
  MPI_Info info_;
  MPI_Comm comm_;
  std::vector<MPI_Request> *requests_;
  xbt_mutex_t mut_;
  msg_bar_t bar_;
  MPI_Win* connected_wins_;
  char* name_;
  int opened_;
  MPI_Group group_;
  int count_; //for ordering the accs

public:
  Win(void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm);
  ~Win();
  void get_name( char* name, int* length);
  void get_group( MPI_Group* group);
  void set_name( char* name);
  int start(MPI_Group group, int assert);
  int post(MPI_Group group, int assert);
  int complete();
  int wait();
  int fence(int assert);
  int put( void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank,
              MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype);
  int get( void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank,
              MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype);
  int accumulate( void *origin_addr, int origin_count, MPI_Datatype origin_datatype, int target_rank,
              MPI_Aint target_disp, int target_count, MPI_Datatype target_datatype, MPI_Op op);
};


}
}

#endif