/* Copyright (c) 2010-2012, 2014-2015. The SimGrid Team.
 * All rights reserved.                                                     */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef JED_SIMGRID_PLATFORM_H_
#define JED_SIMGRID_PLATFORM_H_

#include "simgrid_config.h"
#include "simgrid/forward.h"
#include "xbt/dynar.h"
#include <unordered_map>
#include <vector>
#include <string>

#if HAVE_JEDULE

namespace simgrid {
namespace jedule{
XBT_PUBLIC_CLASS Container {
public:
  Container(std::string name);
  virtual ~Container();
private:
  int last_id;
  int is_lowest = 0;
public:
  std::string name;
  std::unordered_map<const char*, int> name2id;
  Container *parent = nullptr;
  std::vector<Container*> children;
  std::vector<sg_host_t> resource_list;
  void addChild(Container* child);
  void addResources(std::vector<sg_host_t> hosts);
  void createHierarchy(AS_t from_as);
  std::vector<int> getHierarchy();
  std::string getHierarchyAsString();
  void print(FILE *file);
  void printResources(FILE *file);
};

}
}
SG_BEGIN_DECL()
typedef simgrid::jedule::Container * jed_container_t;

/* FIXME: jedule_container should be objectified too */
/** selection of a subset of resources from the original set */
struct jed_res_subset {
  jed_container_t parent;
  int start_idx; // start idx in resource_list of container
  int nres;      // number of resources spanning starting at start_idx
};

typedef struct jed_res_subset s_jed_res_subset_t, *jed_res_subset_t;

/**
 * it is assumed that the host_names in the entire system are unique that means that we don't need parent references
 *
 * subset_list must be allocated
 * host_names is the list of host_names associated with an event
 */
void jed_simgrid_get_resource_selection_by_hosts(xbt_dynar_t subset_list, std::vector<sg_host_t>* host_list);

/*
  global:
      hash host_id -> container
  container:
      hash host_id -> jed_host_id
      list <- [ jed_host_ids ]
      list <- sort( list )
      list_chunks <- chunk( list )   -> [ 1, 3-5, 7-9 ]
*/

SG_END_DECL()

#endif

#endif /* JED_SIMGRID_PLATFORM_H_ */