/* Copyright (c) 2016. The SimGrid Team. All rights reserved.               */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#ifndef SIMGRID_ROUTING_NETZONEIMPL_HPP
#define SIMGRID_ROUTING_NETZONEIMPL_HPP

#include <map>

#include "xbt/graph.h"

#include "simgrid/s4u/NetZone.hpp"
#include "simgrid/s4u/forward.hpp"

namespace simgrid {
namespace kernel {
class EngineImpl;
namespace routing {
class BypassRoute;

/** @ingroup ROUTING_API
 *  @brief Private implementation of the Networking Zones
 *
 * A netzone is a network container, in charge of routing information between elements (hosts and sub-netzones)
 * and to the nearby netzones. In SimGrid, there is a hierarchy of netzones, ie a tree with a unique root
 * NetZone, that you can retrieve with simgrid::s4u::Engine::netRoot().
 *
 * The purpose of the kernel::routing module is to retrieve the routing path between two points in a time- and
 * space-efficient manner. This is done by NetZoneImpl::getGlobalRoute(), called when creating a communication to
 * retrieve both the list of links that the create communication will use, and the summed latency that these
 * links represent.
 *
 * The network model could recompute the latency by itself from the list, but it would require an additional
 * traversal of the link set. This operation being on the critical path of SimGrid, the routing computes the
 * latency on the behalf of the network while constructing the link set.
 *
 * Finding the path between two nodes is rather complex because we navigate a hierarchy of netzones, each of them
 * being a full network. In addition, the routing can declare shortcuts (called bypasses), either within a NetZone
 * at the route level or directly between NetZones. Also, each NetZone can use a differing routing algorithm, depending
 * on its class. @ref FullZone have a full matrix giving explicitly the path between any pair of their
 * contained nodes, while @ref DijkstraZone or @ref FloydZone rely on a shortest path algorithm. @ref VivaldiZone
 * does not even have any link but only use only coordinate information to compute the latency.
 *
 * So NetZoneImpl::getGlobalRoute builds the path recursively asking its specific information to each traversed NetZone
 * with NetZoneImpl::getLocalRoute, that is redefined in each sub-class.
 * The algorithm for that is explained in http://hal.inria.fr/hal-00650233/ (but for historical reasons, NetZones are
 * called Autonomous Systems in this article).
 *
 */
XBT_PUBLIC_CLASS NetZoneImpl : public s4u::NetZone
{
  friend simgrid::kernel::EngineImpl; // it destroys netRoot_

protected:
  explicit NetZoneImpl(NetZone * father, const char* name);
  virtual ~NetZoneImpl();

public:
  /** @brief Make an host within that NetZone */
  simgrid::s4u::Host* createHost(const char* name, std::vector<double>* speedPerPstate, int coreAmount);
  /** @brief Creates a new route in this NetZone */
  void addBypassRoute(sg_platf_route_cbarg_t e_route) override;

protected:
  /**
   * @brief Probe the routing path between two points that are local to the called NetZone.
   *
   * @param src where from
   * @param dst where to
   * @param into Container into which the traversed links and gateway informations should be pushed
   * @param latency Accumulator in which the latencies should be added (caller must set it to 0)
   */
  virtual void getLocalRoute(NetCard * src, NetCard * dst, sg_platf_route_cbarg_t into, double* latency) = 0;
  /** @brief retrieves the list of all routes of size 1 (of type src x dst x Link) */
  /* returns whether we found a bypass path */
  bool getBypassRoute(routing::NetCard * src, routing::NetCard * dst,
                      /* OUT */ std::vector<surf::Link*> * links, double* latency);

public:
  /* @brief get the route between two nodes in the full platform
   *
   * @param src where from
   * @param dst where to
   * @param links Accumulator in which all traversed links should be pushed (caller must empty it)
   * @param latency Accumulator in which the latencies should be added (caller must set it to 0)
   */
  static void getGlobalRoute(routing::NetCard * src, routing::NetCard * dst,
                             /* OUT */ std::vector<surf::Link*> * links, double* latency);

  virtual void getGraph(xbt_graph_t graph, xbt_dict_t nodes, xbt_dict_t edges) = 0;
  enum class RoutingMode {
    unset = 0, /**< Undefined type                                   */
    base,      /**< Base case: use simple link lists for routing     */
    recursive  /**< Recursive case: also return gateway information  */
  };
  /* FIXME: protect the following fields once the construction madness is sorted out */
  RoutingMode hierarchy_ = RoutingMode::unset;

private:
  std::map<std::pair<NetCard*, NetCard*>, BypassRoute*> bypassRoutes_; // src x dst -> route
  routing::NetCard* netcard_ = nullptr;                                // Our representative in the father NetZone
};
}
}
}; // Namespace simgrid::kernel::routing

#endif /* SIMGRID_ROUTING_NETZONEIMPL_HPP */