/* Copyright (c) 2013-2017. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "simgrid/s4u/Engine.hpp"
#include "src/instr/instr_private.h"
#include "src/plugins/vm/VirtualMachineImpl.hpp"

XBT_LOG_EXTERNAL_DEFAULT_CATEGORY(surf_kernel);

/*********
 * TOOLS *
 *********/

extern double NOW;

void surf_presolve()
{
  double next_event_date = -1.0;
  tmgr_trace_event_t event          = nullptr;
  double value = -1.0;
  simgrid::surf::Resource *resource = nullptr;

  XBT_DEBUG ("Consume all trace events occurring before the starting time.");
  while ((next_event_date = future_evt_set->next_date()) != -1.0) {
    if (next_event_date > NOW)
      break;

    while ((event = future_evt_set->pop_leq(next_event_date, &value, &resource))) {
      if (value >= 0)
        resource->apply_event(event, value);
    }
  }

  XBT_DEBUG ("Set every models in the right state by updating them to 0.");
  for (auto model : *all_existing_models)
      model->updateActionsState(NOW, 0.0);
}

double surf_solve(double max_date)
{
  double time_delta = -1.0; /* duration */
  double next_event_date = -1.0;
  double model_next_action_end = -1.0;
  double value = -1.0;
  simgrid::surf::Resource *resource = nullptr;
  tmgr_trace_event_t event          = nullptr;

  if (max_date > 0.0) {
    xbt_assert(max_date > NOW,"You asked to simulate up to %f, but that's in the past already", max_date);

    time_delta = max_date - NOW;
  }

  /* Physical models MUST be resolved first */
  XBT_DEBUG("Looking for next event in physical models");
  double next_event_phy = surf_host_model->nextOccuringEvent(NOW);
  if ((time_delta < 0.0 || next_event_phy < time_delta) && next_event_phy >= 0.0) {
    time_delta = next_event_phy;
  }
  if (surf_vm_model != nullptr) {
    XBT_DEBUG("Looking for next event in virtual models");
    double next_event_virt = surf_vm_model->nextOccuringEvent(NOW);
    if ((time_delta < 0.0 || next_event_virt < time_delta) && next_event_virt >= 0.0)
      time_delta = next_event_virt;
  }

  XBT_DEBUG("Min for resources (remember that NS3 don't update that value): %f", time_delta);

  XBT_DEBUG("Looking for next trace event");

  while (1) { // Handle next occurring events until none remains
    next_event_date = future_evt_set->next_date();
    XBT_DEBUG("Next TRACE event: %f", next_event_date);

    if (not surf_network_model->nextOccuringEventIsIdempotent()) { // NS3, I see you
      if (next_event_date!=-1.0 && time_delta!=-1.0) {
        time_delta = MIN(next_event_date - NOW, time_delta);
      } else {
        time_delta = MAX(next_event_date - NOW, time_delta); // Get the positive component
      }

      XBT_DEBUG("Run the NS3 network at most %fs", time_delta);
      // run until min or next flow
      model_next_action_end = surf_network_model->nextOccuringEvent(time_delta);

      XBT_DEBUG("Min for network : %f", model_next_action_end);
      if(model_next_action_end>=0.0)
        time_delta = model_next_action_end;
    }

    if (next_event_date < 0.0) {
      XBT_DEBUG("no next TRACE event. Stop searching for it");
      break;
    }

    if ((time_delta == -1.0) || (next_event_date > NOW + time_delta))
      break; // next event occurs after the next resource change, bail out

    XBT_DEBUG("Updating models (min = %g, NOW = %g, next_event_date = %g)", time_delta, NOW, next_event_date);

    while ((event = future_evt_set->pop_leq(next_event_date, &value, &resource))) {
      if (resource->isUsed() || xbt_dict_get_or_null(watched_hosts_lib, resource->cname())) {
        time_delta = next_event_date - NOW;
        XBT_DEBUG("This event invalidates the next_occuring_event() computation of models. Next event set to %f", time_delta);
      }
      // FIXME: I'm too lame to update NOW live, so I change it and restore it so that the real update with surf_min will work
      double round_start = NOW;
      NOW = next_event_date;
      /* update state of the corresponding resource to the new value. Does not touch lmm.
         It will be modified if needed when updating actions */
      XBT_DEBUG("Calling update_resource_state for resource %s", resource->cname());
      resource->apply_event(event, value);
      NOW = round_start;
    }
  }

  /* FIXME: Moved this test to here to avoid stopping simulation if there are actions running on cpus and all cpus are with availability = 0.
   * This may cause an infinite loop if one cpu has a trace with periodicity = 0 and the other a trace with periodicity > 0.
   * The options are: all traces with same periodicity(0 or >0) or we need to change the way how the events are managed */
  if (time_delta < 0) {
    XBT_DEBUG("No next event at all. Bail out now.");
    return -1.0;
  }

  XBT_DEBUG("Duration set to %f", time_delta);

  // Bump the time: jump into the future
  NOW = NOW + time_delta;

  // Inform the models of the date change
  for (auto model : *all_existing_models) {
    model->updateActionsState(NOW, time_delta);
  }
  simgrid::s4u::onTimeAdvance(time_delta);

  TRACE_paje_dump_buffer (0);

  return time_delta;
}

/*********
 * MODEL *
 *********/
static surf_action_t ActionListExtract(simgrid::surf::ActionList* list)
{
  if (list->empty())
    return nullptr;
  surf_action_t res = &list->front();
  list->pop_front();
  return res;
}

surf_action_t surf_model_extract_done_action_set(surf_model_t model)
{
  return ActionListExtract(model->getDoneActionSet());
}

surf_action_t surf_model_extract_failed_action_set(surf_model_t model){
  return ActionListExtract(model->getFailedActionSet());
}

int surf_model_running_action_set_size(surf_model_t model){
  return model->getRunningActionSet()->size();
}

surf_action_t surf_host_open(sg_host_t host, const char* fullpath){
  return host->pimpl_->open(fullpath);
}

surf_action_t surf_host_close(sg_host_t host, surf_file_t fd){
  return host->pimpl_->close(fd);
}

int surf_host_unlink(sg_host_t host, surf_file_t fd){
  return host->pimpl_->unlink(fd);
}

size_t surf_host_get_size(sg_host_t host, surf_file_t fd){
  return host->pimpl_->getSize(fd);
}

surf_action_t surf_host_read(sg_host_t host, surf_file_t fd, sg_size_t size){
  return host->pimpl_->read(fd, size);
}

surf_action_t surf_host_write(sg_host_t host, surf_file_t fd, sg_size_t size){
  return host->pimpl_->write(fd, size);
}

xbt_dynar_t surf_host_get_info(sg_host_t host, surf_file_t fd){
  return host->pimpl_->getInfo(fd);
}

size_t surf_host_file_tell(sg_host_t host, surf_file_t fd){
  return host->pimpl_->fileTell(fd);
}

int surf_host_file_seek(sg_host_t host, surf_file_t fd,
                               sg_offset_t offset, int origin){
  return host->pimpl_->fileSeek(fd, offset, origin);
}

int surf_host_file_move(sg_host_t host, surf_file_t fd, const char* fullpath){
  return host->pimpl_->fileMove(fd, fullpath);
}

sg_size_t surf_storage_get_size(surf_storage_t resource)
{
  return static_cast<simgrid::surf::StorageImpl*>(resource)->size_;
}

sg_size_t surf_storage_get_free_size(surf_storage_t resource)
{
  return static_cast<simgrid::surf::StorageImpl*>(resource)->getFreeSize();
}

sg_size_t surf_storage_get_used_size(surf_storage_t resource)
{
  return static_cast<simgrid::surf::StorageImpl*>(resource)->getUsedSize();
}

xbt_dict_t surf_storage_get_properties(surf_storage_t resource)
{
  return static_cast<simgrid::surf::StorageImpl*>(resource)->getProperties();
}

const char* surf_storage_get_host(surf_storage_t resource)
{
  return static_cast<simgrid::surf::StorageImpl*>(resource)->attach_;
}

const char* surf_storage_get_name(surf_storage_t resource)
{
  return static_cast<simgrid::surf::StorageImpl*>(resource)->cname();
}

void surf_cpu_action_set_bound(surf_action_t action, double bound) {
  static_cast<simgrid::surf::CpuAction*>(action)->setBound(bound);
}

surf_file_t surf_storage_action_get_file(surf_action_t action){
  return static_cast<simgrid::surf::StorageAction*>(action)->file_;
}
