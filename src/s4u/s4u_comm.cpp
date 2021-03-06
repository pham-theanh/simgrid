/* Copyright (c) 2006-2017. The SimGrid Team. All rights reserved.          */

/* This program is free software; you can redistribute it and/or modify it
 * under the terms of the license (GNU LGPL) which comes with this package. */

#include "xbt/log.h"
#include "src/msg/msg_private.h"

#include "simgrid/s4u/Comm.hpp"
#include "simgrid/s4u/Mailbox.hpp"

XBT_LOG_NEW_DEFAULT_SUBCATEGORY(s4u_comm,s4u_activity,"S4U asynchronous communications");

namespace simgrid {
namespace s4u {
Comm::~Comm()
{
  if (state_ == started && not detached_ && (pimpl_ == nullptr || pimpl_->state == SIMIX_RUNNING)) {
    XBT_INFO("Comm %p freed before its completion. Detached: %d, State: %d", this, detached_, state_);
    if (pimpl_ != nullptr)
      XBT_INFO("pimpl_->state: %d", pimpl_->state);
    else
      XBT_INFO("pimpl_ is null");
    xbt_backtrace_display_current();
  }
  if (pimpl_)
    pimpl_->unref();
}

s4u::CommPtr Comm::send_init(s4u::MailboxPtr chan)
{
  CommPtr res   = CommPtr(new s4u::Comm());
  res->sender_ = SIMIX_process_self();
  res->mailbox_ = chan;
  return res;
}

s4u::CommPtr Comm::recv_init(s4u::MailboxPtr chan)
{
  CommPtr res    = CommPtr(new s4u::Comm());
  res->receiver_ = SIMIX_process_self();
  res->mailbox_ = chan;
  return res;
}

void Comm::setRate(double rate) {
  xbt_assert(state_==inited);
  rate_ = rate;
}

void Comm::setSrcData(void * buff) {
  xbt_assert(state_==inited);
  xbt_assert(dstBuff_ == nullptr, "Cannot set the src and dst buffers at the same time");
  srcBuff_ = buff;
}
void Comm::setSrcDataSize(size_t size){
  xbt_assert(state_==inited);
  srcBuffSize_ = size;
}
void Comm::setSrcData(void * buff, size_t size) {
  xbt_assert(state_==inited);

  xbt_assert(dstBuff_ == nullptr, "Cannot set the src and dst buffers at the same time");
  srcBuff_ = buff;
  srcBuffSize_ = size;
}
void Comm::setDstData(void ** buff) {
  xbt_assert(state_==inited);
  xbt_assert(srcBuff_ == nullptr, "Cannot set the src and dst buffers at the same time");
  dstBuff_ = buff;
}
size_t Comm::getDstDataSize(){
  xbt_assert(state_==finished);
  return dstBuffSize_;
}
void Comm::setDstData(void ** buff, size_t size) {
  xbt_assert(state_==inited);

  xbt_assert(srcBuff_ == nullptr, "Cannot set the src and dst buffers at the same time");
  dstBuff_ = buff;
  dstBuffSize_ = size;
}

void Comm::start() {
  xbt_assert(state_ == inited);

  if (srcBuff_ != nullptr) { // Sender side
    pimpl_ = simcall_comm_isend(sender_, mailbox_->getImpl(), remains_, rate_,
        srcBuff_, srcBuffSize_,
        matchFunction_, cleanFunction_, copyDataFunction_,
        userData_, detached_);
  } else if (dstBuff_ != nullptr) { // Receiver side
    pimpl_ = simcall_comm_irecv(receiver_, mailbox_->getImpl(), dstBuff_, &dstBuffSize_,
        matchFunction_, copyDataFunction_,
        userData_, rate_);

  } else {
    xbt_die("Cannot start a communication before specifying whether we are the sender or the receiver");
  }
  state_ = started;
}
void Comm::wait() {
  xbt_assert(state_ == started || state_ == inited);

  if (state_ == started)
    simcall_comm_wait(pimpl_, -1/*timeout*/);
  else { // state_ == inited. Save a simcall and do directly a blocking send/recv
    if (srcBuff_ != nullptr) {
      simcall_comm_send(sender_, mailbox_->getImpl(), remains_, rate_,
          srcBuff_, srcBuffSize_,
          matchFunction_, copyDataFunction_,
          userData_, -1 /*timeout*/);
    } else {
      simcall_comm_recv(receiver_, mailbox_->getImpl(), dstBuff_, &dstBuffSize_,
          matchFunction_, copyDataFunction_,
          userData_, -1/*timeout*/, rate_);
    }
  }
  state_ = finished;
  if (pimpl_)
    pimpl_->unref();
}

void Comm::wait(double timeout) {
  xbt_assert(state_ == started || state_ == inited);

  if (state_ == started) {
    simcall_comm_wait(pimpl_, timeout);
    state_ = finished;
    pimpl_->unref();
    return;
  }

  // It's not started yet. Do it in one simcall
  if (srcBuff_ != nullptr) {
    simcall_comm_send(sender_, mailbox_->getImpl(), remains_, rate_,
        srcBuff_, srcBuffSize_,
        matchFunction_, copyDataFunction_,
        userData_, timeout);
  } else { // Receiver
    simcall_comm_recv(receiver_, mailbox_->getImpl(), dstBuff_, &dstBuffSize_,
        matchFunction_, copyDataFunction_,
        userData_, timeout, rate_);
  }
  state_ = finished;
  if (pimpl_)
    pimpl_->unref();
}

void Comm::send_detached(MailboxPtr dest, void* data, int simulatedSize)
{
  s4u::CommPtr res = CommPtr(s4u::Comm::send_init(dest));
  res->setRemains(simulatedSize);
  res->srcBuff_     = data;
  res->srcBuffSize_ = sizeof(void*);
  res->detached_    = true;
  res->start();
}

s4u::CommPtr Comm::send_async(MailboxPtr dest, void* data, int simulatedSize)
{
  s4u::CommPtr res = CommPtr(s4u::Comm::send_init(dest));
  res->setRemains(simulatedSize);
  res->srcBuff_     = data;
  res->srcBuffSize_ = sizeof(void*);
  res->start();
  return res;
}

s4u::CommPtr Comm::recv_async(MailboxPtr dest, void** data)
{
  s4u::CommPtr res = CommPtr(s4u::Comm::recv_init(dest));
  res->setDstData(data, sizeof(*data));
  res->start();
  return res;
}

void Comm::cancel()
{
  simgrid::kernel::activity::CommImpl* commPimpl = static_cast<simgrid::kernel::activity::CommImpl*>(pimpl_);
  commPimpl->cancel();
}

bool Comm::test() {
  xbt_assert(state_ == inited || state_ == started || state_ == finished);
  
  if (state_ == finished) 
    xbt_die("Don't call test on a finished comm.");
  
  if (state_ == inited) {
    this->start();
  }
  
  if(simcall_comm_test(pimpl_)){
    state_ = finished;
    pimpl_->unref();
    return true;
  }
  return false;
}

void intrusive_ptr_release(simgrid::s4u::Comm* c)
{
  if (c->refcount_.fetch_sub(1, std::memory_order_release) == 1) {
    std::atomic_thread_fence(std::memory_order_acquire);
    delete c;
  }
}
void intrusive_ptr_add_ref(simgrid::s4u::Comm* c)
{
  c->refcount_.fetch_add(1, std::memory_order_relaxed);
}
}
} // namespaces
