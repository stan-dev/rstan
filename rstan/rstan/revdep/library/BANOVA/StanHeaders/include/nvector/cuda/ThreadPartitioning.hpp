/*
 * -----------------------------------------------------------------
 * Programmer(s): Slaven Peles @ LLNL
 * -----------------------------------------------------------------
 * LLNS Copyright Start
 * Copyright (c) 2014, Lawrence Livermore National Security
 * This work was performed under the auspices of the U.S. Department
 * of Energy by Lawrence Livermore National Laboratory in part under
 * Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS Copyright End
 * -----------------------------------------------------------------
 */



#ifndef _THREAD_PARTITIONING_HPP_
#define _THREAD_PARTITIONING_HPP_

#include <iostream>
#include <cuda_runtime.h>

namespace suncudavec
{

template<class T, class I>
class StreamPartitioning
{
public:
  StreamPartitioning(I N, unsigned block)
: block_(block),
  grid_((N + block - 1) / block)
{
}

  explicit StreamPartitioning(StreamPartitioning<T, I>& p)
  : block_(p.block_),
    grid_(p.grid_)
  {
  }

  unsigned grid() const
  {
    return grid_;
  }

  unsigned block() const
  {
    return block_;
  }

private:
  unsigned block_;
  unsigned grid_;
};


template<class T, class I=int>
class ReducePartitioning
{
public:
  ReducePartitioning(I N, unsigned block)
: block_(block),
  grid_((N + (block_ * 2 - 1)) / (block_ * 2)),
  shMemSize_(block_*sizeof(T))
{
    allocateBuffer();
}

  explicit ReducePartitioning(StreamPartitioning<T, I>& p)
  : block_(p.block_),
    grid_(p.grid_),
    shMemSize_(p.shMemSize_)
  {
    allocateBuffer();
  }

  ~ReducePartitioning()
  {
    cudaError_t err;
    if (bufferSize_ > 0)
      free(h_buffer_);
    if (bufferSize_ > 0)
    {
      err = cudaFree(d_buffer_);
      if(err != cudaSuccess)
        std::cerr << "Failed to free device vector (error code " << err << ")!\n";
    }
  }

  int setPartitioning(I N, unsigned& grid, unsigned& block, unsigned& shMemSize)
  {
    block = block_;
    grid  = (N + (block * 2 - 1)) / (block * 2);
    shMemSize = block * sizeof(T);

    return 0;
  }

  unsigned grid() const
  {
    return grid_;
  }

  unsigned block() const
  {
    return block_;
  }

  unsigned shmem() const
  {
    return shMemSize_;
  }

  unsigned int buffSize()
  {
    return bufferSize_;
  }

  T* devBuffer()
  {
    return d_buffer_;
  }

  const T* devBuffer() const
  {
    return d_buffer_;
  }

  T* hostBuffer()
  {
    return h_buffer_;
  }

  const T* hostBuffer() const
  {
    return h_buffer_;
  }

  void copyFromDevBuffer(unsigned int n) const
  {
    cudaError_t err = cudaMemcpy(h_buffer_, d_buffer_, n*sizeof(T), cudaMemcpyDeviceToHost);
    if(err != cudaSuccess)
      std::cerr << "Failed to copy vector from device to host (error code " << err << ")!\n";
  }

private:
  int allocateBuffer()
  {
    bufferSize_ = grid_ * sizeof(T);
    h_buffer_ = static_cast<T*>(malloc(bufferSize_));
    if(h_buffer_ == NULL)
      std::cerr << "Failed to allocate host vector!\n";

    cudaError_t err;
    err = cudaMalloc((void**) &d_buffer_, bufferSize_);
    if(err != cudaSuccess)
      std::cerr << "Failed to allocate device vector (error code " << err << ")!\n";

    return 0;
  }

private:
  unsigned block_;
  unsigned grid_;
  unsigned shMemSize_;
  T* d_buffer_;
  T* h_buffer_;
  unsigned bufferSize_;

};


} // namespace suncudavec

#endif // _THREAD_PARTITIONING_HPP_
