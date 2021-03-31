/**
 * @file MPIWrapper.cpp
 * @author Oleksandr Ostrenko <oleksandr.ostrenko@tu-dresden.de>
 * @author Pietro Incardona <incardon@mpi-cbg.de>
 * @author Rajesh Ramaswamy <rrajesh@pks.mpg.de>
 * @version 1.0.0
 * @date Mon, 10 Feb 2017
 * @section LICENSE
 * 
 * The GNU LGPL v3 or any later version is applied to this software, see the LICENSE.txt file.
 * 
 * @section DESCRIPTION
 *
 * Implementation of the MPI interface for the \c SimulationInfo class.
 */

#include "../../include/stdheaders.h"
// #include "../../include/util/IO.hpp"
#include "../../include/util/MPIWrapper.h"
#include "../../include/datamodel/SimulationInfo.h"

#ifdef HAVE_MPI

#if __GNUC__ > 4 || \
    (__GNUC__ == 4 && (__GNUC_MINOR__ >= 2))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#endif

namespace pssalib
{
namespace mpi
{
  ////////////////////////////////
  // Global functions
  bool checkFilterMpiStdOut(void *)
  {
    return pssalib::getMPIWrapperInstance().isMaster();
  }

  bool checkFilterMpiStdErr(void *)
  {
    return pssalib::getMPIWrapperInstance().isMaster();
  }

  ////////////////////////////////
  // Enumerators

  typedef enum tagMPISpreadFlags {
    msfChunksCalculated  = 0x1,
    msfChunksDistributed = 0x2,
    msfChunksCompleted   = 0x4
  } MPISpreadFlags;

  ////////////////////////////////
  // Constructors

  //! Default constructor : initialize MPI
  MPIWrapper::MPIWrapper()
    : nRank(0)
    , nPoolSize(1)
  {
// std::cerr << "MPIWrapper::MPIWrapper()" << std::endl;
    int flag = 0;
    MPI_Initialized(&flag);
    if(!flag)
      MPI_Init(NULL,NULL);

    MPI_Comm_size(MPI_COMM_WORLD,&nPoolSize);
    MPI_Comm_rank(MPI_COMM_WORLD,&nRank);
  }

  // Default destructor : finalize MPI
  MPIWrapper::~MPIWrapper()
  {
    MPI_Finalize();
  }

  /**
   * Calculate the chunk size that is assigned to this process.
   * 
   * @param size Total size of the data.
   * @return Size of the chunk assigned to this process.
   */
  INTEGER MPIWrapper::getDataChunkSize(INTEGER size)
  {
    if(size > nPoolSize)
    {
      INTEGER block = (INTEGER)(size / nPoolSize);

      // Increase the block, for the ramaining of trajectories
      if (nRank < (size % nPoolSize))
        block++;

      return block;
    }
    else
    {
      if(nRank < size)
        return 1;
      return 0;
    }
  }

  ////////////////////////////////
  // Methods

  /**
   * Prepare the process context for a parallel spread.
   * 
   * @param ptrSimInfo Pointer to the @link SimulationInfo object associated with this run.
   * @param size Total size of data to be processed.
   * @return @true if context initialisation was successful, or @false otherwise.
   */
  bool MPIWrapper::pre_spread(datamodel::SimulationInfo * ptrSimInfo)
  {
    if(ptrSimInfo->unFlags & msfChunksDistributed)
    {
      PSSA_ERROR(ptrSimInfo, << "invalid state." << std::endl);
      return false;
    }

    // reset flags
    ptrSimInfo->unFlags &= ~(msfChunksCalculated | msfChunksDistributed | msfChunksCompleted);

    if(ptrSimInfo->unSamplesTotal <= nPoolSize)
    {
      ptrSimInfo->m_nBlockStart = nRank;

      if(ptrSimInfo->unSamplesTotal > nRank)
        ptrSimInfo->m_nBlockSize = 1;
      else
        ptrSimInfo->m_nBlockSize = 0;
    }
    else
    {
      ptrSimInfo->m_nBlockStart = 0;

      INTEGER block = (INTEGER)(ptrSimInfo->unSamplesTotal / nPoolSize),
              threshold = ptrSimInfo->unSamplesTotal % nPoolSize;
      for(INTEGER n = 0; n < nRank; n++)
      {
        ptrSimInfo->m_nBlockStart += block;
        if (n < threshold)
          ptrSimInfo->m_nBlockStart++;
      }

      ptrSimInfo->m_nBlockSize = block;
      if (nRank < threshold)
        ptrSimInfo->m_nBlockSize++;
    }

    ptrSimInfo->unFlags |= msfChunksCalculated;

    PSSA_INFO(ptrSimInfo, << "BlockStart = " << ptrSimInfo->m_nBlockStart
      << "; BlockSize = " << ptrSimInfo->m_nBlockSize << std::endl);

    return true;
  }

  /**
   * Allocate a block matching the current process context for a parallel spread.
   * 
   * @param ptrSimInfo Pointer to the @link SimulationInfo object associated with this run.
   * @param sizeType Data type size as returned by a call to <code>sizeof(<data_type>)</code>.
   */
  void * MPIWrapper::spread_alloc(datamodel::SimulationInfo * ptrSimInfo, size_t sizeType) const
  {
    if((0 == (ptrSimInfo->unFlags & msfChunksCalculated))||
       (0 == ptrSimInfo->m_nBlockSize))
      return NULL;

    return new char[sizeType*ptrSimInfo->m_nBlockSize];
  }

  /**
   * Spread in a parallel for fashion.
   * 
   * @param counter circular counter (cycle length is set by \b MPIWrapper::pre_spread()) 
   * @return @true if there is still work to do, or @false otherwise
   */
  bool MPIWrapper::spread(datamodel::SimulationInfo * ptrSimInfo, UINTEGER &counter)
  {
    if((0 == (ptrSimInfo->unFlags & msfChunksCalculated))||
       (ptrSimInfo->unFlags & msfChunksCompleted))
    {
      PSSA_ERROR(ptrSimInfo, << "invalid state." << std::endl);
      return false;
    }

//     if(0 == ptrSimInfo->m_nBlockSize)
//     {
//       PSSA_INFO(ptrSimInfo, << "skipping par-for." << std::endl);
//       ptrSimInfo->unFlags &= ~(msfChunksDistributed);
//       return false;
//     }

    if (0 == (ptrSimInfo->unFlags & msfChunksDistributed))
    {
      counter = ptrSimInfo->m_nBlockStart;
      ptrSimInfo->unFlags |= msfChunksDistributed;
      PSSA_INFO(ptrSimInfo, << "starting par-for." << std::endl);
    }
    else
      ++counter;

    if (counter >= (ptrSimInfo->m_nBlockStart + ptrSimInfo->m_nBlockSize))
    {
      PSSA_INFO(ptrSimInfo, << "terminating par-for." << std::endl);
      ptrSimInfo->unFlags &= ~(msfChunksDistributed);
      return false;
    }
    else
      return true;
  }

  /**
   * Query whether this process has to wait for other processes due to
   * uneven spreading. Useful if inter=process communication takes place
   * within the parallel for loop.
   * 
   * @param ptrSimInfo Pointer to the @link SimulationInfo object associated with this run.
   * @return @true if the process has to wait for other in the pool, @false othewise.
   */
  bool MPIWrapper::post_spread(datamodel::SimulationInfo * ptrSimInfo)
  {
    if((0 == (ptrSimInfo->unFlags & msfChunksCalculated))||
      (ptrSimInfo->unFlags & msfChunksCompleted))
    {
      PSSA_ERROR(ptrSimInfo, << "invalid state." << std::endl);
      return false;
    }

    // check if called before finishing the chunk
    if(ptrSimInfo->unFlags & msfChunksDistributed)
    {
      PSSA_INFO(ptrSimInfo, << "premature call." << std::endl);
      ptrSimInfo->unFlags &= ~(msfChunksDistributed);
      return false;
    }

    ptrSimInfo->unFlags |= msfChunksCompleted;

//     UINTEGER size = 0;
//     if(MPI_SUCCESS != MPI_Allreduce(&ptrSimInfo->m_nBlockSize, &size, sizeof(UINTEGER), MPI_CHAR, MPI_SUM, MPI_COMM_WORLD))
//       return false;
//     PSSA_INFO(ptrSimInfo, << "TotalSize = " << size << std::endl);

    if(ptrSimInfo->unSamplesTotal <= nPoolSize)
    {
      return (nRank >= ptrSimInfo->unSamplesTotal);
    }
    else
    {
      return (nRank >= (ptrSimInfo->unSamplesTotal % nPoolSize));
    }
  }

  /**
   * Collect results from all processes into a buffer owned by the master process.
   * Since block size can be different, MPI_Gather was not suitable for this.
   *
   * @param sbuf Buffer to be sent to the master process
   * @param rbuf Pointer to a pointer of the receiving buffer (the pointer itself must be NULL) [OUT]
   * @param sizeType Data type size as returned by a call to <code>sizeof(<data_type>)</code>.
   */
  bool MPIWrapper::spread_collect(datamodel::SimulationInfo * ptrSimInfo, void * sbuf, void ** rbuf, size_t sizeType)
  {
    if(((sbuf == NULL)&&(0 != ptrSimInfo->m_nBlockSize))||
      (0 == (ptrSimInfo->unFlags & msfChunksCalculated)))
      return false;

    if(0 == nRank)
    {
      bool bAllOK = true;

      if((rbuf == NULL)||(NULL != *rbuf))
        bAllOK = false;

      if(!sync_results(bAllOK))
        return false;

      INTEGER recvSize = 0;
      boost::scoped_array<INTEGER> arSize;

      try
      {
        arSize.reset(new INTEGER[nPoolSize]);
      }
      catch(std::bad_alloc& e)
      {
        PSSA_ERROR(ptrSimInfo, << e.what() << ": Unable to allocate memory." << std::endl);
        bAllOK = false;
      }

      if(!sync_results(bAllOK))
        return false;

      if(MPI_SUCCESS != MPI_Gather(&ptrSimInfo->m_nBlockSize, sizeof(INTEGER), MPI_CHAR,
         arSize.get(), sizeof(INTEGER), MPI_CHAR, 0, MPI_COMM_WORLD))
        return false;

      for(INTEGER i = 0; i < nPoolSize; i++)
        recvSize += arSize[i];
// PSSA_INFO(ptrSimInfo, << "TotalSize = " << recvSize << std::endl);

      // receving buffer
      char * tbuf = new char[sizeType*recvSize];

      // copy master data
      memcpy(tbuf, sbuf, ptrSimInfo->m_nBlockSize*sizeType);

      // Receive all chunks
      for (UINTEGER i = 1, offs = ptrSimInfo->m_nBlockSize*sizeType;
           i < nPoolSize; offs += arSize[i]*sizeType, i++)
      {
        if(arSize[i] > 0)
          MPI_Recv(&((char *)tbuf)[offs],arSize[i]*sizeType,MPI_CHAR,i,0,MPI_COMM_WORLD,NULL);
//PSSA_INFO(ptrSimInfo, << "Receive buffer size after communication with proc #" << i << " [" << recvSize << "] :" << std::endl);
      }

      // store result
      *rbuf = tbuf;
    }
    else
    {
      // check input arguments
      if(!sync_results(true))
        return false;

      // dynaic allocation status
      if(!sync_results(true))
        return false;

      // Report block sizes to master
      if(MPI_SUCCESS != MPI_Gather(&ptrSimInfo->m_nBlockSize, sizeof(INTEGER), MPI_CHAR,
         NULL, sizeof(INTEGER), MPI_CHAR, 0, MPI_COMM_WORLD))
        return false;

//       PSSA_INFO(ptrSimInfo, << "Send buffer size before communication with master [" << ptrSimInfo->m_nBlockSize << "] :" << std::endl);
      // Send the chunk
      if(ptrSimInfo->m_nBlockSize > 0)
        MPI_Send((char *)sbuf,ptrSimInfo->m_nBlockSize*sizeType,MPI_CHAR,0,0,MPI_COMM_WORLD);
    }

//     PSSA_INFO(ptrSimInfo, << "spread_collect - all OK" << std::endl);
    return true;
  }

  /**
   * Sync the result of an operation in the process pool.
   * 
   * @param result Result of the operation
   * @return @true if the operation succeeded on all processors, @false otherwise.
   */
  bool MPIWrapper::sync_results(bool result) const
  {
    char in = (char)result;
    char sync = 0;

    if(MPI_SUCCESS == MPI_Allreduce(&in, &sync, 1, MPI_CHAR, MPI_BAND, MPI_COMM_WORLD))
      return (sync != 0);
    return false;
  }

  /**
   * Broadcast a buffer from the one of the processes (by default the master process)
   * to the rest of the process pool.
   * 
   * @param sbuf Buffer to be sent to the master process
   * @param size Buffer size as returned by a call to <code>sizeof(<buffer_type>)</code>
   */
  bool MPIWrapper::broadcast(void * buf, int size, int source) const
  {
    if ((NULL == buf) || (0 == size)) return false;

    return (MPI_SUCCESS == MPI_Bcast(buf, size, MPI_CHAR, source, MPI_COMM_WORLD)); 
  }

  /**
   * Perform a reduce opertaion on all send buffers and transmit the result to the process pool.
   * 
   * @param sbuf Buffer to be sent
   * @param rbuf Buffer to receive result
   * @param size Buffer size as returned by a call to <code>sizeof(<buffer_type>)</code>
   * @param op Reduce operation to be perfromed (see @link MPI_Allreduce)
   */
  bool MPIWrapper::allreduce(void * sbuf, void * rbuf, int size, MPI_Op op) const
  {
    if ((NULL == sbuf) || (NULL == rbuf) || (0 == size)) return false;

    return (MPI_SUCCESS == MPI_Allreduce(sbuf, rbuf, size, MPI_CHAR, op, MPI_COMM_WORLD)); 
  }

  /**
   * Get the size of sent buffer..
   * 
   * @param source Process that sent the buffer
   * @param tag Tag assigned to the buffer
   * @return size as returned by a call to <code>sizeof(<buffer_type>)</code>
   *         or -1 if the respective calls failed.
   */
  int MPIWrapper::getSentBufSize(int source, int tag)
  {
    int msglen = -1;
    static MPI_Status stat;

    if(MPI_SUCCESS == MPI_Probe(source, tag, MPI_COMM_WORLD, &stat))
      if(MPI_SUCCESS != MPI_Get_count(&stat, MPI_CHAR, &msglen))
        msglen = -1;

    return msglen;
  }

  /**
   * Initialize output streams
   */
  void MPIWrapper::setupOutput()
  {
    // Null stream buffer
    static pssalib::io::null_streambuf< STRING::value_type > nullBuffer;
    // Null stream
    static std::ostream        nullStream(&nullBuffer);

    fosMpiCout.reset();
    fosMpiCout.push(pssalib::io::onoff_tee_filter<std::ostream>
      (std::cout,&checkFilterMpiStdOut,NULL));
    fosMpiCout.push(nullStream);

    fosMpiCerr.reset();
    fosMpiCerr.push(pssalib::io::onoff_tee_filter<std::ostream>
      (std::cerr,&checkFilterMpiStdErr,NULL));
    fosMpiCerr.push(nullStream);
  }
} // close mpi namespace

  /**
   * Get the global MPI wrapper.
   * @return Global MPIWrapper object.
   */
  mpi::MPIWrapper & getMPIWrapperInstance()
  {
//     std::cout << "getMPIWrapperInstance()" << std::endl;
    static mpi::MPIWrapper mpiWrapperInstance;
    return mpiWrapperInstance;
  }

} // close pssalib namespace

#if __GNUC__ > 4 || \
    (__GNUC__ == 4 && (__GNUC_MINOR__ >= 2))
#pragma GCC diagnostic pop
#endif

#endif /* HAVE_MPI */
