/**
 * @file MPIWrapper.h
 * @author Oleksandr Ostrenko <oleksandr.ostrenko@tu-dresden.de>
 * @author Pietro Incardona <incardon@mpi-cbg.de>
 * @author Rajesh Ramaswamy <rrajesh@pks.mpg.de>
 * @version 1.0.0
 * @date Mon, 10 Feb 2017
 * @section LICENSE
 * 
 * The GPLv2 or any later version is applied to this software, see the LICENSE.txt file.
 * 
 * @section DESCRIPTION
 *
 * Small class to spread simulation over multiple indipendent computing processes and collect result.
 */

#include "../typedefs.h"

#ifdef HAVE_MPI

#ifndef PSSALIB_UTIL_MPI_WRAPPER_H_
#define PSSALIB_UTIL_MPI_WRAPPER_H_

#include <boost/iostreams/filtering_stream.hpp>

#include "IO.hpp"

// Forward declarations
namespace pssalib
{
  namespace datamodel
  {
    class SimulationInfo;
  }
}

namespace pssalib
{
namespace mpi
{
  ////////////////////////////////////////
  //! \brief Custom wrapper for MPI routines
  class MPIWrapper
  {
  ////////////////////////////////
  // Attributes
  protected:
    //! Current processor rank
    int  nRank;
    //! Number of processors
    int  nPoolSize;

  ////////////////////////////////
  // Attributes
  public:
    //! Standard output stream wrapper
    boost::iostreams::filtering_ostream fosMpiCout;

    //! Standard error stream wrapper
    boost::iostreams::filtering_ostream fosMpiCerr;

  ////////////////////////////////
  // Constructors
  public:

    /**
     * Default constructor that initialises MPI
     * and acquires processor's rank and pool size.
     */
    MPIWrapper();

    /**
     * Default destructor that finalises MPI.
     */
    ~MPIWrapper();

  //////////////////////////////
  // Methods
  public:
    /**
     * Get the current process's rank in the global process pool.
     *
     * @return Current process's rank in the pool.
     */
    inline int getRank() const { return nRank; };

    /**
     * Get the current process pool size.
     *
     * @return Current process pool size.
     */
    inline int getPoolSize() const { return nPoolSize; };

    /**
     * Check if current process is the master process.
     *
     * @return @true if current process is the master process, @false otherwise.
     */
    inline bool isMaster() const { return (0 == nRank); };

    /**
     * Calculate the chunk size that is assigned to this process.
     *
     * @param size Total size of the data.
     * @return Size of the chunk assigned to this process.
     */
    INTEGER getDataChunkSize(INTEGER size);

    /**
    * Prepare the process context for a parallel spread.
    *
    * @param ptrSimInfo Pointer to the @link SimulationInfo object associated with this run.
    * @param size Total size of data to be processed.
     * @return @true if context initialisation was successful, or @false otherwise.
    */
    bool pre_spread(pssalib::datamodel::SimulationInfo * ptrSimInfo);//, INTEGER size);

    /**
     * Spread in a parallel for fashion.
     *
     * @param counter circular counter (cycle length is set by \b MPIWrapper::pre_spread()) 
     * @return @true if there is still work to do, or @false otherwise
     */
    bool spread(pssalib::datamodel::SimulationInfo * ptrSimInfo, UINTEGER &counter);

    /**
     * Allocate a block matching the current process context for a parallel spread.
     *
     * @param ptrSimInfo Pointer to the @link SimulationInfo object associated with this run.
     * @param sizeType Data type size as returned by a call to <code>sizeof(<data_type>)</code>.
     */
    void * spread_alloc(datamodel::SimulationInfo * ptrSimInfo, size_t sizeType) const;

    /**
    * Query whether this process has to wait for other processes due to
    * uneven spreading. Useful if inter=process communication takes place
    * within the parallel for loop.
    *
    * @param ptrSimInfo Pointer to the @link SimulationInfo object associated with this run.
    * @return @true if the process has to wait for other in the pool, @false othewise.
    */
    bool post_spread(pssalib::datamodel::SimulationInfo * ptrSimInfo);

    /**
     * Collect results from all processes into a buffer owned by the master process.
     * Since block size can be different, MPI_Gather was not suitable for this.
     *
     * @param sbuf Buffer to be sent to the master process
     * @param rbuf Pointer to a pointer of the receiving buffer (the pointer itself must be NULL) [OUT]
     * @param sizeType Data type size as returned by a call to <code>sizeof(<data_type>)</code>.
     */
    bool spread_collect(pssalib::datamodel::SimulationInfo * ptrSimInfo, void * sbuf,
                        void ** rbuf, size_t sizeType);

    /**
     * Sync the result of an operation in the process pool.
     *
     * @param result Result of the operation
     * @return @true if the operation succeeded on all processors, @false otherwise.
     */
    bool sync_results(bool result) const;

    /**
     * Broadcast a buffer from the one of the processes (by default the master process)
     * to the rest of the process pool.
     *
     * @param sbuf Buffer to be sent to the master process
     * @param size Buffer size as returned by a call to <code>sizeof(<buffer_type>)</code>
     */
    bool broadcast(void * buf, int size, int source = 0) const;

    /**
     * Perform a reduce opertaion on all send buffers and transmit the result to the process pool.
     *
     * @param sbuf Buffer to be sent
     * @param rbuf Buffer to receive result
     * @param size Buffer size as returned by a call to <code>sizeof(<buffer_type>)</code>
     * @param op Reduce operation to be perfromed (see @link MPI_Allreduce)
     */
    bool allreduce(void * sbuf, void * rbuf, int size, MPI_Op op) const;

    /**
     * Get the size of sent buffer.
     *
     * @param source Process that sent the buffer
     * @param tag Tag assigned to the buffer
     * @return size as returned by a call to <code>sizeof(<buffer_type>)</code>
     *         or -1 if the respective calls fail.
     */
    int getSentBufSize(int source = MPI_ANY_SOURCE, int tag = MPI_ANY_TAG);

    /**
     * Initialize output streams
     */
    void setupOutput();
  };
} // close mpi namespace

  /**
   * Get the global MPI wrapper.
   *
   * @return Global MPIWrapper object.
   */
  mpi::MPIWrapper &getMPIWrapperInstance();

} // close pssalib namespace

#endif /* PSSALIB_UTIL_MPI_WRAPPER_H_ */

#ifndef PSSALIB_MPI_IO_INIT
#  define PSSALIB_MPI_IO_INIT pssalib::getMPIWrapperInstance().setupOutput();
#endif

#ifndef PSSALIB_MPI_RANK
#  define PSSALIB_MPI_RANK pssalib::getMPIWrapperInstance().getRank()
#endif

#ifndef PSSALIB_MPI_POOL_SIZE
#  define PSSALIB_MPI_POOL_SIZE pssalib::getMPIWrapperInstance().getPoolSize()
#endif

#ifndef PSSALIB_MPI_IS_MASTER
#  define PSSALIB_MPI_IS_MASTER pssalib::getMPIWrapperInstance().isMaster()
#endif

#ifndef PSSALIB_MPI_COUT_OR_NULL
#  define PSSALIB_MPI_COUT_OR_NULL pssalib::getMPIWrapperInstance().fosMpiCout
#endif

#ifndef PSSALIB_MPI_CERR_OR_NULL
#  define PSSALIB_MPI_CERR_OR_NULL pssalib::getMPIWrapperInstance().fosMpiCerr
#endif

#else

#ifndef PSSALIB_MPI_IO_INIT
#  define PSSALIB_MPI_IO_INIT
#endif

#ifndef PSSALIB_MPI_RANK
#  define PSSALIB_MPI_RANK ((int)0)
#endif

#ifndef PSSALIB_MPI_POOL_SIZE
#  define PSSALIB_MPI_POOL_SIZE ((int)1)
#endif

#ifndef PSSALIB_MPI_IS_MASTER
#  define PSSALIB_MPI_IS_MASTER (true)
#endif

#ifndef PSSALIB_MPI_COUT_OR_NULL
#  define PSSALIB_MPI_COUT_OR_NULL std::cout
#endif

#ifndef PSSALIB_MPI_CERR_OR_NULL
#  define PSSALIB_MPI_CERR_OR_NULL std::cerr
#endif

#endif /* HAVE_MPI */
