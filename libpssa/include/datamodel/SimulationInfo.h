/**
 * @file SimulationInfo.h
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
 * Declares a datatype that is used to marshall simulation data 
 * between library modules and provide mean to output results
 */

#ifndef PSSALIB_DATAMODEL_SIMULATIONINFO_H_
#define PSSALIB_DATAMODEL_SIMULATIONINFO_H_

#include "./detail/VolumeDecomposition.hpp"
#include "./../PSSA.h"
#include "./../util/MPIWrapper.h"

namespace pssalib
{
namespace datamodel
{
  ////////////////////////////////////////
  //! @brief Simulation data marshalling
  class SimulationInfo
  {
  /////////////////////////////////
  // Friend classes
  public:
#ifdef HAVE_MPI
    friend class pssalib::mpi::MPIWrapper;
#endif
    friend class pssalib::PSSA;

  /////////////////////////////////
  // Constructors
  public:
    // Constructor
    SimulationInfo();

    // Copy constructor
    SimulationInfo(SimulationInfo & right);

    // Destructor
    ~SimulationInfo();

  /////////////////////////////////
  // Data Types
  public:
    // Data type for statistic output
    typedef enum tagSimulationFlags
    {
      sfDelayedUpdate = 0x08,//!<"trial" mode: Statistics (std. deviation & mean)
                             //!<"pdf" mode: Statistics (one-dimensional histogram)
      sfUPDConsuming  = 0x10,//!<"update populations after a consuming reaction?
      sfUPDProducts   = 0x20,//!<Analyse available simulation data only
      sfUPDAll        = 0x40 //!<Output flags are preset externally
    } SimulationFlags;

    typedef enum tagOutputFlags {
      ofNone=0x0000,

      // Log message mask
      ofError=0x0001,      //!<Error messages
      ofWarning=0x0002,    //!<Warning messages
      ofInfo=0x0004,       //!<Information messages
      ofTrace=0x0008,      //!<Debug messages
      ofMaskLog=0x000F,    //!<All logging messages

      // Output streams
      ofStatus=0x0010,     //!<Output to status file
      ofLog=0x0020,        //!<Output to log file
      ofTrajectory=0x0040, //!<Output trajectories
      ofRawTrajectory=0x0080, //!<Output raw trajectories
      ofFinalPops=0x0100,  //!<Output data for PDFs
      ofRawFinalPops=0x0200, //!<Output raw data for PDFs
      ofTimePoints=0x0400, //!<Output time points
      ofTiming=0x0800,     //!<Output timing
      ofSpeciesIDs=0x1000, //!<Output species ids
      ofMaskFile=0x1FF0,   //!<All file output flags

      ofMaskAll=0x1FFF,    //!<All output flags

      // extended output flags
      eofModuleGrouping=0x10000,
      eofModuleSampling=0x20000,
      eofModuleUpdate=0x40000,

      eofModuleAll=0x70000 //!<All modules
    } OutputFlags;

  /////////////////////////////////
  // Attributes
  protected:
    //! Pointer to simulation engine class.
    pssalib::PSSA         *m_ptrPSSA;

#ifdef __linux__
    struct timespec       m_trialStart, m_trialEnd; //!< Internal timing variables
#elif defined(__MACH__)
    UINT64                m_trialStart, m_trialEnd; //!< Internal timing variables
    REAL_EXT              m_tickRatio;              //!< Internal timing variables
#elif defined(_WIN32)
    LARGE_INTEGER         m_trialStart, m_trialEnd; //!< Internal timing variables
    REAL_EXT              m_tickRatio;              //!< Internal timing variables
#endif

    //! Samples [RESERVED]
    UINTEGER              m_unSampleCurrent;

    //! Dimensions [RESERVED]
    BYTE                  m_uDims;
    UINTEGER              *m_arunDims;

#ifdef HAVE_MPI
    //! Process block size [RESERVED]
    INTEGER               m_nBlockSize;
    //! Process block offset [RESERVED]
    INTEGER               m_nBlockStart;
#endif

    //! Internal & external output streams
    FILESTREAMBUFFER      *m_arPtrFileBuffers[7];
    STREAMBUFFER          *m_arPtrExternalBuffers[7];

    //! Array of species indices for output [RESERVED]
    std::vector<UINTEGER> m_arSpeciesIdx;

    //! Data structure describing the current model
    detail::Model         m_Model;

    //! @internal Pointer to the output string buffer
    STRING::pointer       m_ptrOutputLine;

    //! @internal Pointer to the current position in trajectory buffer
    UINTEGER              *m_ptrRawTrajectory,
                          *m_ptrCurrPopulation;

    //! @internal Current index in the output
    UINTEGER              m_unOutputIdx,
                          m_unOutputMax;

  /////////////////////////////////
  // Attributes
  public:

    //! Filename templates
    static const STRING  arFileNames[];

    //! Generic flags
    UINTEGER             unFlags;

    //! Output flags [IN]
    UINTEGER             unOutputFlags;

    //! Output path (directory) [IN MANDATORY]
    STRING               strOutput;

    //! Array of species to output for the probability density [IN OPTIONAL, default: "all species"]
    std::vector<STRING>  *pArSpeciesIds;

    //! Number of samples the simulation [IN MANDATORY]
    UINTEGER             unSamplesTotal;

    // Simulation timing
    REAL dTimeCheckpoint, //!<last output time [RESERVED]
         dTimeStart,      //!<initial output time [IN OPTIONAL, default = 0.0]
         dTimeStep,       //!<time interval between subsequent outputs [IN MANDATORY]
         dTimeEnd,        //!<end time of the simulation [IN MANDATORY]
         dTimeSimulation; //!<current simulation time [RESERVED]

    //! Boundary description, see tagBoundaryConditionsType [IN MANDATORY]
    detail::BoundaryConditionsType eBoundaryConditions;

    //! Initial Population Type see tagInitialPopulationType [IN MANDATORY]
    detail::InitialPopulationType  eInitialPopulation;

    //! Population initializer [IN OPTIONAL, dafault: NULL]
    FCN_POPULATION_INITIALIZER     ptrPopulationInitializer;

    //! Pointer to population initializer user data [IN OPTIONAL, dafault: NULL]
    void                           *ptrPopulationInitializerUserData;

    //! Pointer to raw trajectory storage
    UINTEGER                       *ptrarRawPopulations;

#ifdef BOOST_NO_CXX11_HDR_ATOMIC
    //! Atomic flag to interrupt execution
    boost::atomic<bool>            bInterruptRequested;
#else
    //! Atomic flag to interrupt execution
    std::atomic<bool>              bInterruptRequested;
#endif

  //////////////////////////////
  // Methods
  protected:

    //! Triggered after attaching a simulation engine to this instance
    bool OnAttach()
    {
      setupTiming();
      return true;
    }

  //////////////////////////////
  // Methods
  public:

    /////////////////////////////////
    // Simulation engine interface

    // Attach to simulation engine
    bool attachPSSA(PSSA * pPSSA) { m_ptrPSSA = pPSSA; return this->OnAttach(); };

    // Detach from simulation engine
    PSSA * detachPSSA() { PSSA * pPSSA = m_ptrPSSA; m_ptrPSSA = NULL; resetOutput(); return pPSSA; };

    // Check if this instance is attache to a simulation engine
    bool isAttached() const { return (bool)(NULL != m_ptrPSSA); };

    // Check if this instance is valid
    bool isValid();

    // Process user settings
    bool processSettings();
#ifdef HAVE_LIBSBML
    /////////////////////////////////
    // SBML interface

    //! load a new model from SBML file
    bool readSBMLFile(const STRING & strFilePath)
    {
      // Try to load SBML file
      LIBSBML_CPP_NAMESPACE::SBMLDocument *
        ptrSBMLDocument = LIBSBML_CPP_NAMESPACE::readSBML(strFilePath.c_str());

      return parseSBMLDocument(ptrSBMLDocument);
    }

    /**
     * Convert an SBML model into its pSSAlib representation
     * 
     * @param ptrSBMLDocument a pointer to an @c SBMLDocument object
     * @return @true on success, @false otherwise
     */
    bool parseSBMLDocument(LIBSBML_CPP_NAMESPACE::
      SBMLDocument * ptrSBMLDocument);
#endif
    /////////////////////////////////
    // Getters & setters

    /**
     * Get number of spatial dimensions
     * 
     * @return Number of spatial dimensions defined or 0
     */
    BYTE getDimsCount() const
    {
      return m_uDims;
    };

    /**
     * Get length of all spatial dimensions in an array
     * 
     * @return Array of spatial dimension lengths
     */
    UINTEGER * getDims() const
    {
      return m_arunDims;
    };

    /**
     * Set the number of spatial dimensions and their lengths
     * 
     * @param dims Number of spatial dimensions
     * @param ... Length of each dimension (number of arguments must match @c dims)
     */
    void setDims(BYTE dims, ...)
    {
      va_list vl;
      if(NULL != m_arunDims)
        delete [] m_arunDims;
      m_uDims = dims;
      m_arunDims = new UINTEGER[m_uDims];
      va_start(vl, dims);
      for(BYTE di = 0; di < m_uDims; ++di)
        m_arunDims[di] = va_arg(vl, UINTEGER);
      va_end(vl);
    }

    /**
     * Set the number of spatial dimensions and their lengths
     * 
     * @param dims Number of spatial dimensions
     * @param arDims Array of spatial dimension lengths
     */
    void setDims(BYTE dims, UINTEGER *arDims)
    {
      if(NULL != m_arunDims)
        delete [] m_arunDims;
      m_uDims = dims;
      m_arunDims = new UINTEGER[m_uDims];
      for(BYTE di = 0; di < m_uDims; ++di)
        m_arunDims[di] = arDims[di];
    }

    //! Return number of species
    UINTEGER getNumSpecies()
    {
      if (NULL != getDataModel())
        return getDataModel()->getSpeciesCount();
      else 
      {
        PSSA_ERROR(this, << "DataModel not loaded"<< std::endl);
        return 0;
      }
    };

#ifdef HAVE_MPI
    //! Get the first index in the block
    INTEGER getBlockStart() const { return m_nBlockStart; }

    //! Get block size
    INTEGER getBlockSize() const { return m_nBlockSize; }
#endif

    //! Check whether the update mode is set externally
    bool getDelayedUpdate() const { return (bool)(unFlags & sfDelayedUpdate); };

    //! Get the underlying PSSA object
    PSSA* getPSSA() const { return m_ptrPSSA; };

    //! Get the data model associated with PSSA object
    DataModel* getDataModel() const { return m_ptrPSSA->ptrData; };

    /**
     * Get the @c datamodel::detail::Model object containing current model definition
     * 
     * @return A reference to the @c datamodel::detail::Model object containing current model definition
     */
    detail::Model & getModel()
    {
      return m_Model;
    }

    /////////////////////////////////
    // Functions used within the
    // simulation loop

    //! True if the simulation has not reached dTimeEnd
    bool isRunning() { return (dTimeSimulation < dTimeEnd); };

    // Output the results to file
    void doOutput(void);

    // FIXME : callback to update delayed reactions (consuming AND non-consuming(!) )
    bool UpdateCallback(void);

    /////////////////////////////////
    // Output

    //! Reset internal output streams
    void resetOutput();

    //! Converts output flags to stream index
    static USHORT outputFlagToStreamIndex(const OutputFlags of)
    {
      switch(of)
      {
      case ofLog:
        return 0;
        break;
      case ofStatus:
        return 1;
        break;
      case ofTrajectory:
        return 2;
        break;
      case ofFinalPops:
        return 3;
        break;
      case ofTimePoints:
        return 4;
        break;
      case ofTiming:
        return 5;
        break;
      case ofSpeciesIDs:
        return 6;
        break;
      case ofMaskFile:
        return 7;
        break;
      default:
        return std::numeric_limits<USHORT>::max();
      break;
      }
    }

    //! Get the output stream corresponding to a given type of output
    OSTREAM & getOutputStream(const OutputFlags of);

    //! Set an output stream corresponding for an output type
    bool setOutputStreamBuf(const OutputFlags of, STREAMBUFFER * sb)
    {
      if((NULL != sb)&&(of > ofMaskLog)&&(of < ofMaskFile))
      {
        SHORT idx = outputFlagToStreamIndex(of);

        // free the internal file stream buffer, if any
        resetOutputStream(of);

        // return a pointer to the previous external stream buffer
        m_arPtrExternalBuffers[idx] = sb;

        return true;
      }

      // fail
      return false;
    }

    //! Reset internal output stream for a given type of output
    void resetOutputStream(const OutputFlags of)
    {
      if((of > ofMaskLog)&&(of < ofMaskFile))
      {
        SHORT idx = outputFlagToStreamIndex(of);

        if(NULL != m_arPtrFileBuffers[idx])
        {
          m_arPtrFileBuffers[idx]->pubsync();
          delete m_arPtrFileBuffers[idx];
          m_arPtrFileBuffers[idx] = NULL;
        }
      }
    }

    //! Report a message via the output streams
    OSTREAM & report()
    {
      return getOutputStream(ofLog);
    }

    //! Check if logging is enabled for this message type
    bool isLoggingOn(const UINTEGER of) const
    {
      return ((unOutputFlags & of) == of);
    }

    /////////////////////////////////
    // Timing

    /**
     * Setup variables for timing routines.
     */
    void setupTiming();

    /**
     * Begin timing of the trial.
     * @param sample Current sample number
     * @param total Total number of samples in the series
     */
    bool beginTrial(UINTEGER sample);

    /**
     * Calculate time since last call to beginTrial.
     * @return time since last call to beginTrial in seconds
     */
    REAL_EXT endTrial();
  };

} } // close namespaces datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_SIMULATIONINFO_H_ */
