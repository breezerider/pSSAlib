/**
 * @file SimulationInfo.cpp
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
 * Implements the marshalling datatype
 */

#include "../../include/PSSA.h"

#include "../../include/datamodel/DataModel.h"
#include "../../include/datamodel/DataModel_DM.h"

#include "../../include/update/UpdateModule.h"

#include "../../include/datamodel/SimulationInfo.h"

#include "../../include/util/MPIWrapper.h"
#include "../../include/util/FileSystem.h"
#include "../../include/util/Timing.h"
#include "../../include/util/IO.hpp"

#ifdef HAVE_LIBSBML
#include "../../include/datamodel/detail/SBMLHelper.hpp"
#endif

#include <boost/array.hpp>
#include <boost/math/special_functions/next.hpp>

namespace pssalib
{
namespace datamodel
{
  typedef pssalib::io::null_streambuf< STRING::value_type > NULLSTREAMBUFFER;

  /////////////////////////////////
  // Static class members
  const STRING SimulationInfo::arFileNames[] = {
    STRING(PSSALIB_FILENAME_LOG),
    STRING(PSSALIB_FILENAME_STATUS),
    STRING(PSSALIB_FILENAME_TRAJECTORY),
    STRING(PSSALIB_FILENAME_FINAL_TIME_POINT_POPULATIONS),
    STRING(PSSALIB_FILENAME_TIME_POINTS),
    STRING(PSSALIB_FILENAME_TIMING),
    STRING(PSSALIB_FILENAME_SPECIES_IDS)
  };

  /////////////////////////////////
  // Constructors

  //! Constructor
  SimulationInfo::SimulationInfo()
    : m_ptrPSSA(NULL)
    , m_unSampleCurrent(0)
    , m_uDims(0)
    , m_arunDims(NULL)
#ifdef HAVE_MPI
    , m_nBlockSize(0)
    , m_nBlockStart(0)
#endif
    , m_ptrOutputLine(NULL)
    , m_ptrRawTrajectory(NULL)
    , m_ptrCurrPopulation(NULL)
    , m_unOutputIdx(0)
    , m_unOutputMax(0)
    , unFlags(0)
#ifdef DEBUG
    , unOutputFlags(ofError|ofWarning|ofInfo|ofSpeciesIDs|ofStatus|ofLog)
#else
    , unOutputFlags(ofError|ofWarning|ofInfo|ofSpeciesIDs|ofStatus)
#endif
    , pArSpeciesIds(NULL)
    , unSamplesTotal(0)
    , dTimeCheckpoint(0.0)
    , dTimeStart(0.0)
    , dTimeStep(0.0)
    , dTimeEnd(0.0)
    , dTimeSimulation(0.0)
    , eBoundaryConditions(detail::BC_Invalid)
    , eInitialPopulation(detail::IP_Invalid)
    , ptrPopulationInitializer(NULL)
    , ptrPopulationInitializerUserData(NULL)
    , ptrarRawTrajectory(NULL)
    , bInterruptRequested(false)
  {
    memset(m_arPtrFileBuffers, 0, 7*sizeof(FILESTREAMBUFFER *));
    memset(m_arPtrExternalBuffers, 0, 7*sizeof(STREAMBUFFER *));
  }

  //! Copy constructor
  SimulationInfo::SimulationInfo(SimulationInfo & right)
    : m_ptrPSSA(right.m_ptrPSSA)
    , m_unSampleCurrent(right.m_unSampleCurrent)
    , m_uDims(0)
    , m_arunDims(NULL)
#ifdef HAVE_MPI
    , m_nBlockSize(right.m_nBlockSize)
    , m_nBlockStart(right.m_nBlockStart)
#endif
    , m_ptrOutputLine(right.m_ptrOutputLine)
    , m_ptrRawTrajectory(right.m_ptrRawTrajectory)
    , m_ptrCurrPopulation(right.m_ptrCurrPopulation)
    , m_unOutputIdx(right.m_unOutputIdx)
    , m_unOutputMax(right.m_unOutputMax)
    , unFlags(right.unFlags)
    , unOutputFlags(right.unOutputFlags)
    , pArSpeciesIds(right.pArSpeciesIds)
    , unSamplesTotal(right.unSamplesTotal)
    , dTimeCheckpoint(right.dTimeCheckpoint)
    , dTimeStart(right.dTimeStart)
    , dTimeStep(right.dTimeStep)
    , dTimeEnd(right.dTimeEnd)
    , dTimeSimulation(right.dTimeSimulation)
    , eBoundaryConditions(right.eBoundaryConditions)
    , eInitialPopulation(right.eInitialPopulation)
    , ptrPopulationInitializer(right.ptrPopulationInitializer)
    , ptrPopulationInitializerUserData(right.ptrPopulationInitializerUserData)
    , ptrarRawTrajectory(right.ptrarRawTrajectory)
    , bInterruptRequested(right.bInterruptRequested.load())
  {
    setDims(right.m_uDims, right.m_arunDims);
    memset(m_arPtrFileBuffers, 0, SimulationInfo::outputFlagToStreamIndex(ofMaskFile)*sizeof(FILESTREAMBUFFER *));
    memset(m_arPtrExternalBuffers, 0, SimulationInfo::outputFlagToStreamIndex(ofMaskFile)*sizeof(STREAMBUFFER *));
  }

  //! Destructor
  SimulationInfo::~SimulationInfo()
  {
    if(NULL != m_arunDims)
    {
      delete [] m_arunDims;
      m_arunDims = NULL;
      m_uDims = 0;
    }
    if(NULL != m_ptrOutputLine)
    {
      delete [] m_ptrOutputLine;
      m_ptrOutputLine = NULL;
    }
    if(NULL != m_ptrCurrPopulation)
    {
      delete [] m_ptrCurrPopulation;
      m_ptrCurrPopulation = NULL;
    }
    if(NULL != pArSpeciesIds)
    {
      delete pArSpeciesIds;
      pArSpeciesIds = NULL;
    }
    // free file buffers
    resetOutput();
  }

  //////////////////////////////
  // Methods
#ifdef HAVE_LIBSBML
  /**
    * Convert an SBML model into its pSSAlib representation
    * 
    * @param ptrSBMLDocument a pointer to an @c SBMLDocument object
    * @return @true on success, @false otherwise
    */
  bool SimulationInfo::parseSBMLDocument(LIBSBML_CPP_NAMESPACE::
    SBMLDocument * ptrSBMLDocument)
  {
    pssalib::datamodel::detail::SBMLHelper helper;

    if ((NULL == ptrSBMLDocument)||(NULL == ptrSBMLDocument->getModel()))
    {
      PSSA_ERROR(this, << "Illegal arguments: empty SBML document!" << std::endl);
      return false;
    }
    else if (ptrSBMLDocument->getNumErrors() > 0)
    {
      ptrSBMLDocument->printErrors(getOutputStream(ofLog));
      return false;
    }

    // Parse the model
    bool bResult = m_Model.assign(ptrSBMLDocument->getModel(), helper);

    pssalib::datamodel::detail::SBMLHelper::MSG_RANGE
      range = helper.getMessagesRange();

    for(; range.first != range.second; range.first++)
    {
      switch((*range.first).getType())
      {
      case pssalib::datamodel::detail::SBMLParserMessage::prtError:
      {
        PSSA_ERROR(this, << " at line " << (*range.first).getSourceLineNumber() << " : " << (*range.first).getErrorMessage() << std::endl);
        bResult = false;
      }
      break;
      case pssalib::datamodel::detail::SBMLParserMessage::prtWarning:
      {
        PSSA_WARNING(this, << " at line " << (*range.first).getSourceLineNumber() << " : " << (*range.first).getErrorMessage() << std::endl);
      }
      break;
      case pssalib::datamodel::detail::SBMLParserMessage::prtInfo:
      {
        PSSA_INFO(this, << " at line " << (*range.first).getSourceLineNumber() << " : " << (*range.first).getErrorMessage() << std::endl);
      }
      break;
      case pssalib::datamodel::detail::SBMLParserMessage::prtTrace:
      {
        PSSA_TRACE(this, << " at line " << (*range.first).getSourceLineNumber() << " : " << (*range.first).getErrorMessage() << std::endl);
      }
      break;
      default:
        PSSA_ERROR(this, << "<unknown error type> at line " << (*range.first).getSourceLineNumber() << " : " << (*range.first).getErrorMessage() << std::endl);
      break;
      }
    }

    if(!bResult)
      PSSA_ERROR(this, << "processing SBML model '" << m_Model.getId() << "' failed!\n");

    return bResult;
  }
#endif
  // Check this instance
  bool SimulationInfo::isValid()
  {
    // Check if this instance is attached
    if(NULL == m_ptrPSSA)
    {
      PSSA_ERROR(this, << "method should not be called on a detached instance." << std::endl);
      return false;
    }

    // Check output path
    if(!pssalib::util::checkPath(strOutput))
    {
      PSSA_WARNING(this, << "invalid output path: '" << strOutput << "', no information "
        "is output unless you redirect the respective streams." << std::endl);
      strOutput.clear();
    }

    // Check if number of trial is valid
    if(0 == unSamplesTotal)
    {
      PSSA_ERROR(this, << "number of trials must be a positive integer( > 0)." << std::endl);
      return false;
    }

    // Check if output parameters are valid
    if(isLoggingOn(ofTrajectory))
    {
      if(dTimeEnd < dTimeStart)
      {
        PSSA_ERROR(this, << "initial time point must preceed the final one" << std::endl);
        return false;
      }
#ifndef PSSALIB_ENGINE_CHECK
      if(0 >= dTimeStep)
      {
        PSSA_ERROR(this, << "time increment must be positive" << std::endl);
        return false;
      }
      if((dTimeEnd - dTimeStart) <= dTimeStep)
      {
        PSSA_ERROR(this, << "time increment must be less than the output time span" << std::endl);
        return false;
      }
#endif
    }

    return true;
  }

  // Process user settings
  bool SimulationInfo::processSettings()
  {
    STRINGSTREAM ssTemp;
    boost::iostreams::filtering_ostream os;
    os.push(pssalib::io::prefix_filter("(INFO) : "));
    os.push(ssTemp);

    os << "initializing species ids: \n";
    m_arSpeciesIdx.clear();
    if(NULL != pArSpeciesIds)
    {
      if(0 == pArSpeciesIds->size())
      {
        os << "output no species ids.\n";
      }
      else
      {
        os << "output some species ids.\n";
//         grouping::GroupingModule::MAP_SP2IDX::iterator	it;
        std::vector<STRING>::iterator it = pArSpeciesIds->begin();
        for(; it != pArSpeciesIds->end(); ++it)
        {
          if(PSSALIB_MPI_IS_MASTER)
          {
            os << "Processing element #" << it - pArSpeciesIds->begin() << " : " << (*it) << "\t";
          }

          UINTEGER s;
          try
          {
            s = m_Model.getSpeciesIndexById((*it));
          }
          catch(std::runtime_error & e)
          {
            if(PSSALIB_MPI_IS_MASTER)
            {
              os << "species identifier '" << (*it)
                << "' not found in model; position "
                << (it - pArSpeciesIds->begin())
                << std::endl;
            }
            return false;
          }

          if(PSSALIB_MPI_IS_MASTER)
          {
            os << "found species with identifier '" 
              << (*it) << "' ==> " << s << std::endl;
          }
          m_arSpeciesIdx.push_back(s);
        }
      }
    }
    else
    {
      os << "output all species ids.\n";
      for(UINTEGER i = 0; i < m_Model.getSpeciesCount(); i++)
      {
        //PSSA_INFO(ptrSimInfo, << " " << i << std::endl);
        m_arSpeciesIdx.push_back(i);
      }
    }

    os.flush();
    PSSA_INFO(this, << ssTemp.rdbuf() << std::flush;);

    return true;
  }

  //! Reset internal output streams
  void SimulationInfo::resetOutput()
  {
    for(SHORT of = ofMaskLog + 1; (of & ofMaskFile) > 0; of <<= 1 )
      resetOutputStream((OutputFlags)of);

    memset(m_arPtrFileBuffers, 0, outputFlagToStreamIndex(ofMaskFile) * sizeof(FILESTREAMBUFFER *));
  }

  //! Get an output stream corresponding for an output type
  OSTREAM & SimulationInfo::getOutputStream(const OutputFlags of)
  {
    // Null stream buffer
    static NULLSTREAMBUFFER nullBuffer;
    // Output stream
    static OSTREAM os(&nullBuffer);

    STREAMBUFFER * buffer = &nullBuffer;

    if((of > ofMaskLog)&&(of < ofMaskFile))
    {
      USHORT idx = outputFlagToStreamIndex(of);

      // is an external stream buffer assigned for this output?
      if(NULL != m_arPtrExternalBuffers[idx])
        buffer = m_arPtrExternalBuffers[idx];
      else if(!strOutput.empty())
      {
        // file stream buffer not yet allocated
        if(NULL == m_arPtrFileBuffers[idx])
        {
          STRING strFileName(arFileNames[idx]), strFilePath;
          switch(of)
          {
          case ofTrajectory:
            strFileName = (BOOSTFORMAT(strFileName) %
              m_unSampleCurrent).str();
          break;
#ifdef HAVE_MPI
          case ofStatus:
          case ofLog:
            strFileName = (BOOSTFORMAT(strFileName) %
              getMPIWrapperInstance().getRank()).str();
          break;
#endif
          default: // avoid compiler warnings
          break;
          }

          util::makeFilePath(strOutput, strFileName, strFilePath);

          m_arPtrFileBuffers[idx] = new FILESTREAMBUFFER();
          if((NULL != m_arPtrFileBuffers[idx])&&
             (NULL != m_arPtrFileBuffers[idx]->open(strFilePath.c_str(),
                                                  std::ios_base::out | std::ios_base::trunc)))
          {
            PSSA_INFO(this, << "file '" << strFilePath << "' was successfully created." << std::endl);
            buffer = m_arPtrFileBuffers[idx];
          }
          else
          {
            PSSA_ERROR(this, << "could not open file '" << strFilePath << "' for writing." << std::endl);
            if(m_arPtrFileBuffers[idx])
              delete m_arPtrFileBuffers[idx];
            m_arPtrFileBuffers[idx] = NULL;
          }
        }
        else
          buffer = m_arPtrFileBuffers[idx];
      }

    }

    os.rdbuf(buffer);

    return os;
  }

  //! Initialize timing variables
  void SimulationInfo::setupTiming()
  {
#ifdef __MACH__
    mach_timebase_info_data_t    sTimebaseInfo;
    mach_timebase_info(&sTimebaseInfo);
    m_tickRatio = ((REAL_EXT)  sTimebaseInfo.numer) /
      ((REAL_EXT) sTimebaseInfo.denom) / ((REAL_EXT) 1e9);
#elif defined(_WIN32)
    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency);
    m_tickRatio = (REAL_EXT)frequency.QuadPart;
#endif
  }

  //! Initializes timing variables
  bool SimulationInfo::beginTrial(UINTEGER sample)
  {
    // store sampling info
    m_unSampleCurrent = sample;
#ifdef HAVE_MPI
    PSSA_INFO(this, << "Commencing block sample " << m_unSampleCurrent - m_nBlockStart + 1 << " of " 
      << m_nBlockSize << " (trial " << m_unSampleCurrent + 1 << " of " << unSamplesTotal << ")" << std::endl);
#else
    PSSA_INFO(this, << "Commencing trial " << m_unSampleCurrent + 1 << " of " << unSamplesTotal << std::endl);
#endif
    // Allocate memory for output
    if(isLoggingOn(ofTrajectory)||isLoggingOn(ofRawTrajectory))
    {
      if(NULL != m_ptrCurrPopulation) delete [] m_ptrCurrPopulation;

      std::size_t sz = m_arSpeciesIdx.size() * m_ptrPSSA->ptrData->getSubvolumesCount();
      PSSA_INFO(this, << "Allocating " << sz << " bytes for current population buffer.\n");
      try
      {
        m_ptrCurrPopulation = new UINTEGER[sz];
      }
      catch (std::bad_alloc & e)
      {
        PSSA_ERROR(this, << "Error: " << e.what() << std::endl);
        return false;
      }
    }
    else
    {
      if(NULL != m_ptrCurrPopulation) delete [] m_ptrCurrPopulation;
      m_ptrCurrPopulation = NULL;
    }
    if(isLoggingOn(ofTrajectory))
    {
      if(NULL != m_ptrOutputLine) delete [] m_ptrOutputLine;

      std::size_t sz = (STRING(PSSALIB_TEXTOUTPUT_SUBVOLUMES_DELIMITER).size() +
                        (STRING(PSSALIB_TEXTOUTPUT_SPECIES_DELIMITER).size() +
                        boost::lexical_cast<STRING>(std::numeric_limits<UINTEGER>::max()).size()) * 
                        m_arSpeciesIdx.size()) * m_ptrPSSA->ptrData->getSubvolumesCount();
      PSSA_INFO(this, << "Allocating " << sz << " bytes for output line buffer.\n");
      try
      {
        m_ptrOutputLine = new typename STRING::value_type[sz];
      }
      catch (std::bad_alloc & e)
      {
        PSSA_ERROR(this, << "Error: " << e.what() << std::endl);
        return false;
      }
    }
    else
    {
      if(NULL != m_ptrOutputLine) delete [] m_ptrOutputLine;
      m_ptrOutputLine = NULL;
    }
    if(isLoggingOn(ofRawTrajectory))
      m_ptrRawTrajectory = ptrarRawTrajectory;
    else
      m_ptrRawTrajectory = NULL;
    m_unOutputIdx = 0;
    m_unOutputMax = timing::getNumTimePoints(dTimeStart, dTimeEnd, dTimeStep);
    if(m_unOutputMax > 0) --m_unOutputMax;

    // reset the timing
    dTimeSimulation = 0.0;
    dTimeCheckpoint = dTimeStart;

//     // initialize the offset for output
//     if(dTimeStart >= dTimeStep)
//       dTimeCheckpoint = dTimeStart - dTimeStep;
//     else
//     {
//       // Output the initial population
//       UINTEGER unOutputFlagsTemp = unOutputFlags; // to prevent spurious status updates
//       unOutputFlags &= ~ofStatus;
//       dTimeSimulation = dTimeStep; // pseudo-time step to make the output work
//       doOutput();
//       unOutputFlags = unOutputFlagsTemp;
// 
//       // reset the timing again
//       dTimeSimulation = 0.0;
//       dTimeCheckpoint = 0.0;
//     }

#ifdef __linux__
    clock_gettime(CLOCK_MONOTONIC, &m_trialStart);
#elif defined(__MACH__)
    m_trialStart = mach_absolute_time();
#elif defined(_WIN32)
    QueryPerformanceCounter(&m_trialStart);
#endif

    return true;
  }

  //! Measures time spent on last trial (in seconds) and resets the timers.
  REAL_EXT SimulationInfo::endTrial()
  {
#ifdef __linux__
    // get hi res time
    clock_gettime(CLOCK_MONOTONIC, &m_trialEnd);
#elif defined(__MACH__)
    // Get absolute time
    m_trialEnd = mach_absolute_time();
#elif defined(_WIN32)
    QueryPerformanceCounter(&m_trialEnd);
#endif
    // Output any pending data
    doOutput();

#ifdef HAVE_MPI
    PSSA_INFO(this, << "Concluding block sample " << m_unSampleCurrent - m_nBlockStart + 1 << " of "
      << m_nBlockSize << " (trial " << m_unSampleCurrent + 1 << " of " << unSamplesTotal << ")" << std::endl);
#else
    PSSA_INFO(this, << "Concluding trial " << m_unSampleCurrent + 1 << " of " << unSamplesTotal << std::endl);
#endif
    resetOutputStream(ofTrajectory);

    if(NULL != m_ptrOutputLine)
    {
      delete [] m_ptrOutputLine;
      m_ptrOutputLine = NULL;
    }

#ifdef __linux__
    // Output in seconds
    return  ((REAL_EXT) ( m_trialEnd.tv_sec - m_trialStart.tv_sec )) + ((REAL_EXT) ( m_trialEnd.tv_nsec - m_trialStart.tv_nsec )) / ((REAL_EXT) 1e9);
#elif defined(__MACH__)
    // Compute time span
    m_trialEnd -= m_trialStart;

    // Output in seconds
    return ((REAL_EXT) m_trialEnd ) * m_tickRatio;
#elif defined(_WIN32)
    // Output in seconds
    return (m_trialEnd.QuadPart - m_trialStart.QuadPart) / m_tickRatio;
#endif
  }

  //! Provides output to file
  void SimulationInfo::doOutput()
  {
    if(isLoggingOn(ofStatus))
    {
      SHORT percent = (dTimeSimulation * 100.0) / dTimeEnd;
      if(NULL != m_ptrPSSA->ptrProgrCallback)
        (*(m_ptrPSSA->ptrProgrCallback))(m_unSampleCurrent, unSamplesTotal, percent, m_ptrPSSA->ptrProgrCallbackUserData);
      else
      {
        static SHORT prevPercent = std::numeric_limits<SHORT>::max();

        if(prevPercent != percent)
        {
          prevPercent = percent;

          OSTREAM & osStatus = getOutputStream(ofStatus);

          osStatus.seekp(0);
          osStatus << "Progress : "
#ifdef HAVE_MPI
            "block sample " << m_unSampleCurrent - m_nBlockStart + 1 << " of " 
            << m_nBlockSize << " ("
#endif
            "sample " << m_unSampleCurrent + 1 << " of " << unSamplesTotal << 
#ifdef HAVE_MPI
            ")"
#endif
            " is " << percent << '%' << " done...\n";
        }
      }
    }

    typedef boost::array<STRING::value_type, 32> ConversionBufferType;
    static ConversionBufferType cBuf;
    static size_t szSubVolDelim = STRING(PSSALIB_TEXTOUTPUT_SUBVOLUMES_DELIMITER).size(),
                  szSpeciesDelim = STRING(PSSALIB_TEXTOUTPUT_SPECIES_DELIMITER).size();

    if(isLoggingOn(ofTrajectory)||isLoggingOn(ofRawTrajectory))
    {
#ifndef PSSALIB_ENGINE_CHECK
      // is it too early to begin the output?
      if(((dTimeSimulation < dTimeCheckpoint))||
          (dTimeEnd <= dTimeCheckpoint))
        return;

      // Temporary variables
      bool bFirst = (dTimeStart == dTimeCheckpoint);
      REAL dTemp = REAL(m_unOutputIdx)  * dTimeStep + dTimeStart;
      UINTEGER unTemp = std::floor((std::min(dTimeSimulation, dTimeEnd) - dTemp)/dTimeStep);
      if(bFirst) ++unTemp;

      if(dTimeSimulation > dTimeEnd)
      {
PSSA_TRACE(this, << "Beyond simulation time span\n\n\n");
dTimeCheckpoint = dTimeEnd;
        dTemp = dTimeEnd - dTemp - REAL(unTemp) * dTimeStep;
        if(dTemp > std::numeric_limits<REAL>::epsilon())
        {
PSSA_TRACE(this, << "Additional output line for the final timepoint\n\n\n");
          unTemp++;
        }
PSSA_TRACE(this, << "Total output lines = " << m_unOutputIdx + unTemp << "\n\n\n");
      }
      else
      {
        dTimeCheckpoint = dTemp + dTimeStep * REAL(unTemp);
      }

//       if(dTimeSimulation < dTimeEnd)
//       {
//         unTemp = 1 + std::floor((dTimeSimulation - dTemp)/dTimeStep);
// //         unTemp1 = 1 + std::floor((dTimeSimulation - dTimeStart)/dTimeStep);
// //         unTemp2 = m_unOutputIdx;//std::floor((dTimeCheckpoint - dTimeStart)/dTimeStep);
//       }
//       else
//       {
//         dTemp = (dTimeEnd - dTimeStart) - dTemp;
// //         unTemp1 = std::floor(dTemp/dTimeStep);
//         unTemp2 = m_unOutputIdx;/*std::floor((dTimeCheckpoint - dTimeStart)/dTimeStep) +
//           (UINTEGER)((dTimeCheckpoint - dTimeEnd - dTimeStep) >= (dTimeEnd + dTimeStart) * std::numeric_limits<REAL>::epsilon());*/
// //         std::cout << "dTemp=(dTimeEnd - dTimeStart) = " << dTemp
// //                   << "; dTemp/dTimeStep = " << dTemp/dTimeStep
// //                   << "; EPS = " << (dTimeEnd + dTimeStart) * std::numeric_limits<REAL>::epsilon()
// //                   << "; fmod = " << std::fmod(dTemp, dTimeStep)
// //                   << "; dTemp - REAL(unTemp1)*dTimeStep = " << dTemp - REAL(unTemp1)*dTimeStep
// //                   << "\n";
//         if((dTimeEnd + dTimeStart) * std::numeric_limits<REAL>::epsilon() >= (dTemp - REAL(unTemp1)*dTimeStep))
//           unTemp1 += 1;
//         else
//           unTemp1 += 2;
//       }
// 
// //       dTemp = std::min(dTimeSimulation, dTimeEnd) - dTimeStart;
// //       unTemp1 = std::floor(dTemp/dTimeStep);
// //       if((dTimeSimulation > dTimeEnd)&&
// //          (dTimeEnd != (dTimeStart + REAL(unTemp1)*dTimeStep)))
// //         ++unTemp1;
// 
// //       dTemp = dTimeCheckpoint - dTimeStart;
// 
// 
//       if(unTemp1 > unTemp2)
// //       if(dTimeCheckpoint < (dTimeEnd + dTimeStep))
//         unTemp = unTemp1 - unTemp2;
//       else
//       {
// // std::cout << "invalid timing!\n";
// // std::cout << "dTimeSimulation = " << dTimeSimulation << "; dTimeEnd = " << dTimeEnd << "; dTimeCheckpoint = " << dTimeCheckpoint << "; dTemp = " << dTemp << "; unTemp1 = " << unTemp1 << "; unTemp2 = " << unTemp2 << "\n";
//         return;
//       }
//       else if(dTimeCheckpoint == dTimeStart)
//         unTemp = 1;
//       else if((dTimeSimulation >= dTimeEnd)&&
//         (dTimeCheckpoint < dTimeEnd))
//         unTemp = 1;
// std::cout << "dTimeSimulation = " << dTimeSimulation << "; dTimeEnd = " << dTimeEnd << "; dTimeCheckpoint = " << dTimeCheckpoint << "; dTemp = " << dTemp << "; unTemp1 = " << unTemp1 << "; unTemp2 = " << unTemp2 << "\n";
      // Calculate the time difference
//       if(dTimeCheckpoint >= dTimeEnd)
//       {
// //         dTemp = dTimeCheckpoint - boost::math::float_next(dTimeEnd);
// //         std::cout << "float_distance(dTemp, dTimeStep) = " << boost::math::float_distance(dTemp, dTimeStep) << "\n";
// // std::cout << "past end time dTimeCheckpoint = " << dTimeCheckpoint << "; dTimeEnd = " << dTimeEnd << "\n";
// //         if(dTemp < dTimeStep)
// //         {
// // std::cout << "dTemp < dTimeStep is true : " << dTemp << " < " << dTimeStep << "\n";
//           dTemp = dTimeStep;
// //         }
// //         else
// //           dTemp = 0.0;
//       }
//       else
//       {
//         dTemp = std::min(dTimeSimulation, dTimeEnd) - dTimeCheckpoint;
// std::cout << "dTimeSimulation = " << dTimeSimulation << "; dTimeEnd = " << dTimeEnd << "; dTimeCheckpoint = " << dTimeCheckpoint << "; dTemp = " << dTemp << "\n";
//         if((dTimeSimulation >= dTimeEnd)&&(fabs(dTemp) < dTimeStep))
//         {
// // std::cout << "fabs(dTemp)/dTimeStep : " << (fabs(dTemp)/dTimeStep) << "\n";
// // std::cout << "fabs(1 - dTimeCheckpoint/dTimeEnd) : " << fabs(1 - dTimeCheckpoint/dTimeEnd) << "\n";
// //           if((fabs(dTemp)/dTimeStep > std::numeric_limits<REAL>::epsilon())||
// //              (fabs(1 - dTimeCheckpoint/dTimeEnd) < std::numeric_limits<REAL>::epsilon()))
//             if(dTemp > 0.0)
//               unTemp = 2;
//             else
//               unTemp = 1;
//         }
//         else
//         {
// //           if(dTimeSimulation >= dTimeEnd)
// //           {
// //             if(dTemp < dTimeStep)
// //               unTemp = 1;
// //           }
// //           else
// //         dTemp = dTimeStep*(1 +
// //           std::floor((dTemp - dTimeCheckpoint)/dTimeStep));
// std::cout << "dTemp/dTimeStep = " << dTemp/dTimeStep << "\n";
//             unTemp = 1 + std::floor(dTemp/dTimeStep);
//             if((dTimeCheckpoint + dTimeStep)
//         }
//         if(dTimeCheckpoint == dTimeStart)
//           ++unTemp;
//         if((dTimeSimulation >= dTimeEnd)&&
//            (dTemp != std::fmod(dTemp, dTimeStep)))
//         {
// std::cout << "additional output for final time point: dTimeSimulation = "
//   << dTimeSimulation << "; dTemp = " << dTemp << "; fmod = "
//   << std::fmod(dTemp, dTimeStep) << "\n";
// //           ++unTemp;
//         }
//       }

      if(unTemp > 0)
      {
#endif
        STRING::pointer pCurrOL = m_ptrOutputLine;
        UINTEGER        *pCurrPop = m_ptrCurrPopulation;
        size_t          szPop = m_ptrPSSA->ptrData->getSubvolumesCount()*m_arSpeciesIdx.size();

        // Output the population to stream
        for(UINTEGER svi = 0; svi < m_ptrPSSA->ptrData->getSubvolumesCount(); ++svi)
        {
          if ((NULL != pCurrOL)&&(0 != svi))
          {
            strncpy(pCurrOL, PSSALIB_TEXTOUTPUT_SUBVOLUMES_DELIMITER, szSubVolDelim);
            pCurrOL += szSubVolDelim;
          }

          for(UINTEGER si = 0; si < m_arSpeciesIdx.size(); ++si)
          {
            if ((NULL != pCurrOL)&&(0 != si))
            {
              strncpy(pCurrOL, PSSALIB_TEXTOUTPUT_SPECIES_DELIMITER, szSpeciesDelim);
              pCurrOL += szSpeciesDelim;
            }
            *pCurrPop = m_ptrPSSA->ptrData->getSubvolume(svi).population(m_arSpeciesIdx[si]);
            if (NULL != pCurrOL)
            {
              cBuf = boost::lexical_cast<ConversionBufferType>(*(pCurrPop++));
              size_t szSpeciesNum = strlen(cBuf.c_array());
              strncpy(pCurrOL, cBuf.c_array(), szSpeciesNum);
              pCurrOL += szSpeciesNum;
            }
          }
        }
        if (NULL != pCurrOL)
        {
          *(pCurrOL++) = '\n';
          *(pCurrOL++) = '\0';
        }
#ifndef PSSALIB_ENGINE_CHECK
        // Output to file
        for(; (unTemp > 0) && (m_unOutputIdx < m_unOutputMax); --unTemp, ++m_unOutputIdx)
        {
// std::cout << "dTimeCheckpoint = " << dTimeCheckpoint << "; unTemp = " << unTemp << "\n";
          // Advance in time
//           dTimeCheckpoint += dTimeStep;
//           dTemp           -= dTimeStep;
#endif
          if(NULL != m_ptrOutputLine) getOutputStream(ofTrajectory) << m_ptrOutputLine;
          if(NULL != m_ptrRawTrajectory)
          {
            memcpy(m_ptrRawTrajectory, m_ptrCurrPopulation, szPop*sizeof(UINTEGER));
            m_ptrRawTrajectory += szPop;
          }
#ifndef PSSALIB_ENGINE_CHECK
        }
        if(bFirst) --m_unOutputIdx;
// std::cout << "m_unOutputIdx = " << m_unOutputIdx << "\n\n";
//         while(dTemp >= dTimeStep);
      }
#endif
    }
  }

  //! Called within the time-sampling routine to
  //! perform an update for a delayed reaction.
  //! FIXME : test non-consuming reactions
  bool SimulationInfo::UpdateCallback( void )
  {
    // Output the population
    doOutput();

    unFlags |= sfDelayedUpdate;
    bool bResult = m_ptrPSSA->ptrUpdate->doUpdate(this);
    unFlags &= ~sfDelayedUpdate;

    return bResult;
  }

} } // close namespaces datamodel and pssalib
