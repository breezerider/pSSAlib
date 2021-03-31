/**
 * @file analyzer.cpp
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
 * pSSAlib Command Line Interface: Analyzer component
 */

#include <iostream>
#include <boost/scoped_ptr.hpp>

#include <sbml/SBMLTypes.h>

#include "PSSA.h"
#include "util/FileSystem.h"
#include "util/Timing.h"
#include "util/IO.hpp"
#include "util/SimulationDataSource.hpp"
#include "../libpssa/dependencies/sha1/sha1.h"

#include "CmdLineOptions.hpp"
using namespace pssalib::program_options;

#include <boost/uuid/sha1.hpp>

// null stream
pssalib::io::null_streambuf< STRING::value_type > nullBuffer;
OSTREAM nullStream(&nullBuffer);

// Debug backtrace
#ifdef __linux__
#include <execinfo.h>
#include <signal.h>

void handler(int sig)
{
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, 2);
  exit(1);
}
#endif

//
class AnalyzerData
{
public:
  typedef std::pair< std::vector<UINTEGER>::const_iterator,
                     std::vector<UINTEGER>::const_iterator >
                     PAIR_RANGE_UINTEGER;
                     
  typedef std::pair< std::vector<STRING>::const_iterator,
                     std::vector<STRING>::const_iterator >
                     PAIR_RANGE_STRING;

  typedef ProgramOptionsAnalyzer::AnalyzerFormats AnalyzerFormats;

protected:
  std::vector<UINTEGER> m_arunSpecies,
                        m_arunSubvolumes;
  std::vector<STRING>   m_arstrSpecies;

  REAL                  m_dTimeInitial,
                        m_dTimeBegin,
                        m_dTimeEnd,
                        m_dTimeStep;

  UINTEGER              m_unSamples;

  bool                  m_bQuiet,
                        m_bVerbose;

  AnalyzerFormats       m_Format;

private:
  /**
  * Read species ids from a file and store them in a vector
  * 
  * @param strInputFile Path to input file
  * @param arSpecies Output vector to which species ids are appended
  */
  bool intializeSpeciesIds(const STRING & strInputFile,
                          std::vector<STRING> & arSpecies)
  {
    // input species names file
    FILESTREAMBUFFER fsbSpeciesNames;
    if(!fsbSpeciesNames.open(strInputFile.c_str(), std::ios::in))
    {
      PSSALIB_MPI_CERR_OR_NULL << "could not open file '"
        << strInputFile << "'.\n";
      return false;
    }

    // read species ids from file
    STRING strTemp;
    ISTREAM isData(&fsbSpeciesNames);
    while(isData >> strTemp)
    {
      if(!strTemp.empty())
        arSpecies.push_back(strTemp);
    }
    fsbSpeciesNames.close();

    return true;
  }

  void addSpecies(UINTEGER idx, STRING id)
  {
    if(m_arunSpecies.size() != 0)
    {
      std::vector<UINTEGER>::iterator it =
        std::lower_bound(m_arunSpecies.begin(), m_arunSpecies.end(), idx);

      if(idx == *it) return; // skip if already present

      std::iterator_traits<
        std::vector<UINTEGER>::iterator >::difference_type
        d = std::distance(m_arunSpecies.begin(), it);
      m_arunSpecies.insert(it, idx);
      m_arstrSpecies.insert(m_arstrSpecies.begin() + d, id);
    }
    else
    {
      m_arunSpecies.push_back(idx);
      m_arstrSpecies.push_back(id);
    }
  }

  void addSubvolume(UINTEGER idx)
  {
    m_arunSubvolumes.push_back(idx);
  }

public:

  bool initialize(ProgramOptionsAnalyzer & poAnalyzer, ProgramOptionsSimulator & poSimulator)
  {
    if(poAnalyzer.isVerboseSet()&&poAnalyzer.isQuietSet())
    {
      PSSALIB_MPI_CERR_OR_NULL << "Conflicting output definitions: both 'verbose' and "
      "'quiet' flags set, however, the later has priority over the former.\n";
    }
    
    // Output modifiers
    if(poAnalyzer.isVerboseSet())
    {
      m_bVerbose = true;
      m_bQuiet   = false;
    }
    if(poAnalyzer.isQuietSet())
    {
      m_bVerbose = false;
      m_bQuiet   = true;
    }
    
    // Format
    m_Format = poAnalyzer.getFormat();

    // check if number of trials does not exceed that in the simulation dataset
    m_unSamples = poSimulator.getNumSamples();
    if(poAnalyzer.isNumSamplesSet())
    {
      if(poAnalyzer.getNumSamples() > poSimulator.getNumSamples())
      {
        PSSALIB_MPI_CERR_OR_NULL
          << "Requested number of samples " << poAnalyzer.getNumSamples()
          << " exceeds available in the simulation dataset ("
          << poSimulator.getNumSamples()
          << "), terminating the analyzer.\n";
        return false;
      }
      m_unSamples = poAnalyzer.getNumSamples();
    }

    // check timing
    m_dTimeInitial = poSimulator.getTimeBegin();
    m_dTimeStep = poSimulator.getTimeStep();

//     m_dTimeBegin = poSimulator.getTimeBegin();
//     if(poAnalyzer.isTimeBeginSet())
//     {
      if((poAnalyzer.getTimeBegin() < poSimulator.getTimeBegin())||
        (poAnalyzer.getTimeBegin() > poSimulator.getTimeEnd()))
      {
        PSSALIB_MPI_CERR_OR_NULL << "Requested beginning time is not within simulation dataset range: "
                        << poAnalyzer.getTimeBegin() << " not in [" << poSimulator.getTimeBegin()
                        << "; " << poSimulator.getTimeEnd() << "], terminating the analyzer.\n";
        return false;
      }
      m_dTimeBegin = pssalib::timing::getAdjTimePointHi(poSimulator.getTimeBegin(),
                                                        poSimulator.getTimeEnd(),
                                                        poSimulator.getTimeStep(),
                                                        poAnalyzer.getTimeBegin());
//     }

    m_dTimeEnd = poSimulator.getTimeEnd();
    if(poAnalyzer.isTimeEndSet())
    {
      if((poAnalyzer.getTimeEnd() < poSimulator.getTimeBegin())||
        (poAnalyzer.getTimeEnd() > poSimulator.getTimeEnd()))
      {
        PSSALIB_MPI_CERR_OR_NULL << "Requested final time is not within simulation dataset range: "
                        << poAnalyzer.getTimeEnd() << " not in [" << poSimulator.getTimeBegin()
                        << "; " << poSimulator.getTimeEnd() << "], terminating the analyzer.\n";
        return false;
      }
      m_dTimeEnd = pssalib::timing::getAdjTimePointLo(poSimulator.getTimeBegin(),
                                                      poSimulator.getTimeEnd(),
                                                      poSimulator.getTimeStep(),
                                                      poAnalyzer.getTimeEnd());
    }

    // check species
    {
      STRING speciesNamesPath;
      pssalib::util::makeFilePath(poAnalyzer.getInputPath(),
                                  pssalib::datamodel::SimulationInfo::arFileNames[
                                  pssalib::datamodel::SimulationInfo::outputFlagToStreamIndex(
                                    pssalib::datamodel::SimulationInfo::ofSpeciesIDs)],
                                  speciesNamesPath);

      std::vector<STRING> arSpeciesInput;
      if(!intializeSpeciesIds(speciesNamesPath, arSpeciesInput))
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error : failed to load species names '" << speciesNamesPath << "'.\n";
        return false;
      }

      if(NULL != poAnalyzer.getSpecies())
      {
        for(std::vector<STRING>::const_iterator itA = poAnalyzer.getSpecies()->begin();
            itA != poAnalyzer.getSpecies()->end(); itA++)
        {
          std::vector<STRING>::const_iterator itF =
            std::find(arSpeciesInput.begin(), arSpeciesInput.end(), *itA);
          if(arSpeciesInput.end() != itF)
          {
            addSpecies(itF-arSpeciesInput.begin(), *itA);
          }
          else
          {
            PSSALIB_MPI_CERR_OR_NULL << "Species '" << (*itA) << "' is not "
                            "present in the original data set and will be "
                            "ignored by the analyser.\n";
          }
        }
      }
      else
      {
        for(std::vector<STRING>::const_iterator itS = arSpeciesInput.begin();
            itS != arSpeciesInput.end(); itS++)
          addSpecies(itS - arSpeciesInput.begin(), *itS);
      }
    }

    return true;
  }

  UINTEGER getNumSamples() const
  {
    return m_unSamples;
  }

  REAL getTimeInitial() const
  {
    return m_dTimeInitial;
  }

  REAL getTimeBegin() const
  {
    return m_dTimeBegin;
  }

  REAL getTimeEnd() const
  {
    return m_dTimeEnd;
  }

  REAL getTimeStep() const
  {
    return m_dTimeStep;
  }

  AnalyzerFormats getFormat() const
  {
    return m_Format;
  }

  bool isVerbose() const
  {
    return m_bVerbose;
  }

  bool isQuiet() const
  {
    return m_bQuiet;
  }

  const std::vector<UINTEGER> & getSpeciesIdx() const
  {
    return m_arunSpecies;
  }

  const std::vector<STRING> & getSpeciesIds() const
  {
    return m_arstrSpecies;
  }

  const std::vector<UINTEGER> & getSubvolumesIdx() const
  {
    return m_arunSubvolumes;
  }
};

//! GNU Plot output header for a single trajectory
const char *arGnuPlotHeaderTrajectories = \
  "set title '%s'\n" \
  "set ylabel '%s'\n" \
  "set xlabel '%s'\n" \
  "set grid\n";

/**
 * Generate trajectories form an input data set
 * 
 * @param strInputPath Path to the input data set
 * @param strOutputPath Output path
 * @param analyzerData An instance of @c AnalyzerData
 */
bool generateTrajectories(const STRING & strInputPath,
                          const STRING & strOutputPath,
                          const AnalyzerData & analyzerData)
{
  UINTEGER unTimePoints = pssalib::timing::getNumTimePoints(analyzerData.getTimeBegin(),
                                                            analyzerData.getTimeEnd(),
                                                            analyzerData.getTimeStep()),
           unTPSkip     = pssalib::timing::getNumTimePoints(analyzerData.getTimeInitial(),
                                                            analyzerData.getTimeBegin(),
                                                            analyzerData.getTimeStep());

  // Output ranges
  // time
  std::pair<UINTEGER, UINTEGER>     rangeTime(unTPSkip, unTPSkip + unTimePoints);

  // subvolumes and species
  std::pair< const UINTEGER *, const UINTEGER *>
    rangeSpecies(&(*(analyzerData.getSpeciesIdx().begin())),
                 &(*(analyzerData.getSpeciesIdx().end()))),
    rangeSubvolumes(&(*(analyzerData.getSubvolumesIdx().begin())),
                    &(*(analyzerData.getSubvolumesIdx().end())));

  // Output formatter
  OutputFormatter * fmt = NULL;
  bool bGnuPlotCommands = false;
  switch(analyzerData.getFormat())
  {
    case ProgramOptionsAnalyzer::afCSV:
    {
      // Header
      STRINGSTREAM ssTemp;
      ssTemp << "Time,";
      std::copy(analyzerData.getSpeciesIds().begin(),
                analyzerData.getSpeciesIds().end(),
                std::ostream_iterator<STRING>(ssTemp, ","));
      ssTemp << "\n";

      fmt = new CSVOutputFormatter(ssTemp.str(), analyzerData.getTimeStep(), analyzerData.getTimeBegin());
    }
    break;
    case ProgramOptionsAnalyzer::afGnuPlot:
    {
      // Header
      STRINGSTREAM ssTemp;
      ssTemp << "# Time ";
      std::copy(analyzerData.getSpeciesIds().begin(),
                analyzerData.getSpeciesIds().end(),
                std::ostream_iterator<STRING>(ssTemp, " "));
      ssTemp << "\n";

      fmt = new GnuplotOutputFormatter(ssTemp.str(), analyzerData.getTimeStep(), analyzerData.getTimeBegin());
      bGnuPlotCommands = analyzerData.isVerbose(); 
    }
    break;
    case ProgramOptionsAnalyzer::afVTK:
    {
      fmt = new VTKOutputFormatter(0, NULL, analyzerData.getSpeciesIds());
    }
    break;
    default:
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error: unknown formatter specification '" << (UINTEGER)analyzerData.getFormat() << "'\n";
      return false;
    }
  }

  // additional output
  if(bGnuPlotCommands)
  {
    PSSALIB_MPI_COUT_OR_NULL << "GnuPlot script to plot the results:\n"
                             << boost::format(arGnuPlotHeaderTrajectories)
                             % "Population Trajectories" % "Species Population" % "Time, s";
  }

  // Result
  SimulationDataSource sds;
  for(UINTEGER n = 0; n < analyzerData.getNumSamples(); ++n)
  {
    STRING strFileName((BOOSTFORMAT
      (pssalib::datamodel::SimulationInfo::arFileNames[
        pssalib::datamodel::SimulationInfo::outputFlagToStreamIndex(
          pssalib::datamodel::SimulationInfo::ofTrajectory)
      ]) % n).str()),
           strFilePath;

    pssalib::util::makeFilePath(strInputPath,
                                strFileName,
                                strFilePath);

    if(!sds.load(strFilePath, rangeTime, rangeSpecies, rangeSubvolumes))
    {
      PSSALIB_MPI_CERR_OR_NULL << "failed to process file '" << strFilePath  << "'.\n";
      return false;
    }

    strFileName = (boost::format("trajectory_%d") % n).str();

    pssalib::util::makeFilePath(strOutputPath,
                                strFileName,
                                strFilePath);

    if(!sds.store(strFilePath, *fmt))
    {
      PSSALIB_MPI_CERR_OR_NULL << "failed to store trajectory '" << strFilePath  << "'.\n";
      return false;
    }

    if(bGnuPlotCommands)
    {
      PSSALIB_MPI_COUT_OR_NULL << "\nplot ";

      for(std::size_t i = 0; ;)
      {
        PSSALIB_MPI_COUT_OR_NULL << "'" << strFilePath << "' using 1:" << i+2 <<
          " with line lt " << i << " lw 2 title '" << analyzerData.getSpeciesIds().at(i) << "'";
        ++i;
        if(analyzerData.getSpeciesIds().size() != i)
          PSSALIB_MPI_COUT_OR_NULL << ", \\\n";
        else
        {
          PSSALIB_MPI_COUT_OR_NULL << "\n";
          break;
        }
      }
    }
  }

  return true;
}

//! GNU Plot output header for an averaged trajectory
const char *arGnuPlotHeaderAverageTrajectories = \
  "set title '%s'\n" \
  "set ylabel '%s'\n" \
  "set xlabel '%s'\n" \
  "set grid\n" \
  "set style fill transparent solid 0.5\n";

/**
 * Generate trajectories form an input data set
 * 
 * @param strInputPath Path to the input data set
 * @param strOutputPath Output path
 * @param analyzerData An instance of @c AnalyzerData
 */
bool computeAvgTrajectories(const STRING & strInputPath,
                            const STRING & strOutputPath,
                            const AnalyzerData & analyzerData)
{
  UINTEGER unTimePoints = pssalib::timing::getNumTimePoints(analyzerData.getTimeBegin(),
                                                            analyzerData.getTimeEnd(),
                                                            analyzerData.getTimeStep()),
           unTPSkip     = pssalib::timing::getNumTimePoints(analyzerData.getTimeInitial(),
                                                            analyzerData.getTimeBegin(),
                                                            analyzerData.getTimeStep());
  std::pair<UINTEGER, UINTEGER>     rangeTime(unTPSkip, unTimePoints);

  // subvolumes and species
  std::pair< const UINTEGER *, const UINTEGER *>
    rangeSpecies(&(*(analyzerData.getSpeciesIdx().begin())),
                 &(*(analyzerData.getSpeciesIdx().end()))),
    rangeSubvolumes(&(*(analyzerData.getSubvolumesIdx().begin())),
                    &(*(analyzerData.getSubvolumesIdx().end())));

  // Output formatter
  OutputFormatter * fmt = NULL;
  bool bGnuPlotCommands = false;
  switch(analyzerData.getFormat())
  {
    case ProgramOptionsAnalyzer::afCSV:
    {
      // Header
      STRINGSTREAM ssTemp;
      ssTemp << "Time,";
      std::copy(analyzerData.getSpeciesIds().begin(),
                analyzerData.getSpeciesIds().end(),
                std::ostream_iterator<STRING>(ssTemp, "_mean,"));
      std::copy(analyzerData.getSpeciesIds().begin(),
                analyzerData.getSpeciesIds().end(),
                std::ostream_iterator<STRING>(ssTemp, "_stddev,"));
      ssTemp << "\n";
      
      fmt = new CSVOutputFormatter(ssTemp.str(), analyzerData.getTimeStep(), analyzerData.getTimeBegin());
    }
    break;
    case ProgramOptionsAnalyzer::afGnuPlot:
    {
      // Header
      STRINGSTREAM ssTemp;
      ssTemp << "# Time ";
      std::copy(analyzerData.getSpeciesIds().begin(),
                analyzerData.getSpeciesIds().end(),
                std::ostream_iterator<STRING>(ssTemp, "_mean "));
      std::copy(analyzerData.getSpeciesIds().begin(),
                analyzerData.getSpeciesIds().end(),
                std::ostream_iterator<STRING>(ssTemp, "_stddev "));
      ssTemp << "\n";

      fmt = new GnuplotOutputFormatter(ssTemp.str(), analyzerData.getTimeStep(), analyzerData.getTimeBegin());
      bGnuPlotCommands = analyzerData.isVerbose(); 
    }
    break;
    case ProgramOptionsAnalyzer::afVTK:
    default:
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error: unknown formatter specification '" << (UINTEGER)analyzerData.getFormat() << "'\n";
      return false;
    }
  }
 
  // Result
  SimulationDataSource sdsInput, sdsResult(unTimePoints, 2 * analyzerData.getSpeciesIdx().size(), analyzerData.getSubvolumesIdx().size());
  const REAL dInvSamples = REAL(analyzerData.getNumSamples());
  for(UINTEGER n = 0; n < analyzerData.getNumSamples(); ++n)
  {
    STRING strFileName((BOOSTFORMAT
      (pssalib::datamodel::SimulationInfo::arFileNames[
        pssalib::datamodel::SimulationInfo::outputFlagToStreamIndex(
          pssalib::datamodel::SimulationInfo::ofTrajectory)
      ]) % n).str()),
           strFilePath;

    pssalib::util::makeFilePath(strInputPath,
                                strFileName,
                                strFilePath);

    if(!sdsInput.load(strFilePath, rangeTime, rangeSpecies, rangeSubvolumes))
    {
      PSSALIB_MPI_CERR_OR_NULL << "failed to process file '" << strFilePath  << "'.\n";
      return false;
    }

    for(UINTEGER t = 0; t < sdsInput.getTimePoints(); ++t)
    {
      for(UINTEGER sp = 0, nsp = sdsInput.getSpecies(); sp < sdsInput.getSpecies(); ++sp)
      {
        for(UINTEGER sv = 0; sv < sdsInput.getSubvolumes(); ++sv)
        {
          REAL x = sdsInput.at(t, sp, sv),
               x_mean = sdsResult.at(t, sp, sv);
          REAL dlt = x - x_mean;
          x_mean += dlt / REAL(n + 1);
          REAL dlt2 = x - x_mean;

          sdsResult.at(t, sp, sv) = x_mean;
          sdsResult.at(t, sp + nsp, sv) += dlt * dlt2 * dInvSamples;
        }
      }
    }
  }

  STRING strFilePath;;
  pssalib::util::makeFilePath(strOutputPath,
                              STRING("average_trajectory"),
                              strFilePath);

  if(!sdsResult.store(strFilePath, *fmt))
  {
    PSSALIB_MPI_CERR_OR_NULL << "failed to store trajectory '" << strFilePath  << "'.\n";
    return false;
  }

  if(bGnuPlotCommands)
  {
    PSSALIB_MPI_COUT_OR_NULL << "\nplot ";

    for(std::size_t i = 0; ;)
    {
      PSSALIB_MPI_COUT_OR_NULL << "'" << strFilePath << "' using 1:" << i+2 <<
        " with line lt " << i << " lw 2 title '" << analyzerData.getSpeciesIds().at(i) << "'";
      ++i;
      if(analyzerData.getSpeciesIds().size() != i)
        PSSALIB_MPI_COUT_OR_NULL << ", \\\n";
      else
      {
        PSSALIB_MPI_COUT_OR_NULL << "\n";
        break;
      }
    }
  }

  return true;
}

typedef struct tagPDFInfo
{
  UINTEGER idx, cnt;

  tagPDFInfo()
    : idx(0)
    , cnt(0)
  {
    // Do nothing
  }

  tagPDFInfo(UINTEGER i, UINTEGER c)
    : idx(i)
    , cnt(c)
  {
    // Do nothing
  }
} PDFInfo;

/**
 * Generate trajectories form an input data set
 * 
 * @param strInputPath Path to the input data set
 * @param strOutputPath Output path
 * @param arSpecies Species ids to be included in the output
 * @param arColIds Species ids to be included in the output
 * @param unSamples Number of samples to use for averaging
 */
bool computePDF(const STRING & strInputPath,
                const STRING & strOutputPath,
                const AnalyzerData & analyzerData)
{
  // Input
  SimulationDataSource sdsInput;

  // Output formatter
  OutputFormatter * fmt = NULL;
  switch(analyzerData.getFormat())
  {
    case ProgramOptionsAnalyzer::afCSV:
    {
      // Header
      STRINGSTREAM ssTemp;
      ssTemp << "Frequency,";
      std::copy(analyzerData.getSpeciesIds().begin(),
                analyzerData.getSpeciesIds().end(),
                std::ostream_iterator<STRING>(ssTemp, ","));
      ssTemp << "\n";

      fmt = new CSVOutputFormatter(ssTemp.str(), analyzerData.getTimeStep(), analyzerData.getTimeBegin());
    }
    break;
    case ProgramOptionsAnalyzer::afGnuPlot:
    {
      // Header
      STRINGSTREAM ssTemp;
      ssTemp << "# Frequency ";
      std::copy(analyzerData.getSpeciesIds().begin(),
                analyzerData.getSpeciesIds().end(),
                std::ostream_iterator<STRING>(ssTemp, " "));
      ssTemp << "\n";

      fmt = new GnuplotOutputFormatter(ssTemp.str(), analyzerData.getTimeStep(), analyzerData.getTimeBegin());
    }
    break;
    case ProgramOptionsAnalyzer::afVTK:
    default:
    {
      PSSALIB_MPI_CERR_OR_NULL << "Error: unknown formatter specification '" << (UINTEGER)analyzerData.getFormat() << "'\n";
      return false;
    }
  }

  // subvolumes and species
  std::pair< const UINTEGER *, const UINTEGER *>
    rangeSpecies(&(*(analyzerData.getSpeciesIdx().begin())),
                 &(*(analyzerData.getSpeciesIdx().end()))),
    rangeSubvolumes(&(*(analyzerData.getSubvolumesIdx().begin())),
                    &(*(analyzerData.getSubvolumesIdx().end())));

  STRING strFilePath;
  pssalib::util::makeFilePath(strInputPath,
                              PSSALIB_FILENAME_FINAL_TIME_POINT_POPULATIONS,
                              strFilePath);

  if(!sdsInput.load(strFilePath, std::pair<UINTEGER, UINTEGER>(), rangeSpecies, rangeSubvolumes))
  {
    PSSALIB_MPI_CERR_OR_NULL << "failed to process file '" << strFilePath  << "'.\n";
    return false;
  }

//   SHA1 sha1Hasher;
  boost::uuids::detail::sha1 sha1Hasher;
  unsigned int sha1Hash[5] = {0};
  std::array<unsigned int, 5> sha1HashArray;
  std::map<std::array<unsigned int, 5>, PDFInfo> mapPDFInfo;
  boost::scoped_array<REAL> arInputRow(new REAL[sdsInput.getSpecies() * sdsInput.getSubvolumes()]);

  for(UINTEGER i = 0; i < analyzerData.getNumSamples(); ++i)
  {
    for(UINTEGER j = 0; j < sdsInput.getSpecies(); ++j)
      for(UINTEGER k = 0, kj = j * sdsInput.getSubvolumes(); k < sdsInput.getSubvolumes(); ++k)
        arInputRow[k + kj] = sdsInput.at(i, j, k);

    sha1Hasher.reset();
    sha1Hasher.process_bytes(arInputRow.get(), sizeof(REAL) * sdsInput.getSpecies() * sdsInput.getSubvolumes());
    sha1Hasher.get_digest(sha1Hash);
//     unsigned int hash[5] = {0};
//     std::array<unsigned int, 5> sha1HashArray;
    std::copy(std::begin(sha1Hash), std::end(sha1Hash), sha1HashArray.begin());
//     const char * sha1Hash = (const char *)sha1Hasher.getDigest();
//     STRING strTemp(sha1Hash);
//     delete [] sha1Hash;
//     // Back to string
//     char buf[41] = {0};
// 
//     for (int i = 0; i < 5; i++)
//     {
//         std::sprintf(buf + (i << 3), "%08x", hash[i]);
//     }

    std::map<std::array<unsigned int, 5>, PDFInfo>::iterator
      itMap = mapPDFInfo.find(sha1HashArray);

    if(mapPDFInfo.end() == itMap)
    {
      mapPDFInfo[sha1HashArray] = PDFInfo(i, 1);
    }
    else
    {
      itMap->second.cnt++;
    }
  }

  // Result
  SimulationDataSource sdsResult(mapPDFInfo.size(), 1 + sdsInput.getSpecies() * sdsInput.getSubvolumes());

  REAL dInvN = 1.0 / (REAL)analyzerData.getNumSamples(); UINTEGER n = 0;
  for(std::map<std::array<unsigned int, 5>, PDFInfo>::iterator itMap = mapPDFInfo.begin();
      mapPDFInfo.end() != itMap; itMap++, n++)
  {
    sdsResult.at(n, 0) = ((REAL)itMap->second.cnt) * dInvN;
    for(UINTEGER k = 0, i = itMap->second.idx; k < sdsInput.getSubvolumes(); k++)
      for(UINTEGER j = 0, jk = k * sdsInput.getSpecies(); j < sdsInput.getSpecies(); ++j)
      {
        sdsResult.at(n, j+jk+1) = sdsInput.at(i, j, k);
      }
  }

  pssalib::util::makeFilePath(strOutputPath,
                              STRING("pdf"),
                              strFilePath);

  if(!sdsResult.store(strFilePath, *fmt))
  {
    PSSALIB_MPI_CERR_OR_NULL << "Error: failed to write file '" << strFilePath  << "'.\n";
    return false;
  }

  return true;
}

/**
 * Main
 */
int main(int argc, char** argv)
{
#ifdef __linux__
  signal(SIGSEGV, handler);
#endif
  PSSALIB_MPI_IO_INIT;

  // Parse command line and configuration file, if specified
  prog_opt::variables_map vm;
  ProgramOptionsAnalyzer poAnalyzer;
  if(!poAnalyzer.processCmdLineArgs(argc, argv, vm)) {
    return -127;
  }

  /////////////////////

  // process input path
  if(!pssalib::util::checkPath(poAnalyzer.getInputPath()))
  {
    if(!poAnalyzer.isQuietSet())
      PSSALIB_MPI_COUT_OR_NULL << "Could not find the data directory '" << poAnalyzer.getInputPath() 
                     << "'.\n";
    return -126;
  }

  prog_opt::variables_map vm_sim;
  ProgramOptionsSimulator poSimulator;
  STRING cfgPath;
  pssalib::util::makeFilePath(poAnalyzer.getInputPath(), STRING("pssa.cfg"), cfgPath);
  if(!poSimulator.parseConfigFile(cfgPath, vm_sim)) {
    return -125;
  }

  // make output path
  std::vector<STRING> arPath;
  {
    STRING strPathCurr;
    arPath.push_back(poAnalyzer.getOutputPath());
    if(!pssalib::util::makeDir(arPath,strPathCurr))
    {
      if(!poAnalyzer.isQuietSet())
        PSSALIB_MPI_CERR_OR_NULL << "Error : failed to create output directory '" << strPathCurr << "'.\n";
      return -124;
    }
  }

  // check if there are overlap with dataset for requested methods
  std::vector<pssalib::PSSA::EMethod> arMethods;
  if(0 != poAnalyzer.getMethods().size())
  {
    for(std::vector<pssalib::PSSA::EMethod>::const_iterator itA = poAnalyzer.getMethods().begin();
        itA != poAnalyzer.getMethods().end(); itA++)
    {
      std::vector<pssalib::PSSA::EMethod>::const_iterator itF =
        std::find(poSimulator.getMethods().begin(), poSimulator.getMethods().end(), *itA);
      if(poSimulator.getMethods().end() != itF)
      {
        arMethods.push_back(*itA);
      }
      else
      {
        if(!poAnalyzer.isQuietSet())
          PSSALIB_MPI_CERR_OR_NULL << "Method '" << pssalib::PSSA::getMethodName(*itA) << "' is not "
                        "present in the original data set and will be ignored by the analyser.\n";
      }
    }
  }
  else
  {
    arMethods.resize(poSimulator.getMethods().size());
    std::copy(poSimulator.getMethods().begin(), poSimulator.getMethods().end(), arMethods.begin());
  }
  
  AnalyzerData analyzerData;
  if(!analyzerData.initialize(poAnalyzer, poSimulator))
  {
    return -123;
  }

  //////////////////
  // Analyzer output
  if(!poAnalyzer.isQuietSet())
  {
    PSSALIB_MPI_COUT_OR_NULL << "\npSSAlib command line interface : Analyzer\n\n\n";

    PSSALIB_MPI_COUT_OR_NULL << "Analyze simulator output in :\n\tinput path : " 
                     << poAnalyzer.getInputPath() << "\nand store results in :\n\toutput path : " 
                     << poAnalyzer.getOutputPath() << "\nfor each of these methods:\n";
    if(0 != arMethods.size())
    {
      std::transform(arMethods.begin(), arMethods.end(),
                     std::ostream_iterator<STRING>(PSSALIB_MPI_COUT_OR_NULL, "\t"), 
                     &pssalib::PSSA::getMethodName);
      PSSALIB_MPI_COUT_OR_NULL << std::endl;
    }
    else
    {
      PSSALIB_MPI_CERR_OR_NULL << "No matching methods found in the simulation dataset, terminating the analyzer.\n";
      return -120;
    }

    PSSALIB_MPI_COUT_OR_NULL << "Use " << analyzerData.getNumSamples() << " trials and include data points between '"
                     << analyzerData.getTimeBegin() << "' and '" << analyzerData.getTimeEnd() << "' seconds.\n";

    if(0 == poAnalyzer.getResultOutputFlags())
    {
      PSSALIB_MPI_CERR_OR_NULL << "No results specified, terminating the analyzer.\n";
      return -119;
    }

    PSSALIB_MPI_COUT_OR_NULL << "Produce following results:\n";
    poAnalyzer.outputResultFlags(PSSALIB_MPI_COUT_OR_NULL,"\t");
    PSSALIB_MPI_COUT_OR_NULL << std::endl;

    if(NULL != poAnalyzer.getSpecies())
    {
      PSSALIB_MPI_COUT_OR_NULL << "Output results for species with following ids :\n";
      std::copy(analyzerData.getSpeciesIds().begin(),
                analyzerData.getSpeciesIds().end(),
                std::ostream_iterator< STRING >
                (PSSALIB_MPI_COUT_OR_NULL, "\t"));
      PSSALIB_MPI_COUT_OR_NULL << std::endl;
    }
    else
    {
      PSSALIB_MPI_COUT_OR_NULL << "Output results for all species.\n";
    }
  }

  // Analyze simulation results
  for(std::vector<pssalib::PSSA::EMethod>::const_iterator
      mi = arMethods.begin(); mi != arMethods.end(); mi++)
  {
    // construct input path
    STRING inputPathCurr, outputPathCurr;

    if(!poAnalyzer.isQuietSet())
      PSSALIB_MPI_COUT_OR_NULL << "analysing " << pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)*mi)
                     << " data...\n";

    arPath.clear();
    arPath.push_back(poAnalyzer.getInputPath());
    arPath.push_back(pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)*mi));
    if(!pssalib::util::makeDir(arPath,inputPathCurr,true))
    {
      if(!poAnalyzer.isQuietSet())
        PSSALIB_MPI_CERR_OR_NULL << "Error : failed to establish input directory '" << inputPathCurr << "'.\n";
      return -118;
    }

    arPath.clear();
    arPath.push_back(poAnalyzer.getOutputPath());
    arPath.push_back(pssalib::PSSA::getMethodName((pssalib::PSSA::EMethod)*mi));
    if(!pssalib::util::makeDir(arPath,outputPathCurr))
    {
      if(!poAnalyzer.isQuietSet())
        PSSALIB_MPI_CERR_OR_NULL << "Error : failed to create output directory '" << outputPathCurr << "'.\n";
      return -117;
    }

    if(ProgramOptionsAnalyzer::arTrajectory & poAnalyzer.getResultOutputFlags())
    {
      if(ProgramOptionsSimulator::srTrajectory & poSimulator.getResultOutputFlags())
      {
        if(!generateTrajectories(inputPathCurr, outputPathCurr, analyzerData))
        {
          if(!poAnalyzer.isQuietSet())
            PSSALIB_MPI_CERR_OR_NULL << "Error : failed to produce trajectories using data in '" << inputPathCurr << "'.\n";
          return -116;
        }
      }
      else
      {
        if(!poAnalyzer.isQuietSet())
          PSSALIB_MPI_CERR_OR_NULL << "Error : no data for trajectories is available in '" << inputPathCurr << "'.\n";
        return -116;
      }
    }

    if(ProgramOptionsAnalyzer::arAverageTrajectory & poAnalyzer.getResultOutputFlags())
    {
      if(ProgramOptionsSimulator::srTrajectory & poSimulator.getResultOutputFlags())
      {
        if(!computeAvgTrajectories(inputPathCurr, outputPathCurr, analyzerData))
        {
          if(!poAnalyzer.isQuietSet())
            PSSALIB_MPI_CERR_OR_NULL << "Error : failed to produce average trajectories using data in '" << inputPathCurr << "'.\n";
          return -115;
        }
      }
      else
      {
        if(!poAnalyzer.isQuietSet())
          PSSALIB_MPI_CERR_OR_NULL << "Error : no data for average trajectories is available in '" << inputPathCurr << "'.\n";
        return -115;
      }
    }

    if(ProgramOptionsAnalyzer::arPDF & poAnalyzer.getResultOutputFlags())
    {
      if(ProgramOptionsSimulator::srFinalPopulations & poSimulator.getResultOutputFlags())
      {
        if(!computePDF(inputPathCurr, outputPathCurr, analyzerData))
        {
          if(!poAnalyzer.isQuietSet())
            PSSALIB_MPI_CERR_OR_NULL << "Error : failed to produce probability density function using data in '" << inputPathCurr << "'.\n";
          return -114;
        }
      }
      else
      {
        if(!poAnalyzer.isQuietSet())
          PSSALIB_MPI_CERR_OR_NULL << "Error : no data for average trajectories is available in '" << inputPathCurr << "'.\n";
        return -115;
      }
    }
  }

#if defined(_DEBUG) && defined(WIN32)
  std::cin.get();
#endif

  return 0;
}
