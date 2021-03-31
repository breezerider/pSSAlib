/**
 * @file ProgramOptionsBase.hpp
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
 * Auxiliary classes for parsing command line arguments
 */

#ifndef PSSALIB_UTIL_PROGRAM_OPTIONS_BASE_HPP_
#define PSSALIB_UTIL_PROGRAM_OPTIONS_BASE_HPP_

#include <boost/tokenizer.hpp>

#include <boost/program_options.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

namespace prog_opt = boost::program_options;
namespace prop_tree = boost::property_tree;


  //////////////////////////////
  // Type for converting between literals and numeric constants
  template < const STRING::value_type delimiter >
  struct tagCLIOptionDelimitedList
  {
    STRING value;

    tagCLIOptionDelimitedList()
    {
      // Do nothing
    };

    tagCLIOptionDelimitedList(STRING::const_pointer s)
      : value(s)
    {
      // Do nothing
    };

    operator STRING() const
    {
      return value;
    }
    
//     operator std::string() const
//     {
//       return value;
//     }

    /**
     * Parse the value as a comma separated list and convert each entry
     * from literal to result type using a user-defined mapping.
     * @param mapping An @std::map object containing the relations
     *                between valid literal values and result value.
     * @param out Output vector object
     * @param allowMulti if @true, indicates that every non-empty list
     *                   element should added to the output vector. If
     *                   @false, parsing stops adter first match is found.
     * @param allowRepeats @true allows repeating instances of the same result
     *                     to be added to the output vector. If @false, repeats
     *                     are ignored.
     * @param reverseOrder @true will cause the input to be parsed in a
     *                     reverse order.
     */
    template<typename T>
    void parse(std::vector<T> & out,
              bool allowMulti = false, bool allowRepeats = true, bool reverseOrder = false)
    {
      static const STRING::value_type arSep[2] = {',', '\0'};
      boost::char_separator<STRING::value_type> sep(arSep);

      //std::cout << "Parsing string '" << value << "'" << std::endl;

      try
      {
        if(reverseOrder)
        {
          boost::tokenizer<
            boost::char_separator<STRING::value_type>,
            STRING::const_reverse_iterator,
            STRING > tok(value.rbegin(), value.rend(), sep);

          boost::tokenizer<
            boost::char_separator<STRING::value_type>,
            STRING::const_reverse_iterator,
            STRING >::iterator itTok = tok.begin();

          STRING temp;
          for(; (itTok != tok.end()); itTok++)
          {
            if((*itTok).empty())
              continue;

            temp.resize((*itTok).length());
            std::reverse_copy((*itTok).begin(), (*itTok).end(), temp.begin());

            T tempVal = boost::lexical_cast<T>(*itTok);
            if(allowRepeats||(out.end() == std::find(out.begin(), out.end(), tempVal)))
              out.push_back(tempVal);

            if(!allowMulti) break;

            temp.clear();
          }
        }
        else
        {
          boost::tokenizer<
            boost::char_separator<STRING::value_type>,
            STRING::const_iterator,
            STRING > tok(value.begin(), value.end(), sep);

          boost::tokenizer<
            boost::char_separator<STRING::value_type>,
            STRING::const_iterator,
            STRING >::iterator itTok = tok.begin();

          for(; (itTok != tok.end()); itTok++)
          {
            if((*itTok).empty())
              continue;

            T tempVal = boost::lexical_cast<T>(*itTok);
            if(allowRepeats||(out.end() == std::find(out.begin(), out.end(), tempVal)))
              out.push_back(tempVal);

            if(!allowMulti) break;
          }
        }
      }
      catch (boost::bad_lexical_cast &e)
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error parsing comma separated list: " << e.what() << "\n";
      }
    }

    /**
     * Parse the value as a comma separated list and convert each entry
     * from literal to result type using a user-defined mapping.
     * @param mapping An @std::map object containing the relations
     *                between valid literal values and result value.
     * @param out Output vector object
     * @param allowMulti if @true, indicates that every non-empty list
     *                   element should added to the output vector. If
     *                   @false, parsing stops when first match is found.
     * @param allowRepeats @true allows repeating instances of the same result
     *                     to be added to the output vector. If @false, repeats
     *                     are ignored.
     * @param reverseOrder @true will cause the input to be parsed in a
     *                     reverse order.
     */
    template<typename T>
    void parse(std::map<STRING, T> & mapping, std::vector<T> & out,
              bool allowMulti = false, bool allowRepeats = true, bool reverseOrder = false)
    {
      static const STRING::value_type arSep[2] = {',', '\0'};
      boost::char_separator<STRING::value_type> sep(arSep);

      //std::cout << "Parsing string '" << value << "'" << std::endl;

      if(reverseOrder)
      {
        boost::tokenizer<
          boost::char_separator<STRING::value_type>,
          STRING::const_reverse_iterator,
          STRING > tok(value.rbegin(), value.rend(), sep);

        boost::tokenizer<
          boost::char_separator<STRING::value_type>,
          STRING::const_reverse_iterator,
          STRING >::iterator itTok = tok.begin();

        STRING temp;
        for(; (itTok != tok.end()); itTok++)
        {
          if((*itTok).empty())
            continue;

          temp.resize((*itTok).length());
          std::reverse_copy((*itTok).begin(), (*itTok).end(), temp.begin());

          typename std::map<STRING, T>::iterator
            itMap = mapping.find(temp);

          if(mapping.end() == itMap)
          {
            PSSALIB_MPI_CERR_OR_NULL << "Error parsing comma separated list: invalid value '"
              << temp << "'.\n";
            return;
          }

          if(allowRepeats||(out.end() == std::find(out.begin(), out.end(), itMap->second)))
            out.push_back(itMap->second);

          if(!allowMulti) break;

          temp.clear();
        }
      }
      else
      {
        boost::tokenizer<
          boost::char_separator<STRING::value_type>,
          STRING::const_iterator,
          STRING > tok(value.begin(), value.end(), sep);

        boost::tokenizer<
          boost::char_separator<STRING::value_type>,
          STRING::const_iterator,
          STRING >::iterator itTok = tok.begin();

        for(; (itTok != tok.end()); itTok++)
        {
          if((*itTok).empty())
            continue;

          typename std::map<STRING, T>::iterator
            itMap = mapping.find((*itTok));

          if(mapping.end() == itMap)
          {
            PSSALIB_MPI_CERR_OR_NULL << "Error parsing comma separated list: invalid value '"
              << (*itTok) << "'.\n";
            return;
          }

          if(allowRepeats||(out.end() == std::find(out.begin(), out.end(), itMap->second)))
            out.push_back(itMap->second);

          if(!allowMulti) break;
        }
      }
    }
  };

  template < const STRING::value_type delimiter >
  std::istream& operator>>(std::istream& stream, const tagCLIOptionDelimitedList<delimiter> &o)
  {
      stream >> o.value;
      return stream;
  }

  template < const STRING::value_type delimiter >
  std::ostream& operator<<(std::ostream& stream, const tagCLIOptionDelimitedList<delimiter> &o)
  {
      stream << o.value;
      return stream;
  }
  
  template < const STRING::value_type delimiter >
  void validate(boost::any& v, std::vector<std::string> const& xs, tagCLIOptionDelimitedList<delimiter> *, long)
  {
    static std::stringstream ssTemp;
    static const STRING::value_type arSep[2] = {delimiter, '\0'};
    if (v.empty()) v = tagCLIOptionDelimitedList<delimiter>();

    try
    {
      tagCLIOptionDelimitedList<delimiter> & o = boost::any_cast<tagCLIOptionDelimitedList<delimiter> &>(v);
      ssTemp.str(std::string());
      std::copy(xs.begin(), xs.end(), std::ostream_iterator<std::string>(ssTemp, arSep));
      o.value += ssTemp.str();
    }
    catch(const boost::bad_any_cast & e)
    {
      PSSALIB_MPI_CERR_OR_NULL << "tagCLIOptionDelimitedList validator:" << e.what() << "\n";
    }
  }

namespace pssalib
{
namespace program_options
{
  // Helper type: print the first element of an std::pair
  template <class T>
  struct printPairFirst : public std::unary_function<T, void>
  {
    //! output stream
    std::ostream& os;

    //! delimiter
    const char * d;

    printPairFirst(std::ostream& _os, const char *_d) : os(_os), d(_d) {}

    void operator()(const T& e) const
    {
        os << e.first << d;//<< ", " << e.second << d;
    }
  };

  // Helper type: convert between compatible types
  template <typename T, typename K>
  struct converTypes : public std::unary_function<T, K>
  {
    K & operator()(T & t)
    {
        return (K &)t;
    }

    const K & operator()(const T & t) const
    {
        return (K &)t;
    }
  };
  
//   const STRING::value_type cliOptComma      = ',';
//   const STRING::value_type cliOptDimensions = 'x';
  typedef tagCLIOptionDelimitedList<','> CLIOptionCommaSeparatedList;
  typedef tagCLIOptionDelimitedList<'x'> CLIOptionDimensionList;


  //////////////////////////////

  /**
   * @class ProgramOptionsBase
   * @brief Base class handling command line arguments
   */
  class ProgramOptionsBase
  {
  //////////////////////////////
  // Attributes
  protected:
    //! Program options description
    prog_opt::options_description m_poDesc;

    //! Initialization flag
    bool m_bInitDone;

  //////////////////////////////
  // Constructors
  public:

    //! Constructor
    ProgramOptionsBase(const char * name)
      : m_poDesc(name)
      , m_bInitDone(false)
    {
      // Do nothing
    }

    //! Destructor
    virtual ~ProgramOptionsBase()
    {
      // Do nothing
    }

  //////////////////////////////
  // Methods
  protected:

    /**
     * Initiliaze the description object
     * @return @true if initialization is successfull, @false otherwise.
     */
    virtual bool initialize()
    {
      try
      {
        m_poDesc.add_options()
          ("help,h",                                                                  "Produce help message")
          ;

        return true;
      }
      catch (prog_opt::error & e)
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: " << e.what() << std::endl;
        return false;
      }
    }

    /**
     * Store data structures into the stream
     * @return @true if storing is successfull, @false otherwise.
     */
  //   virtual bool serialize(std::ostream &) const = 0;

  //////////////////////////////
  // Methods
  public:

    /**
     * Reset the object with pre-set values.
     */
    virtual void resetToDefault() = 0;

    /**
     * Initialize the object from a variable map built using the description object.
     * @param vm Variables map to store the parser results.
     * @return @true if parser succeeds, @false otherwise.
     */
    virtual bool parseVariableMap(prog_opt::variables_map & vm) = 0;

    /**
     * Get the underlying description object
     * @return program_option::description object initialized with pSSA CLI Simulator options
     */
    prog_opt::options_description & getDescription()
    {
      if(!m_bInitDone)
        m_bInitDone = initialize();
      return m_poDesc;
    }

    /**
     * Process a configuration file using the description object.
     * @param cfgPath Path to configuration file.
     * @param vm Variables map to store the parser results.
     * @return @true if parser succeeds, @false otherwise.
     */
    bool parseConfigFile(STRING &cfgPath, prog_opt::variables_map & vm)
    {
      if(!m_bInitDone)
        m_bInitDone = initialize();
      if(!m_bInitDone)
        return false;

      try
      {
        prog_opt::options_description cfg("Configuration");
        cfg.add(m_poDesc);

        prog_opt::store(prog_opt::parse_config_file<STRING::value_type>(cfgPath.c_str(), cfg), vm);

        // notify end of processing stage
        prog_opt::notify(vm);
      }
      catch (prog_opt::error &e)
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: " << e.what() << std::endl;
        return false;
      }

      return parseVariableMap(vm);
    }

    /**
     * Process command line arguments using the description object. Also, checks
     * if a configuration file was supplied and tries to parse it along with other
     * command line arguments with the latter having a precedence over the former.
     * @param argc number of arguments in the @argv array.
     * @param argv array of individual command line arguments (as provieded by the @main).
     * @param vm Variables map to store the parser results.
     * @return @true if parser succeeds, @false otherwise.
     */
    virtual bool processCmdLineArgs(int argc, char** argv, prog_opt::variables_map & vm)
    {
      if(!m_bInitDone)
        m_bInitDone = initialize();
      if(!m_bInitDone)
        return false;

      try
      {
        prog_opt::options_description cmdline("Command Line Options");
        cmdline.add_options()
          ("config-file,c", prog_opt::value<STRING>(), "Configuration file")
          ;
        cmdline.add(m_poDesc);

        prog_opt::options_description cfgfile("Configuration");
        cfgfile.add(m_poDesc);

        // First parse the command line, to retrieve the config_file parameter.
        prog_opt::store(prog_opt::parse_command_line(argc, argv, cmdline), vm);

        // check for help options
        if(/*(1 == argc)||*/(vm.count("help")))
        {
          PSSALIB_MPI_COUT_OR_NULL << cmdline << std::endl;
          return false;
        }

        // try to parse the config file.
        if(vm.count("config-file") > 0)
        {
          STRING cfgFileName = vm["config-file"].as<STRING>();

          prog_opt::store(prog_opt::parse_config_file<STRING::value_type>(cfgFileName.c_str(), cfgfile), vm);
          // Now read the command line again, to override the config file parameters where applicable.
          prog_opt::store(prog_opt::parse_command_line(argc, argv, cmdline), vm);
        }

        // notify end of processing stage
        prog_opt::notify(vm);
      }
      catch (prog_opt::error &e)
      {
        PSSALIB_MPI_CERR_OR_NULL << "Error processing program options: " << e.what() << std::endl;
        return false;
      }

      return parseVariableMap(vm);
    }
  };

} } // close program_options & pssalib namespaces

#endif /* PSSALIB_UTIL_PROGRAM_OPTIONS_BASE_HPP_ */
