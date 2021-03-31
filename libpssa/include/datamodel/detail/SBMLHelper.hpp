/**
 * @file SBMLHelper.hpp
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
 * Declares a helper class for storing SBML parser messages
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_SBML_HELPER_HPP_
#define PSSALIB_DATAMODEL_DETAIL_SBML_HELPER_HPP_

#include "../../typedefs.h"

namespace pssalib
{
namespace datamodel
{
namespace detail
{

  /**
   * @class SBMLParametersFilter
   * @brief Implements an SBML element filter for finding 
   * SBML LocalParameter or Parameter objects with a given id.
   */
  class SBMLParametersFilter : public LIBSBML_CPP_NAMESPACE::ElementFilter
  {
  private:
    //! Sought-for element id
    STRING strId;
  public:
      //! Constructor
    explicit SBMLParametersFilter(const STRING & id) 
        : ElementFilter()
        , strId(id)
      {
        // Do nothing
      }

      /**
      * Filters out SBML LocalParameter or Parameter
      * objects with a given id.
      * 
      * @copydoc ElementFilter::filter
      */
    virtual bool filter(const LIBSBML_CPP_NAMESPACE::SBase* element)
      {
        // check if it's a valid LocalParameter or Parameter object
        if((NULL != element) &&
           (0 == element->getPackageName().compare("core")) &&
           ((LIBSBML_CPP_NAMESPACE::SBML_LOCAL_PARAMETER==
           element->getTypeCode()) || (LIBSBML_CPP_NAMESPACE::
           SBML_PARAMETER==element->getTypeCode())))
        {
          // check if ids match
          return (0 == strId.compare(STRING(element->getId())));
        }

        return false;
      }
  };

  /**
   * @brief Parser log message.
   */
  class SBMLParserMessage
  {
  /////////////////////////////////////
  // Friend declarations
  public:
    friend class SBMLHelper;

  ////////////////////////////////
  // Data structures
  public:
    //! Flags
    typedef enum tagSBMLParserErrorType
    {
      prtTrace,
      prtInfo,
      prtWarning,
      prtError
    } SBMLParserErrorType;

  ////////////////////////////////
  // Attributes
  protected:
    SBMLParserErrorType errorType;
    UINTEGER unLineNumber;
    STRING strErrorMessage;

  public:
    SBMLParserMessage(SBMLParserErrorType t, UINTEGER n)
      : errorType(t)
      , unLineNumber(n)
    {
      // Do nothing
    };

    SBMLParserErrorType getType() const
    {
      return errorType;
    };

    UINTEGER getSourceLineNumber() const
    {
      return unLineNumber;
    };

    const STRING & getErrorMessage() const
    {
      return strErrorMessage;
    };
  };

  typedef struct tagUnitSpecification
  {
    REAL dExponent;
    REAL dMultiplier;
  } UnitSpecification;

  /**
   * \class AnnotationDefinition
   * \brief Declares a container for defining properties of annotation children of an SBML node.
   *
   */
  class SBMLHelper
  {
  /////////////////////////////////////
  // Data structures
  public:
    typedef std::pair< std::vector<SBMLParserMessage>::const_iterator, std::vector<SBMLParserMessage>::const_iterator > MSG_RANGE;

    typedef std::pair< LIBSBML_CPP_NAMESPACE::UnitKind_t,
    UnitSpecification > UNIT_SPEC;

    typedef std::map<LIBSBML_CPP_NAMESPACE::UnitKind_t,
    UnitSpecification > MAP_UNIT_SPEC;

  ////////////////////////////////
  // Attributes
  protected:

    STRINGSTREAM ssMessage;

    std::vector<SBMLParserMessage> arMessages;

  /////////////////////////////////////
  // Constructors
  public:
    //! Constructor
    SBMLHelper()
    {
      // Do nothing
    };

    //! Destructor
    ~SBMLHelper()
    {
      // Do nothing
    };

  /////////////////////////////////////
  // Methods
  public:

    /**
     * Wrapper over libSBML function
     * 
     * @param us A @UNIT_SPEC object
     * @return String representation of this unit
     */
  static const char * unitToString(const UNIT_SPEC & us)
    {
      return UnitKind_toString(us.first);
    }

    /**
     * Get the whole message sequence associated with this object
     * 
     * @return An @std::pair obeject containg first and last iterators in the sequence
     */
    MSG_RANGE getMessagesRange()
    {
      STRING str = ssMessage.str();
      if((str.length())&&(arMessages.size() > 0))
      {
        arMessages.back().strErrorMessage = str;
        ssMessage.str(STRING());
      }
      return MSG_RANGE(arMessages.begin(), arMessages.end());
    }

    /**
     * Append message to stream buffer
     * 
     * @param pet Message type
     * @param line Associated line number in SBML file
     * @return An @std::stringstream for collecting the error message
     */
    STRINGSTREAM & report(SBMLParserMessage::SBMLParserErrorType pet, UINTEGER line)
    {
      if(arMessages.size() > 0)
      {
        arMessages.back().strErrorMessage = ssMessage.str();
        ssMessage.str(STRING());
      }
      arMessages.push_back(SBMLParserMessage(pet, line));
      return ssMessage;
    };

    /**
     * Process a value in annotation annotation definition
     * 
     * @param str String to process
     * @param nodeBase Pointer to an @SBase object containing this annotation
     * @param result Placeholder for the evaluted value
     * @param units Units specifications that have to be matched
     * @return @true on success, @false otherwise
     */
    bool processValue(const STRING & str, const LIBSBML_CPP_NAMESPACE::SBase * nodeBase, REAL & result, const MAP_UNIT_SPEC & units)
    {
      result = 0.0;
      // first try to parse it as a number
      try
      {
        result = boost::lexical_cast<REAL>(str);
        report(SBMLParserMessage::prtTrace, nodeBase->getLine()) << "'" << str << "' = " << result;
        return true;
      }
      catch (boost::bad_lexical_cast &)
      {
        // not a plain number
      }

      // then assume it is a param id and try to look it up
      if(NULL == nodeBase) // we need a base node for the look-up
        return false;

      // set up the filter
      SBMLParametersFilter filter(str);

      // perform document search
      LIBSBML_CPP_NAMESPACE::List * lofElements = NULL;
      for(BYTE b = 0; b < 2; ++b)
      {
        if(0 == b)
          // Local search
          lofElements = const_cast<LIBSBML_CPP_NAMESPACE::SBase*>(nodeBase)->getAllElements(&filter);
        else
          // Global search
          lofElements = const_cast<LIBSBML_CPP_NAMESPACE::ListOfParameters*>
            (nodeBase->getModel()->getListOfParameters())->getAllElements(&filter);

        for(UINTEGER li = 0; li < lofElements->getSize(); ++li)
        {
          LIBSBML_CPP_NAMESPACE::Parameter * param = reinterpret_cast<LIBSBML_CPP_NAMESPACE::Parameter*>(lofElements->get(li));

          if(NULL != param)
          {
            result = (REAL)param->getValue();

            LIBSBML_CPP_NAMESPACE::UnitDefinition * unitDefSI = NULL;
            if(!param->isSetUnits())
              report(SBMLParserMessage::prtWarning, nodeBase->getLine()) << "'" << str << "' has no units associated with it!";
            else
            {
              unitDefSI = param->getDerivedUnitDefinition();
              if(NULL != unitDefSI)
              {
                unitDefSI = LIBSBML_CPP_NAMESPACE::UnitDefinition::convertToSI(unitDefSI);

                if(NULL != unitDefSI)
                {
                  report(SBMLParserMessage::prtTrace, unitDefSI->getLine()) << "'" << str << "': units '" << LIBSBML_CPP_NAMESPACE::UnitDefinition::printUnits(unitDefSI, true) << "'. Model value : " << result;

                  REAL factor = 1.0;

                  LIBSBML_CPP_NAMESPACE::ListOfUnits * lou = unitDefSI->getListOfUnits();
                  for(UINTEGER ui = 0; ui < lou->size(); ++ui)
                  {
                    LIBSBML_CPP_NAMESPACE::Unit * unit = lou->get(ui);

                    MAP_UNIT_SPEC::const_iterator it = units.find(unit->getKind());
                    if((units.size() > 0)&&(units.end() == it))
                    {
                      STRINGSTREAM ssTemp;
                      std::transform(units.begin(), units.end(), std::ostream_iterator< const char * >(ssTemp, "', '"), SBMLHelper::unitToString);

                      report(SBMLParserMessage::prtWarning, unitDefSI->getLine()) << "'" << str << "': invalid units '" << STRING(LIBSBML_CPP_NAMESPACE::UnitDefinition::printUnits(unitDefSI, true)) << "', can only contain: '" << ssTemp.rdbuf() << "'";
                    }
                    else
                    {
                      factor *= pow(unit->getMultiplier() * pow(10, unit->getScale()), unit->getExponent());
                      if(units.end() != it)
                      {
                        if(unit->getExponent() ==  (*it).second.dExponent)
                          factor *= (*it).second.dMultiplier;
                        else
                          report(SBMLParserMessage::prtWarning, unitDefSI->getLine()) << "'" << str << "': invalid units '" << STRING(LIBSBML_CPP_NAMESPACE::UnitDefinition::printUnits(unitDefSI, true)) << "', expected exponent " << (*it).second.dExponent << " for " << LIBSBML_CPP_NAMESPACE::UnitKind_toString((*it).first);
                      }
                    }
                  }
                  report(SBMLParserMessage::prtTrace, nodeBase->getLine()) << "'" << str << "' has factor = " << factor; 
                  result *= factor;
                }
              }
              else
                report(SBMLParserMessage::prtWarning, nodeBase->getLine()) << "Units associated with '" << str << "' cannot be retrieved.";
            }

            report(SBMLParserMessage::prtTrace, nodeBase->getLine()) << "'" << str << "' = " << result << " " << ((NULL == unitDefSI) ? " " : STRING(LIBSBML_CPP_NAMESPACE::UnitDefinition::printUnits(unitDefSI, true)));

            if(NULL != unitDefSI)
              delete unitDefSI;

            return true;
          }
        }
      }

      report(SBMLParserMessage::prtError, nodeBase->getLine()) << "could not find identifier '"
                << str << "' among species, compartments or local and global parameters.";
      return false;
    };
  };

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_SBML_HELPER_HPP_ */
