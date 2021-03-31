/**
 * @file IO.hpp
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
 * Auxiliary classes for stream-based I/O
 */

#ifndef PSSALIB_UTIL_IO_HPP_
#define PSSALIB_UTIL_IO_HPP_

#include <streambuf> // Standard I/O streams buffer

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/categories.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace pssalib
{
namespace io
{
  /**
   * @class null_streambuf
   * @brief A stream buffer that ignores its input.
   */
  template< typename CharT,
            typename Traits = std::char_traits<CharT> >
  class null_streambuf : public std::basic_streambuf<CharT, Traits>
  {
  /////////////////////////////////////
  // Methods
  private:

    /**
     * @copydoc std::basic_streambuf::overflow(int)
     */
    virtual typename std::basic_streambuf<CharT, Traits>::int_type overflow(typename std::basic_streambuf<CharT, Traits>::int_type c)
    {
      return c;
    }
  };

  /**
   * @class onoff_tee_filter
   * @brief A filtering stream that conditionally forwards its input
   */
  template<typename Device>
  class onoff_tee_filter : public boost::iostreams::tee_filter<Device>
  {
  protected:
    //! On/Off callback
    bool (*cbOnCheck) (void *);

    //! On/Off callback data
    void *cbOnCheckData;

  public:
    //! Constructor
    onoff_tee_filter(Device &dev) :
      boost::iostreams::tee_filter<Device>(dev),
      cbOnCheck(NULL),
      cbOnCheckData(NULL)
    {
      // Do nothing
    }

    //! Constructor
    onoff_tee_filter(Device &dev, bool (*cbOn)(void*), void *cbOnData) :
      boost::iostreams::tee_filter<Device>(dev),
      cbOnCheck(cbOn),
      cbOnCheckData(cbOnData)
    {
      // Do nothing
    }

    //! Copy constructor
    onoff_tee_filter(const onoff_tee_filter &right) :
      boost::iostreams::tee_filter<Device>(right),
      cbOnCheck(right.cbOnCheck),
      cbOnCheckData(right.cbOnCheckData)
    {
      // Do nothing
    }

    // Destructor
    virtual ~onoff_tee_filter()
    {
      // Do nothing
    }

    //! read input data
    template<typename Source>
    std::streamsize read(Source& src, 
      typename boost::iostreams::tee_filter<Device>::char_type* s, std::streamsize n)
    {
      std::streamsize result = boost::iostreams::read(src, s, n);
      if (result != -1) {
        bool passOn = true;
        if(NULL != cbOnCheck) passOn = (*cbOnCheck)(cbOnCheckData);
        if(passOn)
        {
          std::streamsize result2 = boost::iostreams::write(this->component(), s, result);
          (void) result2; // Suppress 'unused variable' warning.
          BOOST_ASSERT(result == result2);
        }
      }
      return result;
    }

    //! write data to output
    template<typename Sink>
    std::streamsize write(Sink& snk, 
      const typename boost::iostreams::tee_filter<Device>::char_type* s, std::streamsize n) 
    {
      std::streamsize result = boost::iostreams::write(snk, s, n);
      bool passOn = true;
      if(NULL != cbOnCheck) passOn = (*cbOnCheck)(cbOnCheckData);
      if(passOn)
      {
        std::streamsize result2 = boost::iostreams::write(this->component(), s, result);
        (void) result2; // Suppress 'unused variable' warning.
        BOOST_ASSERT(result == result2);
      }
      return result;
    }
  };

  /**
   * @class onoff_tee_filter
   * @brief A filtering stream that conditionally forwards its input
   */
  // line_num_filter is a model of the Boost concept OutputFilter which
  // inserts a sequential line number at the beginning of every line.
  class prefix_filter
      : public boost::iostreams::output_filter
  {
  public:
      prefix_filter(const char *p) :
        m_isAtBeginningOfLine(false)
      {
        if(NULL != p)
          m_prefix.assign(p);
        else
          m_prefix.clear();
      };

      /**
       * Put a single character into the output buffer.
       * put() must return true if c was written to snk, or false if not.
       * After returning false, put() with the same c might be tried again later.
       */
      template<typename Sink>
      bool put(Sink& snk, char c)
      {
        // If at the start of a line, prepend the prefix to the buffer.
        if (m_isAtBeginningOfLine)
        {
          if(boost::iostreams::write(snk, m_prefix.data(), m_prefix.size()) != m_prefix.size())
            return false;
//           std::size_t len = m_prefix.length();
//           const char * str= m_prefix.c_str();
//           for (std::size_t i = 0; i < len; ++i)
//           {
//             if (!boost::iostreams::put(snk, str[i]))
//               return false;
//           }
          m_isAtBeginningOfLine = false;
        }

        // If the character copied was a newline, get ready for the next line.
        if ('\n' == c)
            m_isAtBeginningOfLine = true;

        // Copy the actual character of data.
        return boost::iostreams::put(snk, c);
      };

      template<typename Device>
      void close(Device&)
      {
        m_isAtBeginningOfLine = true;
      };

  private:
      std::string m_prefix;
      bool        m_isAtBeginningOfLine;
  };

} } // close io & pssalib namespaces

#endif /* PSSALIB_UTIL_IO_HPP_ */
