/**
 * @file JaggedMatrix.hpp
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
 * Declares a templatized container representing a vector of vectors, all sharing
 * the same data type
 */

#ifndef PSSALIB_DATAMODEL_DETAIL_JAGGEDMATRIX_HPP_
#define PSSALIB_DATAMODEL_DETAIL_JAGGEDMATRIX_HPP_

#include "../../typedefs.h"

namespace pssalib
{
namespace datamodel
{
namespace detail
{
  /**
   * @class JaggedMatrix
   * @brief A templetized matrix with variable row length.
   *
   */
  template< typename A >
  class JaggedMatrix
  {
  ////////////////////////////////
  // Attributes
  protected:
    std::size_t uRows, //!< number of rows in the matrix
                uInc;  //!< capacity increment
    std::size_t *uCols,//!< number of columns in the matrix
                *uCols_alloc;//!< @internal storage capacity

    //! Data
    A** data;

  /////////////////////////////////////
  // Constructors
  public:
    //! Default constructor
    JaggedMatrix<A> () :
      uRows(0),
      uInc(0),
      uCols(NULL),
      uCols_alloc(NULL),
      data(NULL)
    {
      // Do nothing
    };

    //! Copy constructor
    JaggedMatrix<A> (JaggedMatrix<A> & other) :
      uRows(0),
      uInc(0),
      uCols(NULL),
      uCols_alloc(NULL),
      data(NULL)
    {
      copy(other);
    };

    /**
     * Creates an empty varying row length matrix, 
     * pre-allocating uC[i] elements for each row i,
     * i = 0, ... , uM.
     * 
     * @param uM Number of rows
     * @param uC number of columns per row
     */
    //! Creates a varying row length(C[i]) matrix
    JaggedMatrix<A> (std::size_t uM, std::size_t* uC) :
      uRows(uM),
      uInc(0),
      uCols(NULL),
      uCols_alloc(NULL),
      data(NULL)
    {
      resize(uM, uC);
    };

    /**
     * Creates an empty varying row length matrix, 
     * pre-allocating uC elements for each row.
     * 
     * @param uM Number of rows
     * @param uC number of columns per row
     */
    JaggedMatrix<A> (std::size_t uM, std::size_t uC) :
      uRows(uM),
      uInc(0),
      uCols(NULL),
      uCols_alloc(NULL),
      data(NULL)
    {
      reserve(uM, uC);
    };

    //! Destructor
    ~JaggedMatrix<A> ()
    {
      free();
    }

  /////////////////////////////////////
  // Methods
  public:
    /**
     * Assignment operator for objects of the same type
     * 
     * @param rhs right-hand side.
     * @return a deep copy of the rhs.
     */
    inline JaggedMatrix<A> &operator=(JaggedMatrix<A> & rhs)
    {
      copy(rhs);
      return (*this);
    };

    /**
     * Get matrix element (const reference).
     * 
     * @param i row index.
     * @param j column index.
     */
    inline const A & operator()(std::size_t i, std::size_t j) const
    {
      return const_cast<const A &>(
        const_cast<JaggedMatrix<A> *>(this)->operator()(i,j));
    };

    /**
     * Get matrix element.
     * 
     * @param i row index.
     * @param j column index.
     */
    inline A & operator()(std::size_t i, std::size_t j)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if((i >= uRows)&&(j >= uCols[i]))
        throw std::runtime_error("JaggedMatrix<A>::operator() - subscript out of range.");
      else
#endif
        return (data[i])[j];
    };

//     /**
//      * Get matrix element (const reference).
//      * 
//      * @param i row index.
//      * @param j column index.
//      */
//     inline const A & operator[](std::size_t i, std::size_t j) const
// #ifndef PSSALIB_NO_BOUNDS_CHECKS
//       throw(std::runtime_error)
// #endif
//     {
//       return const_cast<const A &>(
//         const_cast<JaggedMatrix<A> *>(this)->operator[](i,j));
//     };
// 
//     /**
//      * Get matrix element.
//      * 
//      * @param i row index.
//      * @param j column index.
//      */
//     inline A & operator[](std::size_t i, std::size_t j)
// #ifndef PSSALIB_NO_BOUNDS_CHECKS
//       throw(std::runtime_error)
// #endif
//     {
// #ifndef PSSALIB_NO_BOUNDS_CHECKS
//       if((i < uRows)&&(j < uCols[i]))
// #endif
//         return (data[i])[j];
// #ifndef PSSALIB_NO_BOUNDS_CHECKS
//       else
//         throw std::runtime_error("JaggedMatrix<A>::operator[] - subscript out of range.");
// #endif
//     };

    /**
     * Creates a deep copy of another class instance into this object.
     * 
     * @param other the class instance to copy from.
     */
    void copy(JaggedMatrix<A> & other)
    {
      free(); // clean up

      A           ** temp_data;
      std::size_t *  temp_uCols,
                  *  temp_uCols_alloc,
                      temp_uRows;

      temp_uRows = other.uRows;
      temp_data = new A*[uRows];

      temp_uCols = new std::size_t[temp_uRows];
      temp_uCols_alloc = new std::size_t[temp_uRows];

      memcpy(temp_uCols, other.uCols, uRows*sizeof(std::size_t));
      memcpy(temp_uCols_alloc, other.uCols_alloc, uRows*sizeof(std::size_t));

      // allocate memory & compute the total number of columns
      std::size_t inc = 0;
      for(std::size_t  i = 0; i < uRows; i++)
      {
        if(0 != temp_uCols[i])
        {
          temp_data[i] = new A[uCols_alloc[i]];

          memcpy(temp_data[i], other.data[i], temp_uCols[i]*sizeof(std::size_t));
        }
        else
          temp_data[i] = NULL;
        inc += temp_uCols[i];
      }

      free();

      uRows = temp_uRows;
      uCols = temp_uCols;
      uCols_alloc = temp_uCols_alloc;
      data = temp_data;

      // Update the growth factor
      uInc = std::max(inc / uRows, std::size_t(1));
    }

    /**
     * Appends a new element at the end of a given row.
     * 
     * @param row row index.
     * @param elem reference to the new element.
     */
    void push_back(std::size_t i, A &elem)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(i >= uRows)
        throw std::runtime_error("JaggedMatrix<A>::push_back() - subscript out of range.");
      else
      {
#endif
        std::size_t uLen = uCols[i];
        if(uLen < uCols_alloc[i])
        {
          // simply add it
          (data[i])[uLen] = elem;
          uCols[i]++;
        }
        else
        {
          // need to allocate more memory
          A *temp = data[i];

          uCols_alloc[i] += uInc;
          data[i] = new A[uCols_alloc[i]];
          if(uLen > 0)
          {
            memcpy(data[i], temp, uLen*sizeof(A));
            delete [] temp;
          }
          (data[i])[uCols[i]++] = elem;
        }
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      }
#endif
    }

    /**
     * Preallocate memory for the matrix with variable row lengths.
     * 
     * @param uM number of rows.
     * @param uC vector of variable row length.
     */
    void resize(const std::size_t uM, const std::size_t *uC)
    {
      if((0 == uM)||(NULL == uC))
        throw std::runtime_error("JaggedMatrix<A>::resize() - invalid arguments.");

      if(uM == uRows)
      {
        for(std::size_t i = 0; i < uRows; i++)
        {
          if(uC[i] != uCols_alloc[i])
          {
            A * temp = data[i];
            if(uC[i] > 0)
            {
              data[i] = new A[uC[i]];
              memcpy(data[i], temp, std::min(uC[i],uCols[i])*sizeof(A));
            }
            delete [] temp;

            uCols_alloc[i] = uC[i];
          }

          uCols[i] = uC[i];
        }
//         if(uC == uCols_alloc)
//         {
//           memcpy(uCols, uCols_alloc, uRows*sizeof(std::size_t));
//           return;
//         }
//         else if(uC == uCols)
//         {
//           for(std::size_t i = 0; i < uRows; i++)
//           {
//             if(uCols[i] != uCols_alloc[i])
//             {
//               A * temp = data[i];
//               data[i] = new A[uCols[i]];
//               memcpy(data[i], temp, uCols[i]*sizeof(A));
//
//               delete [] temp;
//             }
//           }
//           return;
//         }
      }
      else
      {
        A           ** temp_data;
        std::size_t *  temp_uCols,
                    *  temp_uCols_alloc,
                       temp_uRows;

        temp_uRows = uM;
        temp_data = new A*[uM];

        temp_uCols = new std::size_t[temp_uRows]; 
        temp_uCols_alloc = new std::size_t[temp_uRows];

        memcpy(temp_uCols, uC, temp_uRows*sizeof(std::size_t));
        memcpy(temp_uCols_alloc, temp_uCols, temp_uRows*sizeof(std::size_t));

        // allocate memory & compute the total number of columns
        std::size_t inc = 0;
        for(std::size_t  i = 0; i < temp_uRows; i++)
        {
          if(0 != temp_uCols[i])
          {
            temp_data[i] = new A[temp_uCols[i]];

            if((NULL != uCols)&&(i < uRows)&&(0 != uCols[i]))
            {
              if(temp_uCols[i] < uCols[i])
                memcpy(temp_data[i], data[i], temp_uCols[i]*sizeof(std::size_t));
              else
                memcpy(temp_data[i], data[i], uCols[i]*sizeof(std::size_t));
            }
          }
          else
            temp_data[i] = NULL;
          inc += temp_uCols[i];
        }

        // clean up original data
        free();

        // assign new values
        uRows = temp_uRows;
        data = temp_data;

        uCols = temp_uCols;
        uCols_alloc = temp_uCols_alloc;

        // Calculate the growth factor
        uInc = inc / temp_uRows;
        if(0 == uInc) uInc = 1;
      }
    };

    /**
     * Resize the matrix assuming constant row lengths.
     * 
     * @param uR number of rows.
     * @param uC number of columns in each row.
     */
    void resize(std::size_t uM, std::size_t uC)
    {
      if((0 == uM)||(0 == uC))
        throw std::runtime_error("JaggedMatrix<A>::resize() - invalid arguments.");

      free();

      uRows = uM;
      data = new A*[uM];

      uCols = new std::size_t[uM]; uCols_alloc = new std::size_t[uM];
      for(std::size_t i = 0; i < uM; i++)
        uCols[i] = uC;
      memcpy(uCols_alloc, uCols, uM*sizeof(std::size_t));

      //
      // Calculate the total number of columns
      for(std::size_t i = 0; i < uM; i++)
        data[i] = new A[uC];

      // Calculate the growth factor
      uInc = uC;
    };

    /**
     * Allocate memory for the matrix assuming constant row lengths.
     * 
     * @param uM number of rows.
     * @param uN number of columns in each row.
     */
    void reserve(std::size_t uM, std::size_t uC)
    {
      if(0 == uC) uC = 1;
      uInc = uC;

      if(data != NULL)
      {
        // Check if more rows are requested
        if(uRows < uM)
        {
          A            **temp_data = data;
          std::size_t *temp_uCols = uCols, *temp_uCols_alloc = uCols_alloc;

          // Adjust the container
          data = new A*[uM];
          // Adjust & fill helper variables
          uCols = new std::size_t[uM]; uCols_alloc = new std::size_t[uM];

          // Nullify the new arrays
          memset(uCols, (unsigned char)0, uM*sizeof(std::size_t));
          memset(uCols_alloc, (unsigned char)0, uM*sizeof(std::size_t));

          // Copy old values
          memcpy(uCols, temp_uCols, uRows*sizeof(std::size_t));
          memcpy(uCols_alloc, temp_uCols_alloc, uRows*sizeof(std::size_t));

          // Copy the data
          for(std::size_t  i = 0; i < uRows; i++)
          {
            data[i] = temp_data[i];
          }

          uRows = uM;

          // Delete temporary variables
          delete [] temp_data;
          delete [] temp_uCols;
          delete [] temp_uCols_alloc;
        }
        else if(uRows > uM)
          throw std::runtime_error("JaggedMatrix<A>::reserve() - invalid arguments.");

        std::size_t uLen;
        A *         temp;
        for(std::size_t  i = 0; i < uRows; i++)
        {
          uLen = uCols_alloc[i];
          if(uLen < uC)
          {
            temp = data[i];
            data[i] = new A[uC];
            if(uLen > 0)
            {
              memcpy(data[i], temp, uLen*sizeof(A));
              delete [] temp;
            }
            uCols_alloc[i] = uC;
          }
        }
      }
      else
      {
        uRows = uM;
        data = new A*[uRows];
        uCols = new std::size_t[uRows]; uCols_alloc = new std::size_t[uRows];
        for(std::size_t  i = 0; i < uRows; i++)
        {
          data[i]       = new A[uC];
          uCols[i]      = 0;
          uCols_alloc[i]= uC;
        }
      }
    };

    /**
     * Clear data structures.
     */
    inline void clear()
    {
      for(std::size_t i = 0; i < uRows; i++)
      {
        for(std::size_t j = 0; j < uCols[i]; j++)
          (&(data[i])[j])->A::~A();

        memset(data[i], (unsigned char)0, uCols[i]*sizeof(A));
      }
      memset(uCols, (unsigned char)0, uRows*sizeof(std::size_t));
    };

    /**
     * Compact the memory by free'ing the 
     * alloocated memory surplus.
     */
    inline void compact()
    {
      resize(uRows, uCols);
    }

    /**
     * Free allocated resources.
     */
    inline void free()
    {
      if(data != NULL)
      {
        for(std::size_t i = 0; i < uRows; i++)
          if(0 != uCols_alloc[i])
            delete [] (data[i]);

        uRows = 0;  uInc = 0;

        delete [] data;
        data = NULL;

        delete [] uCols;
        uCols = NULL;

        delete [] uCols_alloc;
        uCols_alloc = NULL;
      }
    }

    /**
     * Get number of rows.
     * 
     * @return number of rows.
     */
    inline std::size_t get_rows() const
    {
      return uRows;
    };

    /**
     * Get number of columns for a given row.
     * 
     * @param row row index.
     * @return number of columns for row @row .
     */
    inline std::size_t get_cols(std::size_t row) const
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if(row >= uRows)
        throw std::runtime_error("JaggedMatrix<A>::get_cols() - subscript out of range.");
      else
#endif
        return uCols[row];
    };

    /**
     * Get number of columns per each rows
     * 
     * @return Vector containing number of columns per each row.
     */
    inline std::size_t * get_cols() const
    {
      return uCols;
    };

    /**
     * Sort elements in each row
     */
    inline void sort_cols()
    {
      for(std::size_t i = 0; i < uRows; i++)
      {
        if(uCols[i]>0)
        {
          std::sort(data[i],data[i]+uCols[i]);
        }
      }
    }

    /**
     * Swap rows
     * 
     * @param i1 index of the first row
     * @param i2 index of the second row
     */
    void swap(std::size_t i1, std::size_t i2)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if((i1 >= uRows)&&(i2 >= uRows))
        throw std::runtime_error("JaggedMatrix<A>::swap() - subscript out of range.");
      else
#endif
        std::swap(data[i1], data[i2]);
    }

    /**
     * Swap data within the same row
     * 
     * @param i index of the row
     * @param j1 index of the first column
     * @param j2 index of the second column
     */
    void swap(std::size_t i, std::size_t j1, std::size_t j2)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if((i >= uRows)&&(std::max(j1,j2) >= uCols[i]))
        throw std::runtime_error("JaggedMatrix<A>::swap() - subscript out of range.");
      else
#endif
        std::swap(data[i][j1], data[i][j2]);
    }

    /**
     * Swap data between rows
     * 
     * @param i1 index of the first row
     * @param j1 index of the first column
     * @param i2 index of the second row
     * @param j2 index of the second column
     */
    void swap(std::size_t i1, std::size_t j1, std::size_t i2, std::size_t j2)
    {
#ifndef PSSALIB_NO_BOUNDS_CHECKS
      if((i1 >= uRows)&&(i2 >= uRows)&&(j1 >= uCols[i1])&&(j2 >= uCols[i2]))
        throw std::runtime_error("JaggedMatrix<A>::swap() - subscript out of range.");
      else
#endif
        std::swap(data[i1][j1], data[i2][j2]);
    }

    /**
     * Shift operator for console output
     * 
     * @param output an instance of @link std::ostream
     * @param M an instance of @link JaggedMatrix<A>
     */
    friend std::ostream & operator<<(std::ostream & output, 
                                     const JaggedMatrix<A> & M)
    {
      output << "Array [" << M.uRows << 'x' 
        << "*]" << std::endl << std::setprecision(7);
      for(std::size_t i = 0; i < M.uRows; i++)
      {
        for(std::size_t j = 0; j < M.uCols[i]; j++)
        {
          output << std::setw(9) << M(i,j) << ' ';
        }
        output << std::endl;
      }
      return output;
    };
  };

} } } // close namespaces detail, datamodel & pssalib

#endif /* PSSALIB_DATAMODEL_DETAIL_JAGGEDMATRIX_HPP_ */
