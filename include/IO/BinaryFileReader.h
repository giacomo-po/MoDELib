/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez <ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BINARYFILEREADER_H_
#define model_BINARYFILEREADER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <vector>
#include <cmath>

namespace model
{
    
    /*!	\brief A class template that reads generic binary files containing arrays
     *  of DataType elements.
     *
     * remember that there are 8 bytes in a double so, for example, if
     * double fnum[4] = {9.5, -3.4, 1.0, 2.1};
     * then sizeof(fnum) is 32
     */
    template <typename DataType>
    class BinaryFileReader : public std::vector<DataType>
    {
        
        static constexpr size_t _bytes_in_Type=sizeof (DataType);

        
        std::ifstream::pos_type _bytes_in_file;
        bool _success;
        
    public:
        
        /**********************************************************************/
        BinaryFileReader(const std::string& filename) :
        /* init list */ _bytes_in_file(0)
        /* init list */,_success(false)
        {/*! @param[in] filename The name of the file to be read
          *  The constructor reads the file filename and stores DataType elements
          *  in its internal array.
          */
            
//            this->clear();
            
            // open file with the ios::ate flag in order to set the get pointer at the end of the file.
            std::ifstream file (filename.c_str(), std::ios::in|std::ios::binary|std::ios::ate);
            
            if (file.is_open())
            {
                
                // Since the get pointer is at the end of the file, call file.tellg() to obtain the number of bytes in the file
                _bytes_in_file = file.tellg();
                
                // The number of bytes in the file must be a multiple of the number of bytes in Type
                assert((_bytes_in_file%_bytes_in_Type)==0 && " Incorrect format.");
                
                // Calculate how many element of type Type are in the file
                const size_t _array_size = _bytes_in_file/_bytes_in_Type;
                this->resize(_array_size);
                
                // resize memblock with _array_size
                
                // bring back the get pointer to the beginning of the file
                file.seekg (0, std::ios::beg);
                
                // read the file into memblock
                file.read (reinterpret_cast<char*>(this->data()), _bytes_in_file);
                
                // close the file
                file.close();
                
                _success=true;
            }
            else
            {
                std::cout << "Unable to open file "<<filename<<std::endl;
            }
            
        }
        
        /**********************************************************************/
        const size_t& bytes() const
        {/*! The number of bytes in the file read
          */
            return _bytes_in_file;
        }
        
        /**********************************************************************/
        const bool& success() const
        {/*! The number of bytes in the file read
          */
            return _success;
        }

    };
    
} // namespace model
#endif

//        /**********************************************************************/
//        ~BinaryFileReader()
//        {/*! Destructor
//          */
//            //		delete[] memblock;
//            		delete memblock;
//        }

//        /**********************************************************************/
//        const size_t& size() const
//        {/*! The number of DataType(s) elements read and stored in this
//          */
//            return _array_size;
//        }


//        /**********************************************************************/
//        const DataType& operator[](const size_t& k) const
//        {/*! @param[in] k the position in the BinaryFileReader container
//          *  \returns the k-th DataType element in the BinaryFileReader container
//          */
//            assert( k<_array_size && "Index out of bound.");
//            return memblock[k];
//        }
