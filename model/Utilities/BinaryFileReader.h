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

namespace model {

template <typename Type>
class BinaryFileReader  {
	
	/*!	\brief ??????
	 */
	
	Type* memblock;
	std::ifstream::pos_type _bytes_in_file;
	size_t _array_size;
	const size_t _bytes_in_Type;
	
	bool success;
	
public:
	
	// remember that there are 8 bytes in a double so, for example, if 
	// double fnum[4] = {9.5, -3.4, 1.0, 2.1};
	// then sizeof(fnum) is 32
	

	BinaryFileReader(const std::string& filename) : _bytes_in_Type(sizeof (Type) ){
		
		// open file with the ios::ate flag in order to set the get pointer at the end of the file. 
		std::ifstream file (filename.c_str(), std::ios::in|std::ios::binary|std::ios::ate);
		
		if (file.is_open()){
		
			// Since the get pointer is at the end of the file, call file.tellg() to obtain the number of bytes in the file
			_bytes_in_file = file.tellg();
//			std::cout << "File '"<< filename <<"' contains:"<<std::endl;
//			std::cout << "	bytes in file : "     << _bytes_in_file <<std::endl;
//			std::cout << "	bytes per element:  " << _bytes_in_Type  <<std::endl;
			
			// The number of bytes in the file must be a multiple of the number of bytes in Type 
			assert((_bytes_in_file%_bytes_in_Type)==0 && " Incorrect format.");
			
			// Calculate how many element of type Type are in the file
			_array_size = _bytes_in_file/_bytes_in_Type;
//			std::cout << "	array size:  "<< _array_size  <<std::endl;
			
			// resize memblock with _array_size
			memblock = new Type [_array_size];
			
			// bring back the get pointer to the beginning of the file
			file.seekg (0, std::ios::beg);
			
			// read the file into memblock
			file.read (reinterpret_cast<char*>(memblock), _bytes_in_file);
			
			// close the file
			file.close();
			
			success=true;
		}
		else{
			success=false;
			std::cout << "Unable to open file "<<filename<<std::endl;
		}

	}
	
	~BinaryFileReader(){
//		delete[] memblock;
//		delete memblock;
	}
	
	
	const size_t& size() const{
		return _array_size;
	}
	
	
	const Type& operator[](const size_t& k) const{
		assert( k<_array_size && "Index out of bound.");
		return memblock[k];
	}
	
	
	
};

} // namespace model
#endif
