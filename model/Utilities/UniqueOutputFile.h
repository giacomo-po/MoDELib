/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_UniqueOutputFile_H_
#define model_UniqueOutputFile_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream> 
#include <assert.h> 

namespace model {

		template <char prefix>
		class UniqueOutputFile {
			
			/*! \brief
			 *  A class template for outputing data to files with automatic sequential filename.
			 *  Instances of UniqueOutputFile can be used analogously to std::cout.
			 *  The first template parameter is the prefix of the filename. The second template
			 *  parameter defines the behaviour of the class when an existing file with the same name is found.
			 *
			 *  Example: create 4 files named A_0, A_1, A_2, A_3 containing int 0,1,2,3 respectively.
			 *  \code
			 *  for (int k=0;k<4;++k){
			 *  model::utils::UniqueOutputFile<'A',1> file;
			 *  file<<k;
			 *  }
			 *  \endcode
			 */
			
			// std::endl is a function that takes and returns a reference to a std::basic_ostream 
			typedef std::basic_ostream<char, std::char_traits<char> > StlEndl_IO;
			// define StlEndl as a pointer-to-function taking and returning a reference to StlEndl_IO
			typedef StlEndl_IO& (*StlEndl)(StlEndl_IO&);

			
			std::ofstream* p_ofstream;
			std::string filename;

			
			//////////////////////////////////////////////////////////////
			// Copy Constructor: private
			UniqueOutputFile(const UniqueOutputFile &){
				new_file();
			}
			
			//////////////////////////////////////////////////////////////
			// Assigment operator: private 
			UniqueOutputFile& operator=(const UniqueOutputFile& Other){
				if (this != &Other){
					new_file();
				}
				return *this;
			}
			
			
			//////////////////////////////////////////////////////////////
			void new_file(){
				
				// define the filename 
				std::ostringstream filestream;
//				switch (useTXT){
//					case true:
						filestream << prefix << "/" << prefix << "_0.txt"; // TXT
//						break;
//					default:
//						filestream << prefix << "/" << prefix << "_" << this->sID << ".bin"; // BIN
//				}
				filename=filestream.str();
				
				// eventually delete existing file with same name	
//				if (autoDelete){
//					deleteFile();
//				}
	
				// open a file named filename
				p_ofstream = new std::ofstream;
//				switch (useTXT){
//					case true:
						p_ofstream-> open(filename.c_str(), std::ios::out | std::ios::app ); // TXT
//						break;
//					default:
//						p_ofstream-> open(filename.c_str(), std::ios::out | std::ios::app | std::ios::binary); // BIN
//				}
				
				
				if (p_ofstream->is_open()){
					//std::cout<<"Creating file: "<<filename<<std::endl;
				}
				else{
					std::cout<<"FILE: "<<filename<<" CANNOT BE OPENED" <<std::endl;
					std::cout<<"TRY CREATING THE FOLDER ./"<<prefix<<" FIRST" <<std::endl;
					assert(0);				
				}
				

			}
			
			
		public:
			UniqueOutputFile()
            {
				new_file();
			}
						
			~UniqueOutputFile()
            {
                delete p_ofstream;
			}
			
			template <typename T>
			UniqueOutputFile & operator<<(const T & AnyOutput)
            {// define operator << for template types
				*p_ofstream << AnyOutput;
				return *this;
			}
			
			
			UniqueOutputFile & operator<<(StlEndl manip)
            {// Overload operator << for std::endl
				manip(*p_ofstream);
				return *this;
			}
			
			
		};
		
} // namespace model
#endif
