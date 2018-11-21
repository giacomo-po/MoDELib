/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SequentialBinFile_H_
#define model_SequentialBinFile_H_

//#include <boost/filesystem.hpp> // this would help in creating directories but generates link error

#include <iostream>
#include <fstream>
#include <string>
#include <sstream> 
#include <assert.h> 



#include <StaticID.h>

namespace model {
//	namespace utils {	
		
		
		
		template <char prefix, typename OutputType, bool autoDelete=true>
		class SequentialBinFile : public StaticID<SequentialBinFile<prefix,OutputType,autoDelete> >{
			
			/*! \brief
			 *  A class template for outputing data to files with automatic sequential filename.
			 *  Instances of SequentialBinFile can be used analogously to std::cout.
			 *  The first template parameter is the prefix of the filename. The second template
			 *  parameter defines the behaviour of the class when an existing file with the same name is found.
			 *
			 *  Example: create 4 files named A_0, A_1, A_2, A_3 containing int 0,1,2,3 respectively.
			 *  \code
			 
			 for (int k=0;k<4;++k){
				model::utils::SequentialBinFile<'A',1> file;
				file<<k;
			 }
			 
			 * \endcode
			 */
			
			
			///////////////////////////////
			// Private
		private:
			///////////////////////////////
			
			std::ofstream* p_ofstream;
			std::string filename;

			
			//////////////////////////////////////////////////////////////
			// Copy Constructor: private
			SequentialBinFile(const SequentialBinFile &){
				new_file();
			}
			
			//////////////////////////////////////////////////////////////
			// Assigment operator: private 
			SequentialBinFile& operator=(const SequentialBinFile& Other){
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
//						filestream << prefix << "/" << prefix << "_" << this->sID << ".txt"; // TXT
//						break;
//					default:
						filestream << prefix << "/" << prefix << "_" << this->sID << ".bin"; // BIN
//				}
				filename=filestream.str();
				
				// eventually delete existing file with same name	
				if (autoDelete){
					deleteFile();
				}
	
				// open a file named filename
				p_ofstream = new std::ofstream;
//				switch (useTXT){
//					case true:
//						p_ofstream-> open(filename.c_str(), std::ios::out | std::ios::app ); // TXT
//						break;
//					default:
//						p_ofstream-> open(filename.c_str(), std::ios::out | std::ios::app | std::ios::binary); // BIN
				p_ofstream-> open(filename.c_str(), std::ios::out  | std::ios::binary); // BIN

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
			
			
			//////////////////////////////////////////////////////////////
			// Close file and delete if empty
			void deleteFile(){
				//std::cout<<"Deleting file: "<<filename<<std::endl;
				remove(filename.c_str());	// ENABLE
			}
			
		public:
			//////////////////////////////////////////////////////////////
			// Constructor
			SequentialBinFile()
            {
				new_file();
			}
			
			
			//////////////////////////////////////////////////////////////
			// Destructor
			~SequentialBinFile()
            {
				
				p_ofstream->close();
				
//              bool isEmpty = !p_ofstream->tellp();	// if tellp()==0 then file is empty
//				if( isEmpty)
//                {
//				//	std::cout<<"File is empty. ";
//					deleteFile();
//				} 
				
				delete p_ofstream;  // is this needed ???????
			}
			
			
			//////////////////////////////////////////////////////////////
			void write(const OutputType& nextOutput)
            {
				p_ofstream->write((char *) &nextOutput, (sizeof nextOutput));

			}
			
            void operator()(const OutputType& nextOutput)
            {
                write(nextOutput);
            }
			
			
		};
		
} // namespace model
#endif
