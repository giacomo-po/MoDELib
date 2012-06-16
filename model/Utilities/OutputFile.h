/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_OUTPUTFILE_H_
#define model_OUTPUTFILE_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream> 

#include <model/Utilities/StaticID.h>

namespace model {
//	namespace utils {
		
		//template <bool autoDelete = 0> // the problem with this is that then StaticID<Outputfile<outodelete> > has two distinct counters
		class OutputFile : public StaticID<OutputFile>{
			
			///////////////////////////////
			// Private
		private:
			///////////////////////////////
			
			std::ofstream* p_outputfile;
			std::string filename;
			bool autoDelete;
			
			//////////////////////////////////////////////////////////////
			void default_filename(){
				std::ostringstream filestream;
				filestream << "output_" << this->sID << ".txt";
				filename=filestream.str();
			}
			
			//////////////////////////////////////////////////////////////
			void new_Output(){
				p_outputfile = new std::ofstream;
				
				if (autoDelete){
					deleteFile();
				}
				
				p_outputfile->open(filename.c_str(), std::ios::out | std::ios::app);
			}
			
			//////////////////////////////////////////////////////////////
			// Close file and delete if empty
			void close_Output(){
				bool isEmpty = !p_outputfile->tellp();	// if tellp()==0 then file is empty
				
				p_outputfile->close();
				
				if( isEmpty){
					std::cout<<"File is empty. ";
					deleteFile();
				} 
				
			}
			
			//////////////////////////////////////////////////////////////
			// Close file and delete if empty
			void deleteFile(){
				std::cout<<"Deleting file: "<<filename<<std::endl;
				remove(filename.c_str());
			}
			
			///////////////////////////////
			// Public
		public:
			///////////////////////////////
			
			//////////////////////////////////////////////////////////////
			// Constructor
			OutputFile(const bool & autodelete_in=0) : autoDelete(autodelete_in) {
				//			autoDelete=0;
				default_filename(); 
				new_Output();
			}
			
			//		//////////////////////////////////////////////////////////////
			//		// Constructor with filename
			//		OutputFile(const std::string & filename_in) : autoDelete(0) {
			//			filename=filename_in;
			////			autoDelete=0;
			//			new_Output();
			//		}
			
			//////////////////////////////////////////////////////////////
			// Constructor with filename and autodelete
			OutputFile(const std::string & filename_in, const bool & autodelete_in=0)  : autoDelete(autodelete_in) {
				filename=filename_in;
				//			autoDelete=autodelete_in;
				new_Output();
			}
			
			//////////////////////////////////////////////////////////////
			// Copy Constructor
			OutputFile(const OutputFile &) : StaticID<OutputFile>(), autoDelete(0){
				//autoDelete=0;
				default_filename(); 
				new_Output();
			}
			
			//////////////////////////////////////////////////////////////
			// Assigment operator 
			OutputFile& operator=(const OutputFile& Other){
				if (this != &Other){
					autoDelete=0;
					default_filename(); 
					new_Output();
				}
				return *this;
			}
			
			//////////////////////////////////////////////////////////////
			// Destructor
			~OutputFile(){
				close_Output();
				delete p_outputfile;  // is this needed ???????
			}
			
			//		//////////////////////////////////////////////////////////////
			//		// open
			//		void open(const std::string & filename_in){
			//			open(filename_in,0);
			//		}
			
			//////////////////////////////////////////////////////////////
			// open with outodelete
			void open(const std::string & filename_in, const bool & autodelete_in=0){
				close_Output();
				autoDelete=autodelete_in;
				filename=filename_in;
				new_Output();
			}
			
			//////////////////////////////////////////////////////////////
			// define operator << for template types
			template <typename T>
			OutputFile & operator<<(const T & AnyOutput){
				//std::cout<<"AnyOutput"<<std::endl;
				
				(*p_outputfile) << AnyOutput;
				return *this;
			}

			

			
			//////////////////////////////////////////////////////////////
			// define operator << for std::endl
			
			// std::endl is a function that takes and returns a reference to a std::basic_ostream 
			typedef std::basic_ostream<char, std::char_traits<char> > StlEndl_IO;
			// define StlEndl as a pointer-to-function taking and returning a reference to StlEndl_IO
			typedef StlEndl_IO& (*StlEndl)(StlEndl_IO&);
			
			// Overload << for Std::endl
			OutputFile & operator<<(StlEndl manip){
				manip(*p_outputfile);
				return *this;
			}
			
			
		};
		
		//////////////////////////////////////////////////////////////
//	} // namespace utils
} // namespace model
#endif
