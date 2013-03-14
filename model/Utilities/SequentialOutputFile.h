/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SequentialOutputFile_H_
#define model_SequentialOutputFile_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream> 
#include <assert.h> 
#include <model/Utilities/StaticID.h>

namespace model {
    
    template <char prefix, bool autoDelete=true, bool useTXT=true>
    class SequentialOutputFile : public StaticID<SequentialOutputFile<prefix,autoDelete,useTXT> >,
    /*                        */ public std::ofstream {
        
        /*! \brief
         *  A class template for outputing data to files with automatic sequential filename.
         *  Instances of SequentialOutputFile can be used analogously to std::cout.
         *  The first template parameter is the prefix of the filename. The second template
         *  parameter defines the behaviour of the class when an existing file with the same name is found.
         *
         *  Example: create 4 files named A_0, A_1, A_2, A_3 containing int 0,1,2,3 respectively.
         *  \code
         
         for (int k=0;k<4;++k){
         model::utils::SequentialOutputFile<'A',1> file;
         file<<k;
         }
         
         * \endcode
         */
        
    private:
        
        /* copy constructor(private) ******************************************/
        SequentialOutputFile(const SequentialOutputFile &){
            assert(0 && "SEQUENTIALOUTPUTFILE CANNOT BE COPIED.");
            //                new_file();
        }
        
        /* assignment operator (private) **************************************/
        SequentialOutputFile& operator=(const SequentialOutputFile& Other){
            if (this != &Other){
                assert(0 && "SEQUENTIALOUTPUTFILE CANNOT BE ASSIGNED.");
            }
            return *this;
        }
        
        /* deleteFile *********************************************************/
        void deleteFile() const {
            remove(filename.c_str());
        }
        
        /* getFilename ********************************************************/
        std::string getFilename() const {
            std::ostringstream filestream;
            switch (useTXT){
                case true:
                    filestream << prefix << "/" << prefix << "_" << this->sID << ".txt"; // TXT
                    break;
                default:
                    filestream << prefix << "/" << prefix << "_" << this->sID << ".bin"; // BIN
            }
            return filestream.str();
        }
        
    public:
        
        //! the name of the file
        const std::string filename;
        
        /* Constructor ********************************************************/
        SequentialOutputFile() : filename(getFilename()) {
            
            // possibly delete existing file with same name	
            if (autoDelete){
                deleteFile();
            }
            
            // open a file named filename
            switch (useTXT){
                case true:
                    this->open(filename.c_str(), std::ios::out | std::ios::app ); // TXT
                    break;
                default:
                    this->open(filename.c_str(), std::ios::out | std::ios::app | std::ios::binary); // BIN
            }
            
            
            if (this->is_open()){
                //std::cout<<"Creating file: "<<filename<<std::endl;
            }
            else{
                std::cout<<"FILE: "<<filename<<" CANNOT BE OPENED" <<std::endl;
                std::cout<<"TRY CREATING THE FOLDER ./"<<prefix<<" FIRST" <<std::endl;
                assert(0);				
            }
            
        }
        
        /* Destructor *********************************************************/
        ~SequentialOutputFile(){
            const bool isEmpty(!this->tellp());	// if tellp()==0 then file is empty
            
            this->close();
            
            if( isEmpty){
            //    deleteFile();
            } 
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
} // namespace model
#endif
