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

namespace model
{
    
    template <char prefix>
    class UniqueOutputFile
    {
        
        // std::endl is a function that takes and returns a reference to a std::basic_ostream
        typedef std::basic_ostream<char, std::char_traits<char> > StlEndl_IO;
        // define StlEndl as a pointer-to-function taking and returning a reference to StlEndl_IO
        typedef StlEndl_IO& (*StlEndl)(StlEndl_IO&);
        
        
        std::ofstream* p_ofstream;
        std::string filename;
        
        
        UniqueOutputFile(const UniqueOutputFile &)
        {
            new_file();
        }
        
        //////////////////////////////////////////////////////////////
        // Assigment operator: private
        UniqueOutputFile& operator=(const UniqueOutputFile& Other)
        {
            if (this != &Other){
                new_file();
            }
            return *this;
        }
        
        
        //////////////////////////////////////////////////////////////
        void new_file()
        {
            // define the filename
            std::ostringstream filestream;
            
            filestream << prefix << "/" << prefix << "_0.txt"; // TXT
            filename=filestream.str();
            
            // open a file named filename
            p_ofstream = new std::ofstream;
            p_ofstream-> open(filename.c_str(), std::ios::out | std::ios::app ); // TXT
            
            if (!p_ofstream->is_open())
            {
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
