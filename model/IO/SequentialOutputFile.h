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

#ifdef _MODEL_MPI_
#include <mpi.h>
#endif

namespace model {
    
    /*! \brief
     *  A class template for outputing data to files with automatic sequential filename.
     *  Instances of SequentialOutputFile can be used analogously to std::cout.
     *  The first template parameter is the prefix of the filename. The second template
     *  parameter defines the behaviour of the class when an existing file with the same name is found.
     *
     *  Example: create 4 files named A_0, A_1, A_2, A_3 containing int 0,1,2,3 respectively.
     *  \code
     for (int k=0;k<4;++k)
     {
        model::SequentialOutputFile<'A',1> file;
        file<<k;
     }
     * \endcode
     */
    template <char prefix, bool autoDelete=true>
    class SequentialOutputFile : public StaticID<SequentialOutputFile<prefix,autoDelete> >,
    //    /*                        */ public std::ofstream {
    /*                        */ private std::ofstream
    {
        

        
        typedef std::basic_ostream<char, std::char_traits<char> > StlEndl_IO;
        // define StlEndl as a pointer-to-function taking and returning a reference to StlEndl_IO
        typedef StlEndl_IO& (*StlEndl)(StlEndl_IO&);

        
        /**********************************************************************/
#ifdef _MODEL_MPI_
        int getRank() const
        {
            // Check that MPI is initialized
            int temp(0);
            MPI_Initialized(&temp);
            assert(temp && "MPI NOT INITIALIZED.");
            
            // Write mpi rank into temp
            MPI_Comm_rank(MPI_COMM_WORLD,&temp);
            return temp;
        }
#endif
        
        /* copy constructor(private) ******************************************/
        SequentialOutputFile(const SequentialOutputFile &){
            assert(0 && "SEQUENTIALOUTPUTFILE CANNOT BE COPIED.");
            //                new_file();
        }
        
        /* assignment operator (private) **************************************/
        SequentialOutputFile& operator=(const SequentialOutputFile& Other)
        {
            if (this != &Other){
                assert(0 && "SEQUENTIALOUTPUTFILE CANNOT BE ASSIGNED.");
            }
            return *this;
        }
        
        /* deleteFile *********************************************************/
        void deleteFile() const
        {
            remove(filename.c_str());
        }
        
        /* getFilename ********************************************************/
        std::string getFilename() const
        {
            std::ostringstream filestream;
            filestream << prefix << "/" << prefix << "_" << this->sID << ".txt"; // TXT
            return filestream.str();
        }
        
        /* openFile ***********************************************************/
        void openFile()
        {
            // possibly delete existing file with same name
            if (autoDelete){
                deleteFile();
            }
            
            this->open(filename.c_str(), std::ios::out | std::ios::app ); // TXT
            
            if (this->is_open()){
                //std::cout<<"Creating file: "<<filename<<std::endl;
            }
            else{
                std::cout<<"FILE: "<<filename<<" CANNOT BE OPENED" <<std::endl;
                std::cout<<"TRY CREATING THE FOLDER ./"<<prefix<<" FIRST." <<std::endl;
                assert(0);
            }
        }
        
    public:

#ifdef _MODEL_MPI_
        const int mpiRank;
#endif

        
        //! the name of the file
        const std::string filename;
        
        
        /* Constructor ********************************************************/
        SequentialOutputFile() :
#ifdef _MODEL_MPI_
        /* init list          */ mpiRank(getRank()),
#endif
        /* init list          */ filename(getFilename())
        {/*
          *
          */
#ifdef _MODEL_MPI_
            if(mpiRank==0)
            {
                openFile();
            }
#else
            openFile();
#endif
            
        }
        
        /* Destructor *********************************************************/
        ~SequentialOutputFile()
        {
#ifdef _MODEL_MPI_
//            int mpiRank;
//            MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
            if(mpiRank==0)
            {
                const bool isEmpty(!this->tellp());	// if tellp()==0 then file is empty
                
                this->close();
                
                if( isEmpty){
                    //    deleteFile();
                }
            }
#else
            const bool isEmpty(!this->tellp());	// if tellp()==0 then file is empty
            
            this->close();
            
            if( isEmpty){
                //    deleteFile();
            }
#endif
        }
        
        /**********************************************************************/
        template <typename T>
        SequentialOutputFile<prefix,autoDelete>& operator<< (const T& os)
        {
#ifdef _MODEL_MPI_
//            int mpiRank;
//            MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
            if(mpiRank==0)
            {
                *static_cast<std::ofstream* const>(this)<<os;
            }
#else
            *static_cast<std::ofstream* const>(this)<<os;

#endif
            
            //            std::ofstream::operator<<(os);
            return *this;
        }
        
        /**********************************************************************/
        SequentialOutputFile<prefix,autoDelete>& operator<<(StlEndl manip)
        {/*! Overload << for Std::endl
          */
#ifdef _MODEL_MPI_
//            int mpiRank;
//            MPI_Comm_rank(MPI_COMM_WORLD,&mpiRank);
            if(mpiRank==0)
            {
                manip(*static_cast<std::ofstream* const>(this));
            }
#else
            manip(*static_cast<std::ofstream* const>(this));
            
#endif            
            return *this;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
} // namespace model
#endif
