/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VERTEXREADER_H_
#define model_VERTEXREADER_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <time.h> // clock()
#include <assert.h>
#include <Eigen/Dense>
#include <model/Utilities/BinaryFileReader.h>

namespace model
{

    /*!	\brief A class template to read Vertex data from file.
     *	\param[c]    the data files are c/c_x.txt or c/c_x.bin
     *	\param[cols] is the number of columns in the file.
     *	\param[scalra] the type of data stored in the columns (e.g. double, int, ...)
     *	First column is the Vertex ID (int), other columns are optional.
     */
    template <char c, int cols, typename scalar>
	class VertexReader : public std::map<int, Eigen::Matrix<scalar,1,cols-1>,
	/*                                */ std::less<int>,
	/*                                */ Eigen::aligned_allocator<std::pair<const int, Eigen::Matrix<scalar,1,cols-1> > > >
    {
        
        int currentFrame;
        bool success;
        
        /**********************************************************************/
		void readBIN(const std::string& filename)
        {
			std::cout<<"Reading: "<<filename;
			double t0(clock());
			typedef std::pair<int, Eigen::Matrix<scalar,1,cols-1> > BinVertexType;
			BinaryFileReader<BinVertexType> rV(filename);
            //			double t1(clock());
			for (size_t k=0;k<rV.size();++k)
            {
				//this->insert(std::make_pair(rV[k].first,rV[k].second));
                const bool inserted=this->insert(std::make_pair(rV[k].first,rV[k].second)).second;
                assert(inserted && "COULD NOT INSERT VERTEX AFTER BINARY READ.");
			}
			std::cout<<" ("<<this->size()<<" vertices) ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
		}
        
        
        /**********************************************************************/
        void readTXT(const std::string& filename)
        {
            scalar temp;
            int key;//, row;//, col;
            std::string line;
            Eigen::Matrix<scalar,1,cols-1> vertexPT=Eigen::Matrix<scalar,1,cols-1>::Zero();
            
            std::ifstream ifs ( filename.c_str() , std::ifstream::in );
            
            if (ifs.is_open()) {
                std::cout<<"Reading: "<<filename;
                double t0(clock());
                int row(0);
                while (std::getline(ifs, line))
                {
                    std::stringstream ss(line);
                    
                    int col(0);
                    
                    while (ss >> temp) {
                        if (col==0) {
                            key=static_cast<int>(temp);
                        }
                        else {
                            vertexPT(col-1)=temp;
                        }
                        col++;
                    }
                    
                    this->insert ( std::pair<int, Eigen::Matrix<scalar,1,cols-1> >(key, vertexPT) );
                    
                    row++;
                }
                
                
                ifs.close();
                
                success=true;
                std::cout<<" ("<<this->size()<<" vertices) ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
            }
            else {
                std::cout<<"Unable to  open: "<<filename<<std::endl;
                success=false;
            }

            
        }
        
    public:
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        /*****************************************/
        /* Constructor */
        VertexReader() : currentFrame(-1), success(false)
        {

        }
        
        /*****************************************/
        static std::string getFilename(const int& frameN, const bool& useTXT){
			std::stringstream filename;
			if(useTXT){
				filename << c << "/" << c << "_" << frameN << ".txt";
			}
			else{
				filename << c << "/" << c << "_" << frameN << ".bin";
			}
			return filename.str();
		}
        
        /*****************************************/
        static bool isGood(const int& frameN, const bool& useTXT)
        {
            std::stringstream filename;
            std::ifstream ifs;
            if(useTXT)
            {
                //filename << "V/V_" << frameN << ".txt";; // TXT
                filename << c << "/" << c << "_" << frameN << ".txt";
                ifs.open( filename.str().c_str() , std::ifstream::in );

            }
            else
            {
                filename << c << "/" << c << "_" << frameN << ".bin";
                //filename << "V/V_" << frameN << ".bin";; // BIN
                ifs.open( filename.str().c_str() , std::ifstream::in );

            }
//            switch (useTXT){
//                case true:
//                    //filename << "V/V_" << frameN << ".txt";; // TXT
//                    filename << c << "/" << c << "_" << frameN << ".txt";
//                    ifs.open( filename.str().c_str() , std::ifstream::in );
//                    break;
//                default:
//                    filename << c << "/" << c << "_" << frameN << ".bin";
//                    //filename << "V/V_" << frameN << ".bin";; // BIN
//                    ifs.open( filename.str().c_str() , std::ifstream::in );
//            }
            return ifs.good();
        }
        
        
		/*****************************************/
		bool read(const int& frameN, const bool& useTXT)
        {
			assert(frameN>=0);
			if (frameN!=currentFrame || this->empty())
            {
				currentFrame=frameN;
				
				// clear the content of the map
				this->clear();
				
				if (useTXT){
					readTXT(getFilename(frameN,true));
				}
				else{
					readBIN(getFilename(frameN,false));
				}
			}
            else
            {
                success=true;
                std::cout<<getFilename(frameN,useTXT)<<" already read."<<std::endl;
            }
			return success;
			
		}
        
    };
    
} // namespace model
#endif
