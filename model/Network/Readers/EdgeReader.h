/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EDGEREADER_H_
#define model_EDGEREADER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h> // clock()
#include <map>
#include <assert.h>
#include <utility>
#include <Eigen/Dense>
#include <model/Utilities/BinaryFileReader.h>

namespace model
{
	
    /*!	\brief A class template to read Edge data from file.
     *	\param[cols] is the number of columns in the file.
     *	First column is the source Vertex ID, second column is the sink Vertex ID.
     *  Other columns are optional.
     */
	template <char c, int cols, typename scalar>
	class EdgeReader : public std::map<std::pair<int, int>, Eigen::Matrix<scalar,1,cols-2>,
	/*                              */ std::less<std::pair<int, int> >,   
	/*                              */ Eigen::aligned_allocator<std::pair<const std::pair<int, int>, Eigen::Matrix<scalar,1,cols-2> > > > {	
		
		int currentFrame;
		
        /**********************************************************************/
		bool readTXT(const std::string& filename)
        {
			scalar temp;
			//int key1, key2, row, col;
			std::string line;
			
			Eigen::Matrix<scalar,1,cols-2> edgeNB = Eigen::Matrix<scalar,1,cols-2>::Zero();
			
			std::ifstream ifs ( filename.c_str() , std::ifstream::in );
			
            bool success(false);
            
			if (ifs.is_open())
            {
				std::cout<<"Reading: "<<filename;
				double t0(clock());
				int row = 0;
				while (std::getline(ifs, line))
                {
					std::stringstream ss(line);
					
					int col=0;
					int key1, key2;
					while (ss >> temp)
                    {
						
						switch (col)
                        {
							case 0:
								key1=static_cast<int>(temp);
								break;
							case 1:
								key2=static_cast<int>(temp);
								break;
								
							default:
								edgeNB(col-2)=temp;
								break;
						}
						
						col++;
					}
					
					this->insert ( std::pair<std::pair<int, int>, Eigen::Matrix<scalar,1,cols-2> >(std::pair<int, int>(key1, key2), edgeNB) );
					
					row++;
				}
				
				
				ifs.close();
				
				success=true;
				
                std::cout<<" ("<<this->size()<<" edges) ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
			}
			else
            {
				std::cout<<"Unable to  open:"<<filename<<std::endl;
			}
			
            return success;
		}
		
        /**********************************************************************/
		bool readBIN(const std::string& filename)
        {/*! Reads the binary file filename and stores its data in this
          */
            //bool success(false);
			std::cout<<"Reading: "<<filename;
			double t0(clock());
			typedef std::pair<std::pair<int,int>, Eigen::Matrix<scalar,1,cols-2> > BinEdgeType;
			BinaryFileReader<BinEdgeType> rE(filename);
			for (unsigned int k=0;k<rE.size();++k)
            {
//				this->insert(std::make_pair(rE[k].first,rE[k].second));
                const bool inserted=this->insert(std::make_pair(rE[k].first,rE[k].second)).second;
				assert(inserted && "COULD NOT INSERT EDGE AFTER BINARY READ.");
			}
			std::cout<<" ("<<this->size()<<" edges) ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
            return rE.success();
		}
		
	public:	
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		
        /**********************************************************************/
		EdgeReader () :
        /* init list */ currentFrame(-1)
 //       /* init list */ success(false)
        {/*! Constructor initializes currentFrame to -1 so that the statement
          *  frameN!=currentFrame will initially return false
          */
        }
		
        /**********************************************************************/
		static std::string getFilename(const int& frameN, const bool& useTXT)
        {
			std::stringstream filename;
			if(useTXT){
				filename << c << "/" << c << "_" << frameN << ".txt";
			}
			else{
				filename << c << "/" << c << "_" << frameN << ".bin";
			}
			return filename.str();
		}
		
        /**********************************************************************/
		static bool isGood(const int& frameN, const bool& useTXT){
			/*!	Checks whether the file named "E/E_ \param[frameN] .txt" is good for reading.
			 */
			std::ifstream ifs ( getFilename(frameN,useTXT).c_str() , std::ifstream::in );
			return ifs.good();
		}
		
		/**********************************************************************/
		bool read(const int& frameN, const bool& useTXT)
        {/*! @param[in] frameN the frame number to be read
          *  @param[in] useTXT if true text files (.txt) are read, 
          *  otherwise binary files (.bin) are read
          *  \returns true if the file has been read correctly
          */
			assert(frameN>=0 && "frameN MUST BE >= 0");
			bool success(false);
            if (frameN!=currentFrame || this->empty())
            {
				currentFrame=frameN;
				
				// clear the content of the map
				this->clear();
				
				if (useTXT)
                {
					success=readTXT(getFilename(frameN,true));
				}
				else
                {
					success=readBIN(getFilename(frameN,false));
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
	
	/******************************************************************/
	/******************************************************************/
} /* namespace model */
#endif
