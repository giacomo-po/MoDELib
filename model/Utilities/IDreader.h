/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_IDREADER_H_
#define model_IDREADER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h> // clock()
#include <map>
#include <assert.h>
#include <utility>
//#include <Eigen/Dense>
#include <model/Utilities/BinaryFileReader.h>

namespace model
{
	
    /*!	\brief A class template
     */
	template <char c, int keySize, int valueSize, typename T>
	class IDreader : public std::map<std::array<int, keySize>, std::array<T, valueSize> >
    {
		
        typedef std::array<int, keySize> KeyType;
        typedef std::array<T, valueSize> ValueType;
        
		int currentFrame;
		
        /**********************************************************************/
		bool readTXT(const std::string& filename)
        {
            
			std::ifstream ifs ( filename.c_str() , std::ifstream::in );
			
            bool success(false);
            
			if (ifs.is_open())
            {
				std::cout<<"Reading: "<<filename;
				double t0(clock());
                std::string line;


				int row = 0;
				while (std::getline(ifs, line))
                {
					std::stringstream ss(line);
					
					int col=0;
                    KeyType key;
                    ValueType value;
                    T temp;

					while (ss >> temp)
                    {
						if(col<keySize)
                        {
                            key[col]=static_cast<int>(temp);
                        }
                        else
                        {
                            value[col-keySize]=temp;
                        }
						col++;
					}
					
					this->emplace(key,value);
					
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
		
//        /**********************************************************************/
//		bool readBIN(const std::string& filename)
//        {/*! Reads the binary file filename and stores its data in this
//          */
//            //bool success(false);
//			std::cout<<"Reading: "<<filename;
//			double t0(clock());
//			typedef std::pair<std::pair<int,int>, Eigen::Matrix<scalar,1,cols-2> > BinEdgeType;
//			BinaryFileReader<BinEdgeType> rE(filename);
//			for (unsigned int k=0;k<rE.size();++k)
//            {
////				this->insert(std::make_pair(rE[k].first,rE[k].second));
//				assert(this->insert(std::make_pair(rE[k].first,rE[k].second)).second && "COULD NOT INSERT EDGE AFTER BINARY READ.");
//			}
//			std::cout<<" ("<<this->size()<<" edges) ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
//            return rE.success();
//		}
		
	public:	
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		
        /**********************************************************************/
		IDreader () :
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
		static bool isGood(const int& frameN, const bool& useTXT)
        {
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
                    assert(0 && "NOT IMEPLEMETNTED YET");
					//success=readBIN(getFilename(frameN,false));
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
