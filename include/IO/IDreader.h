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
#include <utility>
#include <chrono>
#include <array>
//#include <Eigen/Dense>
#include <BinaryFileReader.h>

namespace model
{
    template <int keySize>
    struct IDreaderKeySelector
    {
        static_assert(keySize>0,"keySize must be >0");
        typedef long long int SingleKeyType;
        typedef std::array<SingleKeyType, keySize> KeyType;
        
        static SingleKeyType& keyElement(KeyType& key,const size_t& k)
        {
            return key[k];
        }
    };
    
    template <>
    struct IDreaderKeySelector<1>
    {
        typedef long long int SingleKeyType;
        typedef SingleKeyType KeyType;
        
        static SingleKeyType& keyElement(KeyType& key,const size_t&)
        {
            return key;
        }
    };
    
    
    /*!	\brief A class template
     */
	template <int keySize, int valueSize, typename T>
    class IDreader : public std::map<typename IDreaderKeySelector<keySize>::KeyType, std::array<T, valueSize> >
    {
        static_assert(valueSize>=0,"valueSize must be >0");

        typedef typename IDreaderKeySelector<keySize>::SingleKeyType SingleKeyType;
        typedef typename IDreaderKeySelector<keySize>::KeyType KeyType;
        typedef std::array<T, valueSize> ValueType;
        typedef std::pair<KeyType, ValueType > BinType;

		long long int currentFrame;
        std::map<std::string,size_t> labelsMap;
		
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


				size_t row = 0;
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
                            if(temp!=std::round(temp))
                            {
                                throw std::runtime_error("IDreader: key "+std::to_string(temp)+" must be an integer.");
                            }
                            IDreaderKeySelector<keySize>::keyElement(key,col)=static_cast<SingleKeyType>(temp);
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
				
                std::cout<<" ("<<this->size()<<" entries) ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
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
            std::cout<<"Reading: "<<filename<<std::flush<<" (";
            const auto t0=std::chrono::system_clock::now();
			BinaryFileReader<BinType> rE(filename);
            for (const auto& pair : rE)
            {
                //				this->insert(std::make_pair(rE[k].first,rE[k].second));
                const bool success=this->emplace(pair.first,pair.second).second;
                if(!success)
                {
                    throw std::runtime_error("IDreader: cannot insert key="+std::to_string(pair.first)+".");
                }
            }
            std::cout<<this->size()<<"elements in "<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec)"<<std::endl;
            return rE.success();
		}
		
	public:	
		const std::string fileNamePrefix;
		
        /**********************************************************************/
//		IDreader () :
//        /* init list */ currentFrame(-1), folderSuffix(std::string(""))
// //       /* init list */ success(false)
//        {/*! Constructor initializes currentFrame to -1 so that the statement
//          *  frameN!=currentFrame will initially return false
//          */
//        }

        /***Overloaded constructor***********************************************/
		IDreader (const std::string& fileNamePrefix_in) :
        /* init list */ currentFrame(-1),
        /* init list */ fileNamePrefix(fileNamePrefix_in)
 //       /* init list */ success(false)
        {/*! Constructor initializes currentFrame to -1 so that the statement
          *  frameN!=currentFrame will initially return false
          */
        }
        
        /**********************************************************************/
        void readLabelsFile(const std::string& labelsFileName)
        {
            std::ifstream labelsFile(labelsFileName.c_str(), std::ifstream::in);
            if (labelsFile.is_open())
            {
                std::cout<<"reading labels file: "<<labelsFileName;
                double t0(clock());
                std::string line;
                
                std::getline(labelsFile, line); // labels of ID
                
                size_t row = 0;
                while (std::getline(labelsFile, line))
                {
                                        
                    const bool success=labelsMap.emplace(line,row).second;
                    if(!success)
                    {
                        std::cout<<"Unable to insert label "<<line<<std::endl;
                    }
                    row++;
                }
                
                
                labelsFile.close();
                
                
                std::cout<<" ("<<labelsMap.size()<<" labels) ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
            }
            else
            {
                std::cout<<"Unable to read labels file:"<<labelsFileName<<std::endl;
            }
        }
		
        /**********************************************************************/
        const T& operator()(const KeyType& key,const std::string& label)
        {
            const auto rowIter(this->find(key));
            assert(rowIter!=this->end() && "ID not found.");
            const auto colIter(labelsMap.find(label));
            assert(colIter!=labelsMap.end() && "label not found.");
            return rowIter->second.operator[](colIter->second);
        }
        
        /**********************************************************************/
        const T& first(const std::string& label)
        {
            const auto rowIter(this->begin());
            const auto colIter(labelsMap.find(label));
            assert(colIter!=labelsMap.end() && "label not found.");
            return rowIter->second.operator[](colIter->second);
        }
        
        /**********************************************************************/
        const T& last(const std::string& label)
        {
            const auto rowIter(this->rbegin());
            const auto colIter(labelsMap.find(label));
            assert(colIter!=labelsMap.end() && "label not found.");
            return rowIter->second.operator[](colIter->second);
        }
        
        /**********************************************************************/
		std::string getFilename(const long long int& frameN, const bool& useTXT) const
        {
//			std::stringstream filename;
			if(useTXT)
            {
                return fileNamePrefix+"_"+std::to_string(frameN)+".txt";
//				filename << c <<folderSuffix<< "/" << c << "_" << frameN << ".txt";
			}
			else
            {
                return fileNamePrefix+"_"+std::to_string(frameN)+".bin";
//				filename << c <<folderSuffix<< "/" << c << "_" << frameN << ".bin";
			}
//			return filename.str();
		}
		
        /**********************************************************************/
		bool isGood(const long long int& frameN, const bool& useTXT) const
        {
			/*!	Checks whether the file named "E/E_ \param[frameN] .txt" is good for reading.
			 */
			std::ifstream ifs ( getFilename(frameN,useTXT).c_str() , std::ifstream::in );
			return ifs.good();
		}
		
		/**********************************************************************/
		bool read(const long long int& frameN, const bool& useTXT)
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
