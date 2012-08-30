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
#include <map>
#include <assert.h>
#include <utility>
#include <Eigen/Dense>

#include <model/Utilities/BinaryFileReader.h>


namespace model {
	
	/*******************************************************************************************/
	template <char c, int cols, typename scalar=double>
	//	class EdgeReader : public std::map< std::pair<int, int>, Eigen::Matrix<scalar,1,cols-2> > {
	
	class EdgeReader : public std::map<std::pair<int, int>, Eigen::Matrix<scalar,1,cols-2>, 
	/*                              */ std::less<std::pair<int, int> >,   
	/*                              */ Eigen::aligned_allocator<std::pair<const std::pair<int, int>, Eigen::Matrix<scalar,1,cols-2> > > > {	
		
		/*!	\brief A class template to read Edge data from file. 
		 *	\param[cols] is the number of columns in the file.
		 *	First column is the source Vertex ID, second column is the sink Vertex ID. 
		 *  Other columns are optional.
		 */
		
		int currentFrame;
		bool success;
		
		/* readTXT ************************************/
		void readTXT(const std::string& filename){
			scalar temp;
			//int key1, key2, row, col;
			std::string line;
			
			Eigen::Matrix<scalar,1,cols-2> edgeNB = Eigen::Matrix<scalar,1,cols-2>::Zero();
			
			std::ifstream ifs ( filename.c_str() , std::ifstream::in );
			
			if (ifs.is_open()) {
				std::cout<<"Reading: "<<filename;
				double t0(clock());
				int row = 0;
				while (std::getline(ifs, line)) {
					std::stringstream ss(line);
					
					int col=0;
					int key1, key2;
					while (ss >> temp) {
						
						switch (col) {
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
				
				//				currentFrame=frameN;
				success=true;
				
				std::cout<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<"]"<<std::endl;
			}
			else {
				std::cout<<"Unable to  open:"<<filename<<std::endl;
				success=false;
			}
			
		}
		
		/* readBin *******************************/
		void readBIN(const std::string& filename){
			std::cout<<"Reading: "<<filename;
			double t0(clock());
			typedef std::pair<std::pair<int,int>, Eigen::Matrix<scalar,1,cols-2> > BinEdgeType;
			BinaryFileReader<BinEdgeType> rE(filename);
			std::cout<<" ["<<(clock()-t0)/CLOCKS_PER_SEC<<"]";
			double t1(clock());
			for (unsigned int k=0;k<rE.size();++k){
				this->insert(std::make_pair(rE[k].first,rE[k].second));
//				assert(this->insert(std::make_pair(rE[k].first,rE[k].second)).second && "COULD NOT INSERT EDGE AFTER BINARY READ.");
			}
			std::cout<<" ["<<(clock()-t1)/CLOCKS_PER_SEC<<"]"<<std::endl;
		}
		
	public:	
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		
		/*****************************************/
		/* Constructor*/
		EdgeReader () : currentFrame(-1), success(false) {}
		//		EdgeReader () :  success(false) {}
		
		/*****************************************/
		template <bool useTXT>
		static std::string getFilename(const int& frameN){
			std::stringstream filename;
			if(useTXT){
				filename << c << "/" << c << "_" << frameN << ".txt";
			}
			else{
				filename << c << "/" << c << "_" << frameN << ".bin";
			}
			return filename.str();
		}
		
		/* isGood *********************************/
		template <bool useTXT>
		static bool isGood(const int& frameN){
			/*!	Checks whether the file named "E/E_ \param[frameN] .txt" is good for reading.
			 */
			std::ifstream ifs ( getFilename<useTXT>(frameN).c_str() , std::ifstream::in );
			return ifs.good();
		}
		
		/*****************************************/
		/* read */
		template <bool useTXT>
		bool read(const int& frameN){
			assert(frameN>=0);
			if (frameN!=currentFrame){
				currentFrame=frameN;
				
				// clear the content of the map
				this->clear();
				
				if (useTXT){
					readTXT(getFilename<true >(frameN));
				}
				else{
					readBIN(getFilename<false>(frameN));
				}
			} 
			return success;
			
		}
		
	}; /* EdgeReader ******************************************************************************/
	
	/******************************************************************/
	/******************************************************************/
} /* namespace model */
#endif