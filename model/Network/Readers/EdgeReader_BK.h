/* This file is part of mmdl, the Mechanics of Material Defects Library.
 *
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 *
 * mmdl is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef mmdl_EDGEREADER_H_
#define mmdl_EDGEREADER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <assert.h>
#include <utility>
#include <Eigen/Dense>

namespace mmdl {
	
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
		
	public:	
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
				
		/*****************************************/
		/* Constructor*/
		EdgeReader () : currentFrame(-1), success(false) {}
		
		
		
		static bool isGood(const int& frameN){
			/*!	Checks whether the file named "E/E_ \param[frameN] .txt" is good for reading.
			 */
			std::stringstream filename;
			filename << c << "/" << c << "_" << frameN << ".txt";
			std::ifstream ifs ( filename.str().c_str() , std::ifstream::in );
			return ifs.good();
		}
		
		/*****************************************/
		/* read */
		bool read(const int& frameN){
			assert(frameN>=0);
			//=false;
			
			
			if (frameN!=currentFrame){
				
				// clear the content of the map
				this->clear();
				
				std::stringstream filename;
				filename << c << "/" << c << "_" << frameN << ".txt";
				
				scalar temp;
				//int key1, key2, row, col;
				std::string line;
				
				Eigen::Matrix<scalar,1,cols-2> edgeNB = Eigen::Matrix<scalar,1,cols-2>::Zero();
				
				std::ifstream ifs ( filename.str().c_str() , std::ifstream::in );
				
				if (ifs.is_open()) {
					std::cout<<"Reading: "<<filename.str()<<std::endl;
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
					
					currentFrame=frameN;
					success=true;
					
					
				}
				else {
					std::cout<<"Unable to  open:"<<filename.str()<<std::endl;
					success=false;
				}
			} 
			return success;
			
		}
		
		
		
	
	}; /* EdgeReader ******************************************************************************/
	
} /* namespace mmdl */
#endif
