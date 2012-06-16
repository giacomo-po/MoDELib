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
#include <assert.h>
#include <Eigen/Dense>


namespace model {

template <char c, int cols, typename scalar, bool useTXT=true>
//class VertexReader : public std::map<int, Eigen::Matrix<scalar,1,cols-1> > {
	class VertexReader : public std::map<int, Eigen::Matrix<scalar,1,cols-1>,
	/*                                */ std::less<int>, 
	/*                                */ Eigen::aligned_allocator<std::pair<const int, Eigen::Matrix<scalar,1,cols-1> > > > {
	
	/*!	\brief A class template to read Vertex data from file. 
	 *	\param[cols] is the number of columns in the file.
	 *	First column is the Vertex ID, other columns are optional.
	 */
	
	int currentFrame;
	bool success;
	
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	//std::map<int, Eigen::Matrix<double,1,cols-1> > vertexMap;
	
	/*****************************************/
	/* Constructor */
	VertexReader() : currentFrame(-1), success(false) {}
	
	/*****************************************/
	/* isGood */
	static bool isGood(const int& frameN){
		std::stringstream filename;
		std::ifstream ifs;
		switch (useTXT){
			case true:
				//filename << "V/V_" << frameN << ".txt";; // TXT
				filename << c << "/" << c << "_" << frameN << ".txt";
				ifs.open( filename.str().c_str() , std::ifstream::in );
				break;
			default:
				filename << c << "/" << c << "_" << frameN << ".bin";
				//filename << "V/V_" << frameN << ".bin";; // BIN
				ifs.open( filename.str().c_str() , std::ifstream::in );
		}
		return ifs.good();
	}
	
	
	/*****************************************/
	/* read */
	bool read(const int& frameN){
		assert(frameN>=0);
		
		//success=false;
		
		if (frameN!=currentFrame){
			this->clear();
			
			//const std::string& filename = "V/V_"
			
			std::stringstream filename;
			filename << c << "/" << c << "_" << frameN << ".txt";
			//filename << "V/V_" << frameN << ".txt";
			
			scalar temp;
			int key;//, row;//, col;
			std::string line;
			Eigen::Matrix<scalar,1,cols-1> vertexPT=Eigen::Matrix<scalar,1,cols-1>::Zero();
			
			std::ifstream ifs ( filename.str().c_str() , std::ifstream::in );
			
			if (ifs.is_open()) {
				std::cout<<"Reading: "<<filename.str()<<std::endl;
				int row(0);
				while (std::getline(ifs, line)) {
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
				
				currentFrame=frameN;
				success=true;
			}
			else {
				std::cout<<"Unable to  open: "<<filename.str()<<std::endl;
				success=false;
			}
		}
		
		return success;
		
	}
	
	
	
	
};

} // namespace model
#endif
