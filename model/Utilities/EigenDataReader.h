/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EIGENDATAREADER_H_
#define model_EIGENDATAREADER_H_

#include <Eigen/Dense>
#include <model/MPI/MPIcout.h> // defines mode::cout
#include <model/Utilities/ScalarDataReader.h>


namespace model {
	
	
	class EigenDataReader : public ScalarDataReader {
		
	private:
		
		
		
	public:
		
		template <typename Derived>
		bool readVectorInFile(const std::string & fileName_in, const std::string & varName_in, 
							  const size_t & occurrence_in, Eigen::MatrixBase<Derived> & value);
		
		template <typename Derived>
		bool readVectorInFile(const std::string & fileName_in, const std::string & varName_in,
							  Eigen::MatrixBase<Derived> & value);
		
		template <typename Derived>
		bool readMatrixInFile(const std::string & fileName_in, const std::string & varName_in, 
							  const size_t & occurrence_in, Eigen::MatrixBase<Derived> & value);
		
		template <typename Derived>
		bool readMatrixInFile(const std::string & fileName_in, const std::string & varName_in,
							  Eigen::MatrixBase<Derived> & value);
		
		
	};
	//////////////////////////////////////////////////////////////////////////	
	//////////////////////////////////////////////////////////////////////////	
	
	
	template <typename Derived>
	bool EigenDataReader::readVectorInFile(const std::string & fileName_in, const std::string & varName_in,  
										   Eigen::MatrixBase<Derived> & value){
		
		return readVectorInFile(fileName_in, varName_in, 0, value);
	}
	
	template <typename Derived>
	bool EigenDataReader::readVectorInFile(const std::string & fileName_in, const std::string & varName_in, 
										   const size_t & occurrence_in,  Eigen::MatrixBase<Derived> & value){
		
		bool success=0;
		ScalarDataReader::readInFile(fileName_in, varName_in, occurrence_in);
		
		if(get_table().size()==1)
        {
			size_t Ncols=get_table()[0].size();
			value.derived().resize(Ncols);
			for (size_t k=0;k<Ncols;++k)
            {
				value.derived()(k)=get_table()[0][k];
            }
			//&& >1
			
			//value.derived()=Derived::Map(&get_table()[0][0], get_table()[0].size());
			if(value.derived().cols()==1)
            {
			model::cout<<varName_in<<"="<<std::endl<<value.derived().transpose()<<std::endl;
            }
            else
            {
                model::cout<<varName_in<<"="<<std::endl<<value.derived()<<std::endl;
            }
			success=1;
        }
		else
        {
			model::cout<<"Error in reading the "<<occurrence_in<<"-th occurrence of "<<varName_in << " in "<<fileName_in;
			model::cout<<": not a vector."<<std::endl;
			//exit (1);
		}
		
		return success;
	}
	
	//////////////////////////////////////////////////////////////////////////
	template <typename Derived>
	bool EigenDataReader::readMatrixInFile(const std::string & fileName_in, const std::string & varName_in,  
										   Eigen::MatrixBase<Derived> & value){
		
		return readMatrixInFile(fileName_in, varName_in, 0, value);
	}
	
	template <typename Derived>
	bool EigenDataReader::readMatrixInFile(const std::string & fileName_in, const std::string & varName_in, 
										   const size_t & occurrence_in,  Eigen::MatrixBase<Derived> & value){
		
		bool success=0;
		
		ScalarDataReader::readInFile(fileName_in, varName_in, occurrence_in);
		
		size_t Nrows=get_table().size();
		size_t Ncols=get_table()[0].size();
		for (size_t k=0; k<Nrows;++k)
        {
			if(get_table()[k].size()!=Ncols)
            {
				model::cout<<"Number of rows mismatch"<<std::endl;
				exit(1);
            }
		};
		
		
		
		if(Nrows>=1 )
        {
			
			value.derived().resize(Nrows,Ncols);
			
			
			for (size_t r=0;r<Nrows;++r)
            {
				for (size_t c=0;c<Ncols;++c)
                {
					value.derived()(r,c)=get_table()[r][c];
					
				}
            }
			
			
			model::cout<<varName_in<<"="<<std::endl<<value<<std::endl;
			success=1;
		}
		else
        {
			model::cout<<"Error in reading the "<<occurrence_in<<"-th occurrence of "<<varName_in << " in "<<fileName_in;
			model::cout<<": not a matrix."<<std::endl;
			exit (1);
			
		}
		
		
		
		
		return success;
	}
	
	//////////////////////////////////////////////////////////////////////////
}
#endif
