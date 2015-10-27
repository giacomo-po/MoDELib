/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SCALARDATAREADER_H_
#define model_SCALARDATAREADER_H_

#include <iostream>
#include <typeinfo>
//#include <sstream>
//#include <fstream>
#include <model/MPI/MPIcout.h> // defines mode::cout
#include <model/Utilities/DataReader.h>

namespace model  {
	
	
	class ScalarDataReader : virtual public DataReader
    {
	
	private:
				
		
		
	public:
		
		template <typename T>
		bool readScalarInFile(const std::string & fileName_in, const std::string & varName_in, const size_t & occurrence_in, T & value);
		
		template <typename T>
		bool readScalarInFile(const std::string & fileName_in, const std::string & varName_in, T & value);
	//	const std::vector<std::vector<double> > & get_table(){return table;}
	
	};
	//////////////////////////////////////////////////////////////////////////	
	//////////////////////////////////////////////////////////////////////////	

	
	template <typename T>
	bool ScalarDataReader::readScalarInFile(const std::string & fileName_in, const std::string & varName_in, T & value)
    {
	
		return readScalarInFile(fileName_in, varName_in, 0, value);
	}
	
	template <typename T>
	bool ScalarDataReader::readScalarInFile(const std::string & fileName_in, const std::string & varName_in, const size_t & occurrence_in, T & value)
    {
		
		bool success=0;
		DataReader::readInFile(fileName_in, varName_in, occurrence_in);
		
		if(get_table().size()==1 && get_table()[0].size()==1)
        {
			value=get_table()[0][0];
//			model::cout<<varName_in<<"="<<value<<" (scalar "<< typeid(value).name()<<")"<<std::endl;
			model::cout<<varName_in<<"="<<value<<std::endl;
			success=1;
		}
		else{
			model::cout<<"Error in reading the "<<occurrence_in<<"-th occurrence of "<<varName_in << " in "<<fileName_in;
			model::cout<<": not a scalar."<<std::endl;
			exit (1);
		
		}
		
		return success;
	}
	
		
		//////////////////////////////////////////////////////////////////////////
}
#endif
