/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DATAREADER_H_
#define model_DATAREADER_H_
#include <assert.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
//#include <algorithm>
#include <model/MPI/MPIcout.h> // defines mode::cout


namespace model
{
	
	
	class DataReader
    {
		
	private:
		std::ifstream inputFile;
		std::string fileName;
		std::string varName;
		
		std::string line;
		std::string clean_line;
		
		std::istringstream iss;
		
		double temp;
		std::vector<double> row;
		std::vector<std::vector<double> > table;
		
		//size_t occurrence;
		//size_t count;
		
		
		bool menage_occurences(const size_t & occurrence);
		
		void readLines();
        bool check_name_equal(const std::string&);
		
		void line2row();
		
		
		
	public:
		
		bool readInFile(const std::string & fileName_in, const std::string & varName_in, const size_t & occurrence);
		bool readInFile(const std::string & fileName_in, const std::string & varName_in);
		const std::vector<std::vector<double> > & get_table(){return table;}
		
	};
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	
	
	//////////////////////////////////////////////////////////////////////////
	bool DataReader::readInFile(const std::string & fileName_in, const std::string & varName_in, const size_t & occurrence)
    {
		
		table.clear();		// VERY IMPORTANT
		
		fileName=fileName_in;
		varName=varName_in;
		//occurrence=occurrence;
		
		return menage_occurences(occurrence);
	}
	
	//////////////////////////////////////////////////////////////////////////
	bool DataReader::readInFile(const std::string & fileName_in, const std::string & varName_in)
    {
		return readInFile(fileName_in, varName_in, 0);
	}
	
	//////////////////////////////////////////////////////////////////////////
	bool DataReader::menage_occurences(const size_t & occurrence){
		
		bool success=0;
		int  count=-1;
		
		inputFile.open(fileName.c_str(), std::ios::in);
		
		
		if(inputFile.good()){
			
			
			while(count<int(occurrence))
            {
				//! 1- get the next line
				std::getline(inputFile, line);
				
				//! 2- check if the line contains the variable name and eventually increment counter
				count+=check_name_equal(varName);
				//	model::cout<<"count="<<count<<std::endl;
				
				
				//! 3- if "count" reaches "occurrence" read the variable and exit
				if (count==int(occurrence))
                {
//					model::cout<<"Reading "<<occurrence<<"-th occurrence of "<<varName<< " in "<<fileName<<std::endl;
					readLines();
					success=1;
					break;
                }
				
				//! 4- if end of file is reached exit
				if (inputFile.eof())
                {
					model::cout<<"Warning. In reading the "<<occurrence<<"-th occurrence of "<<varName << " in "<<fileName;
					model::cout<<": found only "<<count<<" occurrences."<<std::endl;
					
					model::cout<<line<<std::endl;
					
					//		exit (1);
					break;
                }
				
			}
		}
		else{
			model::cout<<"Cannot read file "<<fileName<<std::endl;
			assert(success);
		}
		
		inputFile.close();
		//	model::cout<<"success="<<success<<std::endl;
		
		
		return success;
	}
	
	
	//////////////////////////////////////////////////////////////////////////
    bool DataReader::check_name_equal(const std::string& temp)
    {
//		return line.find(varName)==0 && line.find("=")!=std::string::npos;
//        std::replace(temp,' ','');
        return line.find(temp)==0 && (line[temp.size()/(sizeof 'a')]=='=' || (line[temp.size()/(sizeof 'a')]==' ' && line.find("=")!=std::string::npos));
	}
	
	//////////////////////////////////////////////////////////////////////////
	void DataReader::readLines()
    {
		
		table.clear();
		
		//! clean line from stuff before "="
		line=line.substr(line.find("=")+1,std::string::npos);
		
		//! Add the current line
		line2row();
		
		//! Keep adding lines until ";" is found
		while(line.find(";")==std::string::npos)
        {  // read lines until you find a ";"
			
			std::getline(inputFile, line);
			line2row();
			
			if (inputFile.eof())
            {
				//model::cout<<"Error in reading the "<<occurrence<<"-th occurrence of "<<varName << " in "<<fileName;
				model::cout<<"Error: file ended before ';'."<<std::endl;
				exit (1);
				break;
            }
		}
		
		
//		for (unsigned int r=0; r<table.size(); ++r){
//			for (unsigned int c=0; c<table[r].size(); ++c){
//				model::cout<<table[r][c]<<" ";}
//			model::cout<<std::endl;
//        }
		
	}
	
	//////////////////////////////////////////////////////////////////////////
	void DataReader::line2row(){
		
		//! clean_line is line without ";"
		clean_line=line.substr(0,line.find(";"));
		
		row.clear();
		
		iss.clear();
		iss.str (clean_line);
		
		int n=0;
		while (!iss.eof()){
			iss >> temp;
			row.push_back(temp);
			n++;
		}
		
		table.push_back(row);
	}
	
	//////////////////////////////////////////////////////////////////////////
}
#endif
