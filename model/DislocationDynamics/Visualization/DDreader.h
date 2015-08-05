/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDreader_H_
#define model_DDreader_H_

#include <iostream>
#include <iomanip>
#include <string>
#include <dirent.h>
#include <map>
#include <set>
#include <vector>
#include <assert.h>

namespace model {
		
	/*************************************************************/
	/* DDreader **************************************************/
	/*************************************************************/
	class DDreader : public std::map<size_t,std::set<std::string> > {

		const std::string defaultColor;		// the default color for the console
		const std::string greenBoldColor;   // a bold green color
		const std::string blueColor;		// a bold blue color
		
		/* readDir *****************************************************/
		void readDir (const std::string& dirName, const std::string& extension)
        {
			
			const std::string prefix(dirName+"_");
		//	const std::string extension(".txt");
			
			std::cout<<"Opening directory: "<<dirName<<std::endl;
			
			DIR* const dp(opendir(dirName.c_str()));
            if(dp == NULL)
            {
                std::cout<<"Undable to open directory "<<dirName<<std::endl;
                assert(0 && "UNABLE TO OPEN DIRECTORY");
            }
			
			
			struct dirent *dirp;
			
			while ((dirp = readdir(dp)) != NULL) {
				std::string filename(dirp->d_name);
				if (filename.find(prefix)!=std::string::npos && filename.find(extension)!=std::string::npos)
                {
					std::string stringID=filename.substr(filename.find(prefix)+prefix.length());
					const size_t extPos(stringID.rfind(extension));
					assert(extPos>0);
					stringID.replace(extPos,extension.length(),"");
					const unsigned int fileID = atoi(stringID.c_str());
					this->operator[](fileID).insert(filename);
				}
			}
			closedir(dp);
		}
		

	public:
		
		/* Constructor *****************************************/
		DDreader() : defaultColor("\033[0m"),
		/*        */ greenBoldColor("\033[1;32m"),
		/*        */ blueColor("\033[0;34m"){
			
			readDir("E","txt");
			readDir("V","txt");
            
            readDir("E","bin");
            readDir("V","bin");
			//readDir("C");
			//readDir("G");
		}
		
		/* list ************************************************/		
		void list() const
        {
			std::cout<<greenBoldColor<<"listing files found:"<<std::endl;
			for (std::map<size_t,std::set<std::string> >::const_iterator fileIter=this->begin();fileIter!=this->end();++fileIter){
				std::cout<<std::setw(5)<<greenBoldColor<<fileIter->first<<blueColor<<" ";
				for (std::set<std::string>::const_iterator nameIter=fileIter->second.begin();nameIter!=fileIter->second.end();++nameIter){
					std::cout<<std::setw(11)<<(*nameIter)<<" ";
				}
				std::cout<<"\n";
			}
			std::cout<<defaultColor<<std::endl;
		}

	};
	
	
} // namespace model
#endif







