/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez <ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TextFileParser_H_
#define model_TextFileParser_H_

#include <iostream>
#include <string>
#include <fstream>
#include <regex>
#include <vector>
#include <Eigen/Dense>



namespace model
{
    
    template<typename T>
    struct StringToScalar
    {
        
        static T toScalar(const std::string& key)
        {
            std::cout<<"Unknown conversion from std::string to "<< typeid(T).name()<<". Exiting."<<std::endl;
            exit(EXIT_FAILURE);
        }
    };
    
    template<>
    struct StringToScalar<int>
    {
        
        static int toScalar(const std::string& str)
        {
            return std::atoi(str.c_str());
        }
    };
    
    template<>
    struct StringToScalar<long>
    {
        
        static long toScalar(const std::string& str)
        {
            return std::atol(str.c_str());
        }
    };
    
    template<>
    struct StringToScalar<long long>
    {
        
        static long long toScalar(const std::string& str)
        {
            return std::atoll(str.c_str());
        }
    };
    
    template<>
    struct StringToScalar<unsigned long>
    {
        
        static unsigned long toScalar(const std::string& str)
        {
            return std::stoul(str.c_str());
        }
    };
    
    template<>
    struct StringToScalar<unsigned long long>
    {
        
        static unsigned long long toScalar(const std::string& str)
        {
            return std::stoull(str.c_str());
        }
    };
    
    template<>
    struct StringToScalar<float>
    {
        
        static float toScalar(const std::string& str)
        {
            return std::stof(str.c_str());
        }
    };
    
    template<>
    struct StringToScalar<double>
    {
        
        static double toScalar(const std::string& str)
        {
            return std::stod(str.c_str());
        }
    };
    
    
    template<>
    struct StringToScalar<long double>
    {
        
        static long double toScalar(const std::string& str)
        {
            return std::stold(str.c_str());
        }
    };
    
    
    class TextFileParser : public std::ifstream
    {
        
        template <typename Scalar>
        using EigenMapType=Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>, 0, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> >;
        
        /**********************************************************************/
        std::pair<std::string,std::string> readKey( const std::string& key
        //                               ,const bool& removeWitespaces=true
        )
        {
            this->seekg (0, this->beg); // reset the position of the next character at beginning for each read

            std::string line;
            std::string read;
            std::string comment;
            bool success=false;
            
            while (std::getline(*this, line))
            {
                
                const std::size_t foundKey=line.find(key);
                const std::size_t foundEqual=line.find("=");
                const std::size_t foundSemiCol=line.find(";");
                const std::size_t foundPound=line.find("#");

                
                if(   foundKey!=std::string::npos
                   && foundEqual!=std::string::npos
                   && foundSemiCol!=std::string::npos
                   && foundKey<foundEqual
                   && foundEqual<foundSemiCol
                   && foundSemiCol<foundPound
                   )
                {
                    read=line.substr(foundEqual+1,foundSemiCol-foundEqual-1);
                    if(foundPound!=std::string::npos)
                    {
                        comment=line.substr(foundPound,line.size()-foundPound);
                    }
                    success=true;
                    break;
                }
                
            }
            
            if(!success)
            {
                std::cout<<"File "<<fileName<<" does not cointain line with format:\n"<< key <<"=...;\n EXITING"<<std::endl;
                exit(EXIT_FAILURE);
            }
            
            //            if(removeWitespaces)
            //            {
            //                std::regex_replace( read, "\s", "" );
            //            }
            //            std::cout<<key<<"="<<read<<std::endl;
            
            return std::make_pair(read,comment);
        }
        
    public:
        
        const std::string fileName;
        
        /**********************************************************************/
        TextFileParser(const std::string& _fileName) :
        /* init */ std::ifstream(_fileName)
        /* init */,fileName(_fileName)
        {
            if(!this->is_open())
            {
                std::cout<<"File "<<fileName<<" cannot be opened. Exiting."<<std::endl;
                exit(EXIT_FAILURE);
            }
        }
        
        /**********************************************************************/
        std::string readString(const std::string& key,const bool&verbose=false)
        {
            const std::pair<std::string,std::string> strPair(readKey(key));
            //const std::string& read(readKey(key));
            if(verbose) std::cout<<key<<"="<<strPair.first<<" "<<strPair.second<<std::endl;
            return strPair.first;
        }
        
        /**********************************************************************/
        template<typename Scalar>
        Scalar readScalar(const std::string& key,const bool&verbose=false)
        {
            if(verbose) std::cout<<key<<"="<<std::flush;
//            const std::string str(readKey(key));
            const std::pair<std::string,std::string> strPair(readKey(key));
            const Scalar read(StringToScalar<Scalar>::toScalar(strPair.first));
            if(verbose) std::cout<<read<<" "<<strPair.second<<std::endl;
            return read;
        }
        
        /**********************************************************************/
        template<typename Scalar>
        std::vector<Scalar> readArray(const std::string& key,const bool&verbose=false)
        {
            this->seekg (0, this->beg); // reset the position of the next character at beginning for each read

            std::string line;
            std::string lines;
            std::vector<Scalar> array;
            bool success=false;
            
            while (std::getline(*this, line))
            {
                
                const std::size_t foundKey=line.find(key);
                const std::size_t foundEqual=line.find("=");
                const std::size_t foundPound=line.find("#");
                
                if(   foundKey!=std::string::npos
                   && foundEqual!=std::string::npos
                   && foundKey<foundEqual
                   && foundPound==std::string::npos
                   )
                {
                    lines+=line;
                    std::size_t foundSemiCol=line.find(";");
                    
                    if(   foundSemiCol!=std::string::npos
                       && foundEqual<foundSemiCol)
                    {
                        success=true;
                    }
                    else
                    {
                        while(std::getline(*this, line))
                        {
                            lines+=(" "+line);
                            foundSemiCol=lines.find(";");
                            if(foundSemiCol!=std::string::npos && foundEqual<foundSemiCol)
                            {
                                success=true;
                                break;
                            }
                        }
                    }
                                        
                    if(success)
                    {
                        Scalar temp(0.0);
                        std::stringstream ss(lines.substr(foundEqual+1,foundSemiCol-foundEqual-1));
                        while (ss >> temp)
                        {
                            array.push_back(temp);
                        }
                        break;
                    }
                }
                
            }
            
            if(!success)
            {
                std::cout<<"File "<<fileName<<" does not cointain line with format:\n"<< key <<"=...;\nEXITING"<<std::endl;
                exit(EXIT_FAILURE);
            }
            
            if(verbose)
            {
                std::cout<<key<<"=";
                for(const auto& val : array)
                {
                    std::cout<<val<<" ";
                }
                std::cout<<std::endl;
                
            }
            
            return array;
        }
        
        /**********************************************************************/
        template<typename Scalar>
         Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> readMatrix(const std::string& key,const size_t& rows,const size_t& cols,const bool&verbose=false)
        {

            const std::vector<Scalar> array=readArray<Scalar>(key,false);
            if(array.size()!=rows*cols)
            {
                std::cout<<"Error in reading matrix "<<key<<std::endl;
                std::cout<<"array.size="<<array.size()<<", is not equal to rows x cols ("<<rows<<"x"<<cols<<"). EXITING"<<std::endl;
                exit(EXIT_FAILURE);
            }
            
            EigenMapType<Scalar> em(array.data(), rows, cols, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(1, cols));
            if(verbose) std::cout<<key<<"=\n"<<em<<std::endl;
            return  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(em);
        }
        
        
    };
    
}
#endif
