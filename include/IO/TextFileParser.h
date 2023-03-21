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

//#include <filesystem>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <regex>
#include <vector>
#include <set>
#include <map>
#include <regex>
#include <Eigen/Dense>



#include <TerminalColors.h>

namespace model
{

template<typename T>
struct StringToScalar
{
    
    static T toScalar(const std::string& key)
    {
        throw std::runtime_error("Unknown conversion from std::string "+key+" to "+typeid(T).name()+".");
        //            std::cout<<"Unknown conversion from std::string "<<key<<" to "<< typeid(T).name()<<". Exiting."<<std::endl;
        //            exit(EXIT_FAILURE);
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
//    std::pair<std::string,std::string> readKey( const std::string& key
//    //                               ,const bool& removeWitespaces=true
//    )
//    {
//        this->seekg (0, this->beg); // reset the position of the next character at beginning for each read
//    this->clear();
//    this->seekg(0);
//
//        std::string line;
//        std::string read;
//        std::string comment;
//        bool success(false);
//
//        while (std::getline(*this, line))
//        {
//            const size_t foundKey=line.find(key);
//            const size_t foundEqual=line.find("=");
//            const std::string keyRead(removeSpaces(line.substr(0,foundEqual)));
//
//            if(keyRead==key)
//            {
//                const size_t foundSemiCol=line.find(";");
//                const size_t foundPound=line.find("#");
//
//                if(   foundKey!=std::string::npos
//                   && foundEqual!=std::string::npos
//                   && foundSemiCol!=std::string::npos
//                   && foundKey<foundEqual
//                   && foundEqual<foundSemiCol
//                   && foundSemiCol<foundPound
//                   )
//                {
//                    read=line.substr(foundEqual+1,foundSemiCol-foundEqual-1);
//                    if(foundPound!=std::string::npos)
//                    {
//                        comment=line.substr(foundPound,line.size()-foundPound);
//                    }
//                    success=true;
//                    break;
//
//                }
//            }
//        }
//
//        if(!success)
//        {
//            throw std::runtime_error("File "+fileName+" does not cointain line with format "+key+"=...;");
//            //                std::cout<<"File "<<fileName<<" does not cointain line with format:\n"<< key <<"=...;\n EXITING"<<std::endl;
//            //                exit(EXIT_FAILURE);
//        }
//
//        //            if(removeWitespaces)
//        //            {
//        //                std::regex_replace( read, "\s", "" );
//        //            }
//        //            std::cout<<key<<"="<<read<<std::endl;
//
//        return std::make_pair(read,comment);
//    }
    
    std::vector<std::pair<std::string,std::string>> readKey( const std::string& key)
    {
        this->clear();
        this->seekg(0);
        std::vector<std::pair<std::string,std::string>> returnVector;
        std::string line;

        while (std::getline(*this, line))
        {
            const size_t foundKey=line.find(key);
            const size_t foundEqual=line.find("=");
            
            if(   foundKey!=std::string::npos
               && foundEqual!=std::string::npos
               && foundKey<foundEqual)
            {
                const std::string keyRead(removeSpaces(line.substr(0,foundEqual)));

                if(keyRead==key)
                {
                    
                    
                    const size_t foundSemiCol=line.find(";");
                    const size_t foundPound=line.find("#");
                                        
                    if(   //foundKey!=std::string::npos
                        foundEqual!=std::string::npos
                       && foundSemiCol!=std::string::npos
                       //&& foundKey<foundEqual
                       && foundEqual<foundSemiCol
                       && foundSemiCol<foundPound
                       )
                    {
                        const std::string read(line.substr(foundEqual+1,foundSemiCol-foundEqual-1));

                        if(foundPound!=std::string::npos)
                        {
                            const std::string comment(line.substr(foundPound,line.size()-foundPound));
                            returnVector.emplace_back(read,comment);
                        }
                        else
                        {
                            returnVector.emplace_back(read,std::string());
                        }
                    }
                }
            }
        }
        if(returnVector.size()==0)
        {
            throw std::runtime_error("File "+fileName+" does not cointain line with format "+key+"=...;");
        }

        return returnVector;
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
            throw std::runtime_error("File "+fileName+" cannot be opened.");
            //                std::cout<<"File "<<fileName<<" cannot be opened. Exiting."<<std::endl;
            //                exit(EXIT_FAILURE);
        }
    }
    
//    TextFileParser(const std::filesystem::path& _fileName) :
//    /* init */ std::ifstream(_fileName.string())
//    /* init */,fileName(_fileName.string())
//    {
//        if(!this->is_open())
//        {
//            throw std::runtime_error("File "+fileName+" cannot be opened.");
//            //                std::cout<<"File "<<fileName<<" cannot be opened. Exiting."<<std::endl;
//            //                exit(EXIT_FAILURE);
//        }
//    }
    
    static std::string removeSpaces(std::string key)
    {
        key.erase(std::remove_if(key.begin(), key.end(), [](unsigned char x) { return std::isspace(x); }), key.end());
        return key;
    }
    
    /**********************************************************************/
    std::string readString(const std::string& key,const bool&verbose=false)
    {
        const std::pair<std::string,std::string> strPair(readKey(key)[0]);
        if(verbose) std::cout<<cyanColor<<key<<"="<<strPair.first<<" "<<strPair.second<<defaultColor<<std::endl;
        return strPair.first;
    }
    
    /**********************************************************************/
    std::vector<std::pair<std::string,std::string>> readStringVector(const std::string& key)
    {
        return readKey(key);
    }
    
    /**********************************************************************/
    template<typename Scalar>
    Scalar readScalar(const std::string& key,const bool&verbose=false)
    {
        if(verbose) std::cout<<cyanColor<<key<<"="<<std::flush;
        const std::pair<std::string,std::string> strPair(readKey(key)[0]);
        const Scalar read(StringToScalar<Scalar>::toScalar(strPair.first));
        if(verbose) std::cout<<read<<" "<<strPair.second<<defaultColor<<std::endl;
        return read;
    }
    
    /**********************************************************************/
    template<typename Scalar>
    std::set<Scalar> readSet(const std::string& key,const bool&verbose=false)
    {
        std::vector<Scalar> tempV(readArray<Scalar>(key,false));
        std::set<Scalar> tempS;
        for(const auto& val : tempV)
        {
            tempS.insert(val);
        }
        if(verbose)
        {
            std::cout<<cyanColor<<key<<"=";
            for(const auto& val : tempS)
            {
                std::cout<<" "<<val;
            }
            std::cout<<"; "<<defaultColor<<std::endl;
        }
        return tempS;
    }
    
    template<typename KeyType,typename Scalar>
    std::tuple<bool,bool,KeyType,std::vector<Scalar>,std::string> splitMapLine(const std::string& line) const
    {
        std::tuple<bool,bool,KeyType,std::vector<Scalar>,std::string> temp(std::make_tuple(false,false,KeyType(),std::vector<Scalar>(),std::string()));
        
        const size_t foundCol=line.find(":");
        const size_t foundSemiCol=line.find(";");
        const size_t foundPound=line.find("#");

        if(foundCol!=std::string::npos
           && (foundPound==std::string::npos || foundPound>foundSemiCol))
        {
            std::get<0>(temp)=true; // valid line
            std::get<1>(temp)=(foundPound!=std::string::npos); // has semicolumn
            std::get<2>(temp)=StringToScalar<KeyType>::toScalar(removeSpaces(line.substr(0,foundCol)));

            const std::string valStr(foundSemiCol==std::string::npos? line.substr(foundCol+1,line.size()-foundCol-1) : line.substr(foundCol+1,foundSemiCol-foundCol-1));
            Scalar tempVal;
            std::stringstream ss(valStr);
            while (ss >> tempVal)
            {
                std::get<3>(temp).push_back(tempVal);
            }
            if(foundPound!=std::string::npos)
            {
                std::get<4>(temp)=line.substr(foundPound+1,line.size()-foundPound-1);
            }
        }
        return temp;
    }
    
    /**********************************************************************/
    template<typename KeyType,typename Scalar>
    std::map<KeyType,std::vector<Scalar>> readArrayMap(const std::string& key,const bool&verbose=false)
    {
//        this->seekg (0, this->beg); // reset the position of the next character at beginning for each read
        this->clear();
        this->seekg(0);

        std::string line;
//        std::string lines;
//        std::string comment;
//        std::vector<Scalar> array;
        bool success=false;
        std::map<KeyType,std::vector<Scalar>> temp;

        while (std::getline(*this, line) || !success)
        {

            const size_t foundKey=line.find(key);
            const size_t foundEqual=line.find("=");

            const std::string keyRead(removeSpaces(line.substr(0,foundEqual)));

            if(keyRead==key
               && foundKey!=std::string::npos
               && foundEqual!=std::string::npos
               && foundKey<foundEqual)
            {
                
                success=true;
//                size_t foundSemiCol=line.find(";");
//                size_t foundPound=line.find("#");

//                if(
////                   && (foundPound==std::string::npos || foundPound>foundSemiCol)
//                   )
//                {
                    
                    auto tup(splitMapLine<KeyType,Scalar>(line.substr(foundEqual+1,line.size()-foundEqual-1)));
                if(std::get<0>(tup))
                {// valid line
                    temp.emplace(std::get<2>(tup),std::get<3>(tup));
                }
                
                if(std::get<1>(tup))
                {// semi column found
                    break;
                }
                else
                {// read following lines
                    while(std::getline(*this, line))
                    {
                        tup=splitMapLine<KeyType,Scalar>(line);
                        
                        if(std::get<0>(tup))
                        {// valid line
                            temp.emplace(std::get<2>(tup),std::get<3>(tup));
                        }
                        
                        if(std::get<1>(tup))
                        {// semi column found
                            break;
                        }
                    }
                }
                
            }
        }

        if(!success)
        {
            throw std::runtime_error("File "+fileName+" does not cointain line with format "+key+"=...;");
        }

        if(verbose)
        {
            std::cout<<cyanColor<<key<<"="<<std::endl;
            for(const auto& pair : temp)
            {
                std::cout<<pair.first<<" : ";
                for(const auto& val : pair.second)
                {
                    std::cout<<val<<" ";
                }
                std::cout<<std::endl;
            }
//            std::cout<<"; "<<comment<<defaultColor<<std::endl;

        }

        return temp;
    }
    
    /**********************************************************************/
    template<typename Scalar>
    std::vector<Scalar> readArray(const std::string& key,const bool&verbose=false)
    {
//        this->seekg (0, this->beg); // reset the position of the next character at beginning for each read
        this->clear();
        this->seekg(0);
        
        std::string line;
        std::string lines;
        std::string comment;
        std::vector<Scalar> array;
        bool success=false;
        
        while (std::getline(*this, line))
        {
            
            const size_t foundKey=line.find(key);
            const size_t foundEqual=line.find("=");
            
            const std::string keyRead(removeSpaces(line.substr(0,foundEqual)));

            if(keyRead==key)
            {
                size_t foundSemiCol=line.find(";");
                size_t foundPound=line.find("#");
                
                if(   foundKey!=std::string::npos
                   && foundEqual!=std::string::npos
                   && foundKey<foundEqual
                   && (foundPound==std::string::npos || foundPound>foundSemiCol)
                   )
                {
                    lines+=line;
                    foundSemiCol=line.find(";");
                    
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
                        foundPound=line.find("#");
                        if(foundPound!=std::string::npos)
                        {
                            comment=lines.substr(foundPound,lines.size()-foundPound);
                        }
                        
                        
                        Scalar temp;
                        std::stringstream ss(lines.substr(foundEqual+1,foundSemiCol-foundEqual-1));
                        while (ss >> temp)
                        {
                            array.push_back(temp);
                        }
                        break;
                    }
                }
            }
        }
        
        if(!success)
        {
            throw std::runtime_error("File "+fileName+" does not cointain line with format "+key+"=...;");
        }
        
        if(verbose)
        {
            std::cout<<cyanColor<<key<<"=";
            for(const auto& val : array)
            {
                std::cout<<" "<<val;
            }
            std::cout<<"; "<<comment<<defaultColor<<std::endl;
            
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
            throw std::runtime_error("Error in reading matrix "+key+": array.size="+std::to_string(array.size())+" is not equal to rows x cols ("+std::to_string(rows)+"x"+std::to_string(cols)+").");
            //                std::cout<<"Error in reading matrix "<<key<<std::endl;
            //                std::cout<<"array.size="<<array.size()<<", is not equal to rows x cols ("<<rows<<"x"<<cols<<"). EXITING"<<std::endl;
            //                exit(EXIT_FAILURE);
        }
        
        EigenMapType<Scalar> em(array.data(), rows, cols, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(1, cols));
        if(verbose) std::cout<<cyanColor<<key<<"=\n"<<em<<defaultColor<<std::endl;
        return  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(em);
    }
    
    /**********************************************************************/
    template<typename Scalar,int rows,int cols>
    Eigen::Matrix<Scalar,rows,cols> readMatrix(const std::string& key,const bool&verbose=false)
    {
        return  readMatrix<Scalar>(key,rows,cols,verbose);
    }
    
    /**********************************************************************/
    template<typename Scalar>
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> readMatrixCols(const std::string& key,const size_t& cols,const bool&verbose=false)
    {
        
        const std::vector<Scalar> array=readArray<Scalar>(key,false);
        if(array.size()%cols!=0)
        {
            throw std::runtime_error("Error in reading matrix "+key+": array.size="+std::to_string(array.size())+" is not a multiple of cols ("+std::to_string(cols)+").");
            //                std::cout<<"Error in reading matrix "<<key<<std::endl;
            //                std::cout<<"array.size="<<array.size()<<", is not a multiple of cols ("<<cols<<"). EXITING"<<std::endl;
            //                exit(EXIT_FAILURE);
        }
        const size_t rows(array.size()/cols);
        EigenMapType<Scalar> em(array.data(), rows, cols, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(1, cols));
        if(verbose) std::cout<<cyanColor<<key<<"=\n"<<em<<defaultColor<<std::endl;
        return  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(em);
    }
    
    /**********************************************************************/
    template<typename Scalar>
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> readMatrixRows(const std::string& key,const size_t& rows,const bool&verbose=false)
    {
        
        const std::vector<Scalar> array=readArray<Scalar>(key,false);
        if(array.size()%rows!=0)
        {
            throw std::runtime_error("Error in reading matrix "+key+": array.size="+std::to_string(array.size())+" is not a multiple of rows ("+std::to_string(rows)+").");
            
            //                std::cout<<"Error in reading matrix "<<key<<std::endl;
            //                std::cout<<"array.size="<<array.size()<<", is not a multiple of rows ("<<rows<<"). EXITING"<<std::endl;
            //                exit(EXIT_FAILURE);
        }
        const size_t cols(array.size()/rows);
        EigenMapType<Scalar> em(array.data(), rows, cols, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(1, cols));
        if(verbose) std::cout<<cyanColor<<key<<"=\n"<<em<<defaultColor<<std::endl;
        return  Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>(em);
    }
    
};

}
#endif
