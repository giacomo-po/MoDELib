/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDauxIO_H_
#define model_DDauxIO_H_

#include <vector>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <cfloat>      // std::ifstream
#include <DDbaseIO.h>
#include <GlidePlaneBoundaryIO.h>
#include <GlidePlaneFactory.h>



namespace model
{
    
    
    template <int dim>
    struct DDauxIO : public DDbaseIO
    /*            */,private std::vector<GlidePlaneBoundaryIO<dim>>
    
    {
        
        DDauxIO(const std::string& suffix="") :
        /* init */ DDbaseIO("evl","ddAux",suffix)
        {
            
        }
        


        
        /**********************************************************************/
        const std::vector<GlidePlaneBoundaryIO<dim>>& glidePlanesBoundaries() const
        {
            return *this;
        }
        
        std::vector<GlidePlaneBoundaryIO<dim>>& glidePlanesBoundaries()
        {
            return *this;
        }
        
        /**********************************************************************/
        void writeTxt(const GlidePlaneFactory<dim>& dn,
                             const long int& runID)
        {
            
            const std::string filename(this->getTxtFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
            if(file.is_open())
            {
                // Write header
                size_t nGP(0);
                for(const auto& glidePlane : dn.glidePlanes())
                {
                    nGP+=glidePlane.second->meshIntersections.size();
                }
                file<<nGP<<"\n";
                
                // Write Nodes
                for(const auto& glidePlane : dn.glidePlanes())
                {
                    for(const auto& seg : glidePlane.second->meshIntersections)
                    {
                        file<<GlidePlaneBoundaryIO<dim>(glidePlane.second->sID,seg)<<"\n";
                    }
                }

                file.close();
            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
        }
        
        /**********************************************************************/
        void writeTxt(const size_t& runID,
                             const std::vector<GlidePlaneBoundaryIO<dim>> gpBoundaries)
        {
            
            const std::string filename(this->getTxtFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
//            std::cout<<"Writing to "<<filename<<std::endl;
            if(file.is_open())
            {
                // Write header
                file<<gpBoundaries.size()<<"\n";
                
                // Write Nodes
                for(const auto& gpBnd : gpBoundaries)
                {
                    file<<gpBnd<<"\n";
                }
                
                file.close();
            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
        }
        
        /**********************************************************************/
        template<typename DislocationNodeType>
        void writeBin(const DislocationNodeType& dn,
                             const long int& runID)
        {

            const std::string filename(this->getBinFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
            if(file.is_open())
            {
                // Write header
                size_t nGP(0);
                for(const auto& glidePlane : dn.glidePlanes())
                {
                    nGP+=glidePlane.second->meshIntersections.size();
                }
                binWrite(file,nGP);
                
                // Write GlidePlaneBoundaries
                for(const auto& glidePlane : dn.glidePlanes())
                {
                    for(const auto& seg : glidePlane->second->meshIntersections)
                    {
                        binWrite(file,GlidePlaneBoundaryIO<dim>(seg));
                    }
                }
                
                file.close();
            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
        }
        
        /**********************************************************************/
        void writeBin(const size_t& runID,
                             const std::vector<GlidePlaneBoundaryIO<dim>> gpBoundaries)
        {
            
            const std::string filename(this->getBinFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
            if(file.is_open())
            {
                // Write header
                const size_t nGP(gpBoundaries.size());
                binWrite(file,nGP);

                
                // Write Nodes
                for(const auto& gpb : gpBoundaries)
                {
                    binWrite(file,gpb);
                }
                
                file.close();
            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
        }
        
        
        /**********************************************************************/
        void readBin(const size_t& runID)
        {
            const std::string filename(this->getBinFilename(runID));
            
            std::ifstream infile (filename.c_str(), std::ios::in|std::ios::binary);
            if(infile.is_open())
            {
                const auto t0=std::chrono::system_clock::now();
                model::cout<<"reading "<<filename<<std::flush;
                
                // Read header
                size_t sizeGP;
                infile.read (reinterpret_cast<char*>(&sizeGP), 1*sizeof(sizeGP));
                
                // Read vertices
                glidePlanesBoundaries().resize(sizeGP);
                infile.read (reinterpret_cast<char*>(glidePlanesBoundaries().data()),glidePlanesBoundaries().size()*sizeof(GlidePlaneBoundaryIO<dim>));

                infile.close();
                model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
                model::cout<<"  "<<glidePlanesBoundaries().size()<<" glidePlanesBoundaries"<<std::endl;

            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
            
        }

        
        /**********************************************************************/
        void readTxt(const size_t& runID)
        {
            const std::string filename(this->getTxtFilename(runID));
            
            std::ifstream infile (filename.c_str(), std::ios::in);
            if(infile.is_open())
            {
                const auto t0=std::chrono::system_clock::now();
                model::cout<<"reading "<<filename<<std::flush;
                
                size_t sizeGP;

                std::string line;
                std::stringstream ss;

                
                std::getline(infile, line);
                ss<<line;
                ss >> sizeGP;
                ss.clear();
                
                for(size_t k=0;k<sizeGP;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    glidePlanesBoundaries().emplace_back(ss);
                    ss.clear();
                }
                
                infile.close();
                model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
                model::cout<<"  "<<glidePlanesBoundaries().size()<<" glidePlanesBoundaries "<<std::endl;

            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
            
        }

        /**********************************************************************/
        void read(const size_t& runID)
        {
            if(isBinGood(runID))
            {
                readBin(runID);
            }
            else
            {
                if(isTxtGood(runID))
                {
                    readTxt(runID);
                }
                else
                {
                    std::cout<<"COULD NOT FIND INPUT FILEs evl/evl_"<<runID<<".bin or evl/evl_"<<runID<<".txt"<<std::endl;
                    assert(0 && "COULD NOT FIND INPUT FILEs.");
                }
            }
        }
        
    };
    
}
#endif

