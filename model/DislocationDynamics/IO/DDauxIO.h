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
#include <PeriodicGlidePlane.h>
#include <DislocationQuadraturePointIO.h>
#include <PeriodicLoopLinkIO.h>
#include <PeriodicLoopNodeIO.h>

namespace model
{
    
    
    template <int dim>
    struct DDauxIO : public DDbaseIO
    /*            */,private std::vector<GlidePlaneBoundaryIO<dim>>
    /*            */,private std::vector<PeriodicPlanePatchIO<dim>>
    /*            */,private std::vector<DislocationQuadraturePointIO<dim>>
    /*            */,private std::vector<PeriodicLoopNodeIO<dim>>
    /*            */,private std::vector<PeriodicLoopLinkIO<dim>>
    {
        
        /**********************************************************************/
        DDauxIO(const std::string& suffix="") :
        /* init */ DDbaseIO("evl","ddAux",suffix)
        {
            
        }
        
        /**********************************************************************/
        template<typename DislocationNodeType>
        DDauxIO(const DislocationNodeType& DN,
                   const std::string& suffix="") :
        /* init */ DDbaseIO("evl","ddAux",suffix)
        {

            if (DN.outputQuadraturePoints)
            {
                for (const auto& link : DN.links())
                {
                    for(const auto& qPoint : link.second->quadraturePoints())
                    {
                        quadraturePoints().push_back(qPoint);
                    }
                }
            }
            
            if(DN.outputGlidePlanes)
            {
                setGlidePlaneBoundaries(DN.glidePlaneFactory);
            }
            
            if(DN.outputPeriodicConfiguration && DN.periodicDislocationLoopFactory)
            {
                for(const auto& pair : *DN.periodicDislocationLoopFactory)
                {// output periodic glide planes too
                    
                    if(!pair.second.expired())
                    {
                        const auto periodicLoop(pair.second.lock());
                        addPeriodicGlidePlane(*periodicLoop->periodicGlidePlane);
                        
                        for(const auto& node : periodicLoop->nodes())
                        {
                            if(!node.second.expired())
                            {
                                periodicLoopNodes().emplace_back(*node.second.lock());
                            }
                        }
                        
                        for(const auto& link : periodicLoop->loopLinks())
                        {
                            periodicLoopLinks().emplace_back(link.second);
                        }
                        
                    }
                }
            }
            
        }
        
        /**********************************************************************/
        void setGlidePlaneBoundaries(const GlidePlaneFactory<dim>& gpf)
        {
            glidePlanesBoundaries().clear();
            for(const auto& pair : gpf.glidePlanes())
            {
                const auto glidePlane(pair.second.lock());
                if(glidePlane)
                {
                    for(const auto& seg : glidePlane->meshIntersections)
                    {
                        glidePlanesBoundaries().emplace_back(glidePlane->sID,*seg);
                    }
                }
            }
        }
        
        /**********************************************************************/
        void addPeriodicGlidePlane(const PeriodicGlidePlane<dim>& pgp)
        {
            for(const auto& patch : pgp.patches())
            {
                periodicGlidePlanePatches().emplace_back(*patch.second);
            }
        }
        
//        /**********************************************************************/
//        void addDislocationQuadraturePoint(const DislocationQuadraturePointIO<dim>& qPoint)
//        {
//            quadraturePoints().push_back(qPoint);
//        }
        
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
        const std::vector<PeriodicPlanePatchIO<dim>>& periodicGlidePlanePatches() const
        {
            return *this;
        }
        
        std::vector<PeriodicPlanePatchIO<dim>>& periodicGlidePlanePatches()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::vector<DislocationQuadraturePointIO<dim>>& quadraturePoints() const
        {
            return *this;
        }
        
        std::vector<DislocationQuadraturePointIO<dim>>& quadraturePoints()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::vector<PeriodicLoopNodeIO<dim>>& periodicLoopNodes() const
        {
            return *this;
        }
        
        std::vector<PeriodicLoopNodeIO<dim>>& periodicLoopNodes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::vector<PeriodicLoopLinkIO<dim>>& periodicLoopLinks() const
        {
            return *this;
        }
        
        std::vector<PeriodicLoopLinkIO<dim>>& periodicLoopLinks()
        {
            return *this;
        }
        
        /**********************************************************************/
        void write(const size_t& runID,const bool& outputBinary)
        {
            if(outputBinary)
            {
                writeBin(runID);
            }
            else
            {
                writeTxt(runID);
            }
        }
        
        /**********************************************************************/
        void writeTxt(const long int& runID)
        {
            const auto t0=std::chrono::system_clock::now();
            const std::string filename(this->getTxtFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
            if(file.is_open())
            {
                std::cout<<"Writing "<<filename<<std::flush;
                // Write header
                file<<quadraturePoints().size()<<"\n";
                file<<glidePlanesBoundaries().size()<<"\n";
                file<<periodicGlidePlanePatches().size()<<"\n";
                file<<periodicLoopNodes().size()<<"\n";
                file<<periodicLoopLinks().size()<<"\n";

                // Write body
                for(const auto& qPoint : quadraturePoints())
                {
                    file<<qPoint<<"\n";
                }
                for(const auto& gpBnd : glidePlanesBoundaries())
                {
                    file<<gpBnd<<"\n";
                }
                for(const auto& patch : periodicGlidePlanePatches())
                {
                    file<<patch<<"\n";
                }
                for(const auto& periodicNode : periodicLoopNodes())
                {
                    file<<periodicNode<<"\n";
                }
                for(const auto& periodicLink : periodicLoopLinks())
                {
                    file<<periodicLink<<"\n";
                }

                file.close();
                std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
        }
        
        /**********************************************************************/
        void writeBin(const size_t& runID)
        {
            const auto t0=std::chrono::system_clock::now();
            const std::string filename(this->getBinFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
            if(file.is_open())
            {
                std::cout<<"Writing "<<filename<<std::flush;
                // Write header
                binWrite(file,quadraturePoints().size());
                binWrite(file,glidePlanesBoundaries().size());
                binWrite(file,periodicGlidePlanePatches().size());
                binWrite(file,periodicLoopNodes().size());
                binWrite(file,periodicLoopLinks().size());

                // Write body
                for(const auto& qPoint : quadraturePoints())
                {
                    binWrite(file,qPoint);
                }
                for(const auto& gpb : glidePlanesBoundaries())
                {
                    binWrite(file,gpb);
                }
                for(const auto& patch : periodicGlidePlanePatches())
                {
                    binWrite(file,patch);
                }
                for(const auto& node : periodicLoopNodes())
                {
                    binWrite(file,node);
                }
                for(const auto& link : periodicLoopLinks())
                {
                    binWrite(file,link);
                }
                
                file.close();
                std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
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
                size_t sizeQP;
                infile.read (reinterpret_cast<char*>(&sizeQP), 1*sizeof(sizeQP));
                size_t sizeGP;
                infile.read (reinterpret_cast<char*>(&sizeGP), 1*sizeof(sizeGP));
                size_t sizePPP;
                infile.read (reinterpret_cast<char*>(&sizePPP), 1*sizeof(sizePPP));
                size_t sizePLN;
                infile.read (reinterpret_cast<char*>(&sizePLN), 1*sizeof(sizePLN));
                size_t sizePLL;
                infile.read (reinterpret_cast<char*>(&sizePLL), 1*sizeof(sizePLL));

                // Read body
                quadraturePoints().resize(sizeQP);
                infile.read (reinterpret_cast<char*>(quadraturePoints().data()),quadraturePoints().size()*sizeof(DislocationQuadraturePointIO<dim>));
                glidePlanesBoundaries().resize(sizeGP);
                infile.read (reinterpret_cast<char*>(glidePlanesBoundaries().data()),glidePlanesBoundaries().size()*sizeof(GlidePlaneBoundaryIO<dim>));
                periodicGlidePlanePatches().resize(sizePPP);
                infile.read (reinterpret_cast<char*>(periodicGlidePlanePatches().data()),periodicGlidePlanePatches().size()*sizeof(PeriodicPlanePatchIO<dim>));
                periodicLoopNodes().resize(sizePLN);
                infile.read (reinterpret_cast<char*>(periodicLoopNodes().data()),periodicLoopNodes().size()*sizeof(PeriodicLoopNodeIO<dim>));
                periodicLoopLinks().resize(sizePLL);
                infile.read (reinterpret_cast<char*>(periodicLoopLinks().data()),periodicLoopLinks().size()*sizeof(PeriodicLoopLinkIO<dim>));

                infile.close();
                printLog(t0);
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
                std::string line;
                std::stringstream ss;

                size_t sizeQP;
                std::getline(infile, line);
                ss<<line;
                ss >> sizeQP;
                ss.clear();
                
                size_t sizeGP;
                std::getline(infile, line);
                ss<<line;
                ss >> sizeGP;
                ss.clear();

                size_t sizePPP;
                std::getline(infile, line);
                ss<<line;
                ss >> sizePPP;
                ss.clear();
                
                size_t sizePLN;
                std::getline(infile, line);
                ss<<line;
                ss >> sizePLN;
                ss.clear();
                
                size_t sizePLL;
                std::getline(infile, line);
                ss<<line;
                ss >> sizePLL;
                ss.clear();
                
                quadraturePoints().clear();
                quadraturePoints().reserve(sizeQP);
                for(size_t k=0;k<sizeQP;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    quadraturePoints().emplace_back(ss);
                    ss.clear();
                }
                
                glidePlanesBoundaries().clear();
                glidePlanesBoundaries().reserve(sizeGP);
                for(size_t k=0;k<sizeGP;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    glidePlanesBoundaries().emplace_back(ss);
                    ss.clear();
                }
                
                periodicGlidePlanePatches().clear();
                periodicGlidePlanePatches().reserve(sizePPP);
                for(size_t k=0;k<sizePPP;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    periodicGlidePlanePatches().emplace_back(ss);
                    ss.clear();
                }
                
                periodicLoopNodes().clear();
                periodicLoopNodes().reserve(sizePLN);
                for(size_t k=0;k<sizePLN;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    periodicLoopNodes().emplace_back(ss);
                    ss.clear();
                }
                
                periodicLoopLinks().clear();
                periodicLoopLinks().reserve(sizePLL);
                for(size_t k=0;k<sizePLL;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    periodicLoopLinks().emplace_back(ss);
                    ss.clear();
                }

                
                infile.close();
                printLog(t0);
            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
            
        }
        
        void printLog(const std::chrono::time_point<std::chrono::system_clock>& t0) const
        {
            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            model::cout<<"  "<<glidePlanesBoundaries().size()<<" glidePlanesBoundaries "<<std::endl;
            model::cout<<"  "<<periodicGlidePlanePatches().size()<<" periodicPlanePatches"<<std::endl;
            model::cout<<"  "<<quadraturePoints().size()<<" quadraturePoints"<<std::endl;
            model::cout<<"  "<<periodicLoopNodes().size()<<" periodicLoopNodes"<<std::endl;
            model::cout<<"  "<<periodicLoopLinks().size()<<" periodicLoopLinks"<<std::endl;

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


//        /**********************************************************************/
//        void writeTxt(const size_t& runID,
//                             const std::vector<GlidePlaneBoundaryIO<dim>> gpBoundaries)
//        {
//
//            const std::string filename(this->getTxtFilename(runID));
//            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
////            std::cout<<"Writing to "<<filename<<std::endl;
//            if(file.is_open())
//            {
//                // Write header
//                file<<gpBoundaries.size()<<"\n";
//
//                // Write Nodes
//                for(const auto& gpBnd : gpBoundaries)
//                {
//                    file<<gpBnd<<"\n";
//                }
//
//                file.close();
//            }
//            else
//            {
//                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
//                assert(false && "CANNOT OPEN FILE.");
//            }
//        }

//        /**********************************************************************/
//        template<typename DislocationNodeType>
//        void writeBin(const DislocationNodeType& dn,
//                             const long int& runID)
//        {
//
//            const std::string filename(this->getBinFilename(runID));
//            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
//            if(file.is_open())
//            {
//                // Write header
//                size_t nGP(0);
//                for(const auto& glidePlane : dn.glidePlanes())
//                {
//                    nGP+=glidePlane.second->meshIntersections.size();
//                }
//                binWrite(file,nGP);
//
//                // Write GlidePlaneBoundaries
//                for(const auto& glidePlane : dn.glidePlanes())
//                {
//                    for(const auto& seg : glidePlane->second->meshIntersections)
//                    {
//                        binWrite(file,GlidePlaneBoundaryIO<dim>(seg));
//                    }
//                }
//
//                file.close();
//            }
//            else
//            {
//                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
//                assert(false && "CANNOT OPEN FILE.");
//            }
//        }
