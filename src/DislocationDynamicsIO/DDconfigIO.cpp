/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2020 by Danny Perez <danny_perez@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDconfigIO_cpp_
#define model_DDconfigIO_cpp_

#include <type_traits>
#include <vector>
#include <map>
#include <chrono>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <cfloat>      // std::ifstream
#include <DDbaseIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopNodeIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationNodeIO.h>
#include <DislocationSegmentIO.h>
#include <EshelbyInclusionIO.h>
//#include <DislocationNetwork.h>
//#include <PeriodicLoopIO.h>

#include <TerminalColors.h>
#include <DDconfigIO.h>


namespace model
{


    template <int dim>
    void DDconfigIO<dim>::make_maps()
    {
        
        // node map
        for(size_t k=0;k<loopNodes().size();++k)
        {
            loopNodeMap().insert(std::make_pair(loopNodes()[k].sID,k));
        }
        
        // node map
        for(size_t k=0;k<nodes().size();++k)
        {
            nodeMap().emplace(nodes()[k].sID,k);
        }
        
        // loop map
        for(size_t k=0;k<loops().size();++k)
        {
            loopMap().emplace(loops()[k].sID,k);
        }
    }



    template <int dim>
    DDconfigIO<dim>::DDconfigIO(const std::string& folderName,const std::string& suffix) :
    /* init */ DDbaseIO(folderName,"evl",suffix)
    {
        
    }


    template <int dim>
    const std::vector<DislocationNodeIO<dim>>& DDconfigIO<dim>::nodes() const
    {
        return *this;
    }

    template <int dim>
    std::vector<DislocationNodeIO<dim>>& DDconfigIO<dim>::nodes()
    {
        return *this;
    }

    template <int dim>
    const std::vector<DislocationLoopNodeIO<dim>>& DDconfigIO<dim>::loopNodes() const
    {
        return *this;
    }

    template <int dim>
    std::vector<DislocationLoopNodeIO<dim>>& DDconfigIO<dim>::loopNodes()
    {
        return *this;
    }

    template <int dim>
    const std::vector<DislocationLoopIO<dim>>& DDconfigIO<dim>::loops() const
    {
        return *this;
    }

    template <int dim>
    std::vector<DislocationLoopIO<dim>>& DDconfigIO<dim>::loops()
    {
        return *this;
    }



    template <int dim>
    const std::vector<DislocationLoopLinkIO<dim>>& DDconfigIO<dim>::loopLinks() const
    {
        return *this;
    }

    template <int dim>
    std::vector<DislocationLoopLinkIO<dim>>& DDconfigIO<dim>::loopLinks()
    {
        return *this;
    }

    template <int dim>
    const std::vector<EshelbyInclusionIO<dim>>& DDconfigIO<dim>::eshelbyInclusions() const
    {
        return *this;
    }

    template <int dim>
    std::vector<EshelbyInclusionIO<dim>>& DDconfigIO<dim>::eshelbyInclusions()
    {
        return *this;
    }


    template <int dim>
    const std::map<size_t,const size_t>& DDconfigIO<dim>::nodeMap() const
    {
        return _nodeMap;
    }

    template <int dim>
    std::map<size_t,const size_t>& DDconfigIO<dim>::nodeMap()
    {
        return _nodeMap;
    }

    template <int dim>
    const std::map<size_t,const size_t>& DDconfigIO<dim>::loopNodeMap() const
    {
        return _loopNodeMap;
    }

    template <int dim>
    std::map<size_t,const size_t>& DDconfigIO<dim>::loopNodeMap()
    {
        return _loopNodeMap;
    }

    template <int dim>
    const std::map<size_t, const size_t>& DDconfigIO<dim>::loopMap() const
    {
        return _loopMap;
    }

    template <int dim>
    std::map<size_t, const size_t>& DDconfigIO<dim>::loopMap()
    {
        return _loopMap;
    }

    template <int dim>
    std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>> DDconfigIO<dim>::segments() const
    {
        std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim> > temp;
        
        for(const auto& link : loopLinks())
        {
            
            
            const auto loopIter=loopMap().find(link.loopID);
            assert(loopIter!=loopMap().end());
            
            const auto loopSourceIter(loopNodeMap().find(link.sourceID));
            assert(loopSourceIter!=loopNodeMap().end());
            
            const auto loopSinkIter(loopNodeMap().find(link.sinkID));
            assert(loopSinkIter!=loopNodeMap().end());
            
            
            //                const bool arePeriodicBoundaryNodes((loopSourceIter->second->P-loopSinkIter->second->P).norm()<FLT_EPSILON
            //                                                    && loopSourceIter->second->edgeID>=0
            //                                                    &&   loopSinkIter->second->edgeID>=0
            //                                                    && loopSourceIter->second->loopID==loopSinkIter->second->loopID);
            
            if(link.hasNetworkLink)
            {
                const auto& loopNodeSource(loopNodes()[loopSourceIter->second]);
                const auto&   loopNodeSink(loopNodes()[  loopSinkIter->second]);
                const auto& loop(loops()[loopIter->second]);
                
                const size_t sourceID(std::min(loopNodeSource.networkNodeID,loopNodeSink.networkNodeID));
                const size_t   sinkID(std::max(loopNodeSource.networkNodeID,loopNodeSink.networkNodeID));
                const auto key=std::make_pair(sourceID,sinkID);
                
                const auto insertPair=temp.insert(std::make_pair(key,DislocationSegmentIO<dim>(sourceID,sinkID)));
                const auto& iter=insertPair.first;
                const bool& success=insertPair.second;
                
                // std::cout<<" Segment "<<sourceID<<"=>"<<sinkID;
                
                if(success)
                {
                    // std::cout<<" Adding link "<<link.sourceID<<"=>"<<link.sinkID<<" with normal "<<loopIter->second->N.transpose()<<"and burgres "<<loopIter->second->B.transpose();
                    iter->second.n=loop.N;
                    // std::cout<<" After adding Normal is "<<iter->second.n.transpose();
                    
                }
                else
                {
                    // std::cout<<" Adding link "<<link.sourceID<<"=>"<<link.sinkID<<" with normal "<<loopIter->second->N.transpose()<<"and burgres "<<loopIter->second->B.transpose();
                    
                    if(iter->second.n.cross(loop.N).norm()>FLT_EPSILON)
                    {
                        iter->second.n.setZero();
                    }
                    
                    // std::cout<<" After adding Normal is "<<iter->second.n.transpose();
                    
                }
                
                //                    if(link.sourceID<link.sinkID)
                if(loopNodeSource.networkNodeID<loopNodeSink.networkNodeID)
                {
                    iter->second.b+=loop.B;
                    // std::cout<<" Burgers vector "<<iter->second.b.transpose()<<std::endl;
                }
                else
                {
                    iter->second.b-=loop.B;
                    // std::cout<<" Burgers vector "<<iter->second.b.transpose()<<std::endl;
                    
                }
                
                if(iter->second.meshLocation==-1)
                {
                    iter->second.meshLocation=link.meshLocation;
                }
                else
                {
                    assert(iter->second.meshLocation==link.meshLocation);
                }
            }
        }
        
        return temp;
    }


    template <int dim>
    void DDconfigIO<dim>::write(const size_t& runID,const bool& outputBinary)
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

    template <int dim>
    void DDconfigIO<dim>::writeTxt(const size_t& runID)
    {
        const auto t0=std::chrono::system_clock::now();
        const std::string filename(this->getTxtFilename(runID));
        std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
        if(file.is_open())
        {
            std::cout<<"Writing "<<filename<<std::flush;
            writeTxtStream(file);
            file.close();
            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        else
        {
            throw std::runtime_error("Cannot open file"+filename);
//            std::cout<<"CANNOT OPEN "<<filename<<std::endl;
//            assert(false && "CANNOT OPEN FILE.");
        }
    }


    template <int dim>
    void DDconfigIO<dim>::writeTxtStream(std::ostream &file)
    {
        
        // Write header
        file<<nodes().size()<<"\n";
        file<<loops().size()<<"\n";
        file<<loopLinks().size()<<"\n";
        file<<loopNodes().size()<<"\n";
        file<<eshelbyInclusions().size()<<"\n";
        
        // Write Nodes
        for(const auto& node : nodes())
        {
            file<<node<<"\n";
        }
        // Write Loops
        for(const auto& loop : loops())
        {
            file<<loop<<"\n";
        }
        
        // Write Edges
        for(const auto& link : loopLinks())
        {
            file<<link<<"\n";
        }
        
        // Write LoopNodes
        for(const auto& loopNode : loopNodes())
        {
            file<<loopNode<<"\n";
        }
        
        // Eshelby inclusions
        for(const auto& ei : eshelbyInclusions())
        {
            file<<ei<<"\n";
        }
    }

    template <int dim>
    void DDconfigIO<dim>::writeBin(const size_t& runID)
    {
        const auto t0=std::chrono::system_clock::now();
        const std::string filename(this->getBinFilename(runID));
        std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
        if(file.is_open())
        {
            std::cout<<"Writing "<<filename<<std::flush;
            
            // Write header
            const size_t nV(nodes().size());
            const size_t nL(loops().size());
            const size_t nE(loopLinks().size());
            const size_t nLN(loopNodes().size());
            const size_t nEI(eshelbyInclusions().size());
            
            binWrite(file,nV);
            binWrite(file,nL);
            binWrite(file,nE);
            binWrite(file,nLN);
            binWrite(file,nEI);
            
            // Write Nodes
            for(const auto& node : nodes())
            {
                binWrite(file,node);
            }
            
            // Write Loops
            for(const auto& loop : loops())
            {
                binWrite(file,loop);
            }
            
            // Write Edges
            for(const auto& link : loopLinks())
            {
                binWrite(file,link);
            }
            
            // Write LoopNodes
            for(const auto& loopNode : loopNodes())
            {
                binWrite(file,loopNode);
            }
            
            // Write EshelbyInclusions
            for(const auto& ei : eshelbyInclusions())
            {
                binWrite(file,ei);
            }
            
            
            file.close();
            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
        }
        else
        {
            std::cout<<"CANNOT OPEN "<<filename<<std::endl;
            assert(false && "CANNOT OPEN FILE.");
        }
    }

    template <int dim>
    void DDconfigIO<dim>::read(const size_t& runID)
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
                throw std::runtime_error("No file evl_"+std::to_string(runID));
//                std::cout<<"COULD NOT FIND INPUT FILEs evl/evl_"<<runID<<".bin or evl/evl_"<<runID<<".txt"<<std::endl;
//                assert(0 && "COULD NOT FIND INPUT FILEs.");
            }
        }
    }

    template <int dim>
    void DDconfigIO<dim>::readBin(const size_t& runID)
    {
        
        nodes().clear();
        loops().clear();
        loopLinks().clear();
        loopNodes().clear();
        eshelbyInclusions().clear();
        nodeMap().clear();
        loopNodeMap().clear();
        loopMap().clear();
        
        
        const std::string filename(this->getBinFilename(runID));
        
        std::ifstream infile (filename.c_str(), std::ios::in|std::ios::binary);
        if(infile.is_open())
        {
            const auto t0=std::chrono::system_clock::now();
            std::cout<<"reading "<<filename<<std::flush;
            
            // Read header
            size_t sizeV;
            infile.read (reinterpret_cast<char*>(&sizeV), 1*sizeof(sizeV));
            size_t sizeL;
            infile.read (reinterpret_cast<char*>(&sizeL), 1*sizeof(sizeL));
            size_t sizeE;
            infile.read (reinterpret_cast<char*>(&sizeE), 1*sizeof(sizeE));
            size_t sizeLN;
            infile.read (reinterpret_cast<char*>(&sizeLN), 1*sizeof(sizeLN));
            size_t sizeEI;
            infile.read (reinterpret_cast<char*>(&sizeEI), 1*sizeof(sizeEI));
            
            // Read vertices
            nodes().resize(sizeV);
            infile.read (reinterpret_cast<char*>(nodes().data()),nodes().size()*sizeof(DislocationNodeIO<dim>));
            // Read loops
            loops().resize(sizeL);
            infile.read (reinterpret_cast<char*>(loops().data()),loops().size()*sizeof(DislocationLoopIO<dim>));
            // Read links
            loopLinks().resize(sizeE);
            infile.read (reinterpret_cast<char*>(loopLinks().data()),loopLinks().size()*sizeof(DislocationLoopLinkIO<dim>));
            // Read loopNodes
            loopNodes().resize(sizeLN);
            infile.read (reinterpret_cast<char*>(loopNodes().data()),loopNodes().size()*sizeof(DislocationLoopNodeIO<dim>));
            // Read Eshlby Inclusions
            eshelbyInclusions().resize(sizeEI);
            infile.read (reinterpret_cast<char*>(eshelbyInclusions().data()),eshelbyInclusions().size()*sizeof(EshelbyInclusionIO<dim>));
            
            
            infile.close();
            make_maps();
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            std::cout<<"  "<<nodes().size()<<" nodes "<<std::endl;
            std::cout<<"  "<<loops().size()<<" loops "<<std::endl;
            std::cout<<"  "<<loopLinks().size()<<" loopLinks "<<std::endl;
            std::cout<<"  "<<loopNodes().size()<<" loopNodes "<<std::endl;
            std::cout<<"  "<<eshelbyInclusions().size()<<" eshelbyInclusions "<<std::endl;
            
        }
        else
        {
            throw std::runtime_error("Cannot open file"+filename);
//
//            std::cout<<"CANNOT OPEN "<<filename<<std::endl;
//            assert(false && "CANNOT OPEN FILE.");
        }
        
    }


    template <int dim>
    void DDconfigIO<dim>::readTxt(const size_t& runID)
    {
        
        const std::string filename(this->getTxtFilename(runID));
        
        std::ifstream infile (filename.c_str(), std::ios::in);
        if(infile.is_open())
        {
            const auto t0=std::chrono::system_clock::now();
            std::cout<<"reading "<<filename<<std::flush;
            
            readTxtStream(infile);
            infile.close();
            
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            std::cout<<"  "<<nodes().size()<<" nodes "<<std::endl;
            std::cout<<"  "<<loops().size()<<" loops "<<std::endl;
            std::cout<<"  "<<loopLinks().size()<<" loopLinks "<<std::endl;
            std::cout<<"  "<<loopNodes().size()<<" loopNodes "<<std::endl;
            std::cout<<"  "<<eshelbyInclusions().size()<<" eshelbyInclusions "<<std::endl;
        }
        else
        {
            throw std::runtime_error("Cannot open file"+filename);
//
//            std::cout<<"CANNOT OPEN "<<filename<<std::endl;
//            assert(false && "CANNOT OPEN FILE.");
        }
        
    }

    template <int dim>
    void DDconfigIO<dim>::readTxtStream(std::istream &infile)
    {
        
        nodes().clear();
        loops().clear();
        loopLinks().clear();
        loopNodes().clear();
        eshelbyInclusions().clear();
        nodeMap().clear();
        loopNodeMap().clear();
        loopMap().clear();
        
        size_t sizeV;
        size_t sizeL;
        size_t sizeE;
        size_t sizeLN;
        size_t sizeEI;
        
        std::string line;
        std::stringstream ss;
        
        
        std::getline(infile, line);
        ss<<line;
        ss >> sizeV;
        ss.clear();
        
        std::getline(infile, line);
        ss<<line;
        ss >> sizeL;
        ss.clear();
        
        std::getline(infile, line);
        ss<<line;
        ss >> sizeE;
        ss.clear();
        
        std::getline(infile, line);
        ss<<line;
        ss >> sizeLN;
        ss.clear();
        
        std::getline(infile, line);
        ss<<line;
        ss >> sizeEI;
        ss.clear();
        
        
        nodes().clear();
        for(size_t k=0; k<sizeV; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            nodes().emplace_back(ss);
            ss.clear();
        }
        
        loops().clear();
        for(size_t k=0; k<sizeL; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            loops().emplace_back(ss);
            ss.clear();
        }
        
        loopLinks().clear();
        for(size_t k=0; k<sizeE; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            loopLinks().emplace_back(ss);
            ss.clear();
        }
        
        loopNodes().clear();
        for(size_t k=0; k<sizeLN; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            loopNodes().emplace_back(ss);
            ss.clear();
        }
        
        eshelbyInclusions().clear();
        for(size_t k=0; k<sizeEI; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            eshelbyInclusions().emplace_back(ss);
            ss.clear();
        }
        
        
        make_maps();
        
    }



    template <int dim>
    void DDconfigIO<dim>::bin2txt(const size_t& runID,const bool& writeSegments)
    {
        readBin(runID);
        writeTxt(runID);
        if(writeSegments)
        {    // std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>
            const std::string segmentsFilename(getTxtSegmentFilename(runID));
            const auto segmts(segments());
            std::ofstream segFile(segmentsFilename.c_str());
            for(const auto& seg : segmts)
            {
                segFile<<seg.second<<"\n";
            }
        }
    }

    template class DDconfigIO<3>;

}
#endif
