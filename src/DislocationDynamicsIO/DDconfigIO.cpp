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
#include <SphericalInclusionIO.h>
#include <PolyhedronInclusionIO.h>
#include <PolyhedronInclusionEdgeIO.h>
#include <PolyhedronInclusionNodeIO.h>

//#include <DislocationNetwork.h>
//#include <PeriodicLoopIO.h>

#include <TerminalColors.h>
#include <DDconfigIO.h>


namespace model
{


    template <int dim>
    void DDconfigIO<dim>::finalize()
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
void DDconfigIO<dim>::clear()
{
    nodes().clear();
    loops().clear();
    loopLinks().clear();
    loopNodes().clear();
    sphericalInclusions().clear();
    polyhedronInclusions().clear();
    polyhedronInclusionNodes().clear();
    polyhedronInclusionEdges().clear();
    nodeMap().clear();
    loopNodeMap().clear();
    loopMap().clear();
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
    const DislocationLoopIO<dim>& DDconfigIO<dim>::loop(const size_t& loopID) const
    {
        const auto loopIter(loopMap().find(loopID));
        if(loopIter!=loopMap().end())
        {
            return loops()[loopIter->second];
        }
        else
        {
            throw std::runtime_error("loopID not found");
            return loops()[0];
        }
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
    const std::vector<SphericalInclusionIO<dim>>& DDconfigIO<dim>::sphericalInclusions() const
    {
        return *this;
    }

    template <int dim>
    std::vector<SphericalInclusionIO<dim>>& DDconfigIO<dim>::sphericalInclusions()
    {
        return *this;
    }

template <int dim>
const std::vector<PolyhedronInclusionIO<dim> >& DDconfigIO<dim>::polyhedronInclusions() const
{
    return *this;
}

template <int dim>
std::vector<PolyhedronInclusionIO<dim> >& DDconfigIO<dim>::polyhedronInclusions()
{
    return *this;
}


template <int dim>
const std::vector<PolyhedronInclusionNodeIO<dim> >& DDconfigIO<dim>::polyhedronInclusionNodes() const
{
    return *this;
}

template <int dim>
std::vector<PolyhedronInclusionNodeIO<dim> >& DDconfigIO<dim>::polyhedronInclusionNodes()
{
    return *this;
}


template <int dim>
const std::vector<PolyhedronInclusionEdgeIO>& DDconfigIO<dim>::polyhedronInclusionEdges() const
{
    return *this;
}


template <int dim>
std::vector<PolyhedronInclusionEdgeIO>& DDconfigIO<dim>::polyhedronInclusionEdges()
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
        file<<sphericalInclusions().size()<<"\n";
        file<<polyhedronInclusions().size()<<"\n";
        file<<polyhedronInclusionNodes().size()<<"\n";
        file<<polyhedronInclusionEdges().size()<<"\n";

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
        for(const auto& ei : sphericalInclusions())
        {
            file<<ei<<"\n";
        }
        
        // polyhedronInclusions
        for(const auto& pi : polyhedronInclusions())
        {
            file<<pi<<"\n";
        }
        
        // polyhedronInclusionsNodes
        for(const auto& pin : polyhedronInclusionNodes())
        {
            file<<pin<<"\n";
        }
        
        // polyhedronInclusionsEdges
        for(const auto& pie : polyhedronInclusionEdges())
        {
            file<<pie<<"\n";
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
            const size_t nEI(sphericalInclusions().size());
            const size_t nPI(polyhedronInclusions().size());
            const size_t nPIn(polyhedronInclusionNodes().size());
            const size_t nPIe(polyhedronInclusionEdges().size());

            binWrite(file,nV);
            binWrite(file,nL);
            binWrite(file,nE);
            binWrite(file,nLN);
            binWrite(file,nEI);
            binWrite(file,nPI);
            binWrite(file,nPIn);
            binWrite(file,nPIe);


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
            
            // Write SphericalInclusions
            for(const auto& ei : sphericalInclusions())
            {
                binWrite(file,ei);
            }
            
            // polyhedronInclusions
            for(const auto& pi : polyhedronInclusions())
            {
                binWrite(file,pi);
            }
            
            // polyhedronInclusionsNodes
            for(const auto& pin : polyhedronInclusionNodes())
            {
                binWrite(file,pin);
            }
            
            // polyhedronInclusionsEdges
            for(const auto& pie : polyhedronInclusionEdges())
            {
                binWrite(file,pie);
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
        
//        nodes().clear();
//        loops().clear();
//        loopLinks().clear();
//        loopNodes().clear();
//        sphericalInclusions().clear();
//        polyhedronInclusions().clear();
//        polyhedronInclusionNodes().clear();
//        polyhedronInclusionEdges().clear();
//        nodeMap().clear();
//        loopNodeMap().clear();
//        loopMap().clear();
  
        clear();
        
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
            size_t sizePI;
            infile.read (reinterpret_cast<char*>(&sizePI), 1*sizeof(sizePI));
            size_t sizePIn;
            infile.read (reinterpret_cast<char*>(&sizePIn), 1*sizeof(sizePIn));
            size_t sizePIe;
            infile.read (reinterpret_cast<char*>(&sizePIe), 1*sizeof(sizePIe));

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
            sphericalInclusions().resize(sizeEI);
            infile.read (reinterpret_cast<char*>(sphericalInclusions().data()),sphericalInclusions().size()*sizeof(SphericalInclusionIO<dim>));
            polyhedronInclusions().resize(sizePI);
            infile.read (reinterpret_cast<char*>(polyhedronInclusions().data()),polyhedronInclusions().size()*sizeof(PolyhedronInclusionIO<dim>));
            polyhedronInclusionNodes().resize(sizePIn);
            infile.read (reinterpret_cast<char*>(polyhedronInclusionNodes().data()),polyhedronInclusionNodes().size()*sizeof(PolyhedronInclusionNodeIO<dim>));
            polyhedronInclusionEdges().resize(sizePIe);
            infile.read (reinterpret_cast<char*>(polyhedronInclusionEdges().data()),polyhedronInclusionEdges().size()*sizeof(PolyhedronInclusionEdgeIO));

            
            infile.close();
            finalize();
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;

            print();
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

            print();
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
void DDconfigIO<dim>::print() const
{
    std::cout<<"  "<<nodes().size()<<" nodes "<<std::endl;
    std::cout<<"  "<<loops().size()<<" loops "<<std::endl;
    std::cout<<"  "<<loopLinks().size()<<" loopLinks "<<std::endl;
    std::cout<<"  "<<loopNodes().size()<<" loopNodes "<<std::endl;
    std::cout<<"  "<<sphericalInclusions().size()<<" SphericalInclusions "<<std::endl;
    std::cout<<"  "<<polyhedronInclusions().size()<<" polyhedronInclusions "<<std::endl;
    std::cout<<"  "<<polyhedronInclusionNodes().size()<<" polyhedronInclusionNodes "<<std::endl;
    std::cout<<"  "<<polyhedronInclusionEdges().size()<<" polyhedronInclusionEdges "<<std::endl;
}

    template <int dim>
    void DDconfigIO<dim>::readTxtStream(std::istream &infile)
    {
        
//        nodes().clear();
//        loops().clear();
//        loopLinks().clear();
//        loopNodes().clear();
//        sphericalInclusions().clear();
//        polyhedronInclusions().clear();
//        polyhedronInclusionNodes().clear();
//        polyhedronInclusionEdges().clear();
//        nodeMap().clear();
//        loopNodeMap().clear();
//        loopMap().clear();
            
        clear();
        size_t sizeV;
        size_t sizeL;
        size_t sizeE;
        size_t sizeLN;
        size_t sizeEI;
        size_t sizePI;
        size_t sizePIn;
        size_t sizePIe;
        
        std::string line;
        std::stringstream ss;
        
        
        std::getline(infile, line);
        ss<<line;
        ss >> sizeV;
        ss.str("");
        ss.clear();
        
        std::getline(infile, line);
        ss<<line;
        ss >> sizeL;
        ss.str("");
        ss.clear();

        std::getline(infile, line);
        ss<<line;
        ss >> sizeE;
        ss.str("");
        ss.clear();

        std::getline(infile, line);
        ss<<line;
        ss >> sizeLN;
        ss.str("");
        ss.clear();

        std::getline(infile, line);
        ss<<line;
        ss >> sizeEI;
        ss.str("");
        ss.clear();
        
        std::getline(infile, line);
        ss<<line;
        ss >> sizePI;
        ss.str("");
        ss.clear();
        
        std::getline(infile, line);
        ss<<line;
        ss >> sizePIn;
        ss.str("");
        ss.clear();
        
        std::getline(infile, line);
        ss<<line;
        ss >> sizePIe;
        ss.str("");
        ss.clear();
        
        nodes().clear();
        for(size_t k=0; k<sizeV; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            nodes().emplace_back(ss);
            ss.str("");
            ss.clear();
        }
        
        loops().clear();
        for(size_t k=0; k<sizeL; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            loops().emplace_back(ss);
            ss.str("");
            ss.clear();
        }
        
        loopLinks().clear();
        for(size_t k=0; k<sizeE; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            loopLinks().emplace_back(ss);
            ss.str("");
            ss.clear();
        }
        
        loopNodes().clear();
        for(size_t k=0; k<sizeLN; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            loopNodes().emplace_back(ss);
            ss.str("");
            ss.clear();
        }
        
        sphericalInclusions().clear();
        for(size_t k=0; k<sizeEI; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            sphericalInclusions().emplace_back(ss);
            ss.str("");
            ss.clear();
        }
        
        polyhedronInclusions().clear();
        for(size_t k=0; k<sizePI; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            polyhedronInclusions().emplace_back(ss);
            ss.str("");
            ss.clear();
        }
        
        polyhedronInclusionNodes().clear();
        for(size_t k=0; k<sizePIn; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            polyhedronInclusionNodes().emplace_back(ss);
            ss.str("");
            ss.clear();
        }
        
        polyhedronInclusionEdges().clear();
        for(size_t k=0; k<sizePIe; ++k)
        {
            std::getline(infile, line);
            ss<<line;
            polyhedronInclusionEdges().emplace_back(ss);
            ss.str("");
            ss.clear();
        }
        
        
        finalize();
        
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

template <int dim>
std::map<size_t,std::vector<size_t>> DDconfigIO<dim>::loopNodeSequence() const
{
    
    std::map<size_t,std::vector<size_t>> temp;
    
    std::map<size_t,std::map<size_t,size_t>> loopLinkMap;
    for(const auto& link : loopLinks())
    {
        loopLinkMap[link.loopID].emplace(link.sourceID,link.sinkID);
    }
    
    
    for(const auto& pair : loopLinkMap)
    {
        if(pair.second.size())
        {
            size_t sourceID(pair.second.begin()->first);
            size_t sinkID(pair.second.begin()->second);
            for(size_t k=0;k<pair.second.size();++k)
            {
                
                temp[pair.first].push_back(sourceID);
                
                const auto mapIter(pair.second.find(sinkID));
                if(mapIter==pair.second.end())
                {
                    throw std::runtime_error("Could not create loopNodeSequence");
                }
                
                sourceID=mapIter->first;
                sinkID=mapIter->second;
            }

            if(sourceID!=temp[pair.first].front())
            {
                throw std::runtime_error("loopNodeSequence is not a loop");
            }
            
        }
    }
    
    return temp;
}


    template class DDconfigIO<3>;

}
#endif
