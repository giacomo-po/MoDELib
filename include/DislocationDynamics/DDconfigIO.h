/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2020 by Danny Perez <danny_perez@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDconfigIO_H_
#define model_DDconfigIO_H_

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
//#include <PeriodicLoopIO.h>

#include <TerminalColors.h>
#include <set>


namespace model
{
    
    
    template <int dim>
    class DDconfigIO : public DDbaseIO
    /*              */,private std::vector<DislocationNodeIO<dim> >
    /*              */,private std::vector<DislocationLoopIO<dim> >
    /*              */,private std::vector<DislocationLoopNodeIO<dim> >
    /*              */,private std::vector<DislocationLoopLinkIO<dim> >
    /*              */,private std::map<size_t,const DislocationLoopNodeIO<dim>* const>
    /*              */,private std::map<size_t,const DislocationNodeIO<dim>* const>
    /*              */,private std::map<size_t, const DislocationLoopIO<dim>* const>
    {
        
        //        static_assert(std::is_pod<DislocationNodeIO<dim>>::value,"DislocationNodeIO<dim> is NOT PLANE OLD DATA");
        
        /**********************************************************************/
        void make_maps()
        {
            
            // node map
            for(const auto& loopNode : loopNodes())
            {
                loopNodeMap().emplace(loopNode.sID,&loopNode);
            }
            
            // node map
            for(const auto& node : nodes())
            {
                nodeMap().emplace(node.sID,&node);
            }
            
            // loop map
            for(const auto& loop : loops())
            {
                loopMap().emplace(loop.sID,&loop);
            }
        }
        
        
    public:
        
        /**********************************************************************/
        DDconfigIO(const std::string& suffix="") :
        /* init */ DDbaseIO("evl","evl",suffix)
        {
            
        }
        
        /**********************************************************************/
        template<typename DislocationNetworkType>
        DDconfigIO(const DislocationNetworkType& dn,
                   const std::string& suffix="") :
        /* init */ DDbaseIO("evl","evl",suffix)
        {
            
            
            // Write Loops
            for(const auto& loop : dn.loops())
            {
                loops().emplace_back(*loop.second.lock());
            }
            

            // Write LoopNodes
            for(const auto& node : dn.loopNodes())
            {
                loopNodes().emplace_back(*node.second.lock());
            }
            

            
            // Write LoopLinks
            for(const auto& link : dn.loopLinks())
            {
                loopLinks().emplace_back(link.second);
            }
            

            
            // Write NetworkNodes
            for(const auto& node : dn.networkNodes())
            {
                nodes().emplace_back(*node.second.lock());
            }
            

            
 
        }
        
        /**********************************************************************/
        const std::vector<DislocationNodeIO<dim> >& nodes() const
        {
            return *this;
        }
        
        std::vector<DislocationNodeIO<dim> >& nodes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::vector<DislocationLoopNodeIO<dim> >& loopNodes() const
        {
            return *this;
        }
        
        std::vector<DislocationLoopNodeIO<dim> >& loopNodes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::map<size_t,const DislocationNodeIO<dim>* const>& nodeMap() const
        {
            return *this;
        }
        
        std::map<size_t,const DislocationNodeIO<dim>* const>& nodeMap()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::map<size_t,const DislocationLoopNodeIO<dim>* const>& loopNodeMap() const
        {
            return *this;
        }
        
        std::map<size_t,const DislocationLoopNodeIO<dim>* const>& loopNodeMap()
        {
            return *this;
        }
        
        
        /**********************************************************************/
        const std::vector<DislocationLoopIO<dim> >& loops() const
        {
            return *this;
        }
        
        std::vector<DislocationLoopIO<dim> >& loops()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::map<size_t, const DislocationLoopIO<dim>* const>& loopMap() const
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<size_t, const DislocationLoopIO<dim>* const>& loopMap()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::vector<DislocationLoopLinkIO<dim> >& loopLinks() const
        {
            return *this;
        }
        
        std::vector<DislocationLoopLinkIO<dim> >& loopLinks()
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim> > segments() const
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
                    const size_t sourceID(std::min(loopSourceIter->second->networkNodeID,loopSinkIter->second->networkNodeID));
                    const size_t   sinkID(std::max(loopSourceIter->second->networkNodeID,loopSinkIter->second->networkNodeID));
                    const auto key=std::make_pair(sourceID,sinkID);
                    
                    const auto insertPair=temp.insert(std::make_pair(key,DislocationSegmentIO<dim>(sourceID,sinkID)));
                    const auto& iter=insertPair.first;
                    const bool& success=insertPair.second;
                    
                    // std::cout<<" Segment "<<sourceID<<"=>"<<sinkID;

                    if(success)
                    {
                        // std::cout<<" Adding link "<<link.sourceID<<"=>"<<link.sinkID<<" with normal "<<loopIter->second->N.transpose()<<"and burgres "<<loopIter->second->B.transpose();
                        iter->second.n=loopIter->second->N;
                        // std::cout<<" After adding Normal is "<<iter->second.n.transpose();

                    }
                    else
                    {
                        // std::cout<<" Adding link "<<link.sourceID<<"=>"<<link.sinkID<<" with normal "<<loopIter->second->N.transpose()<<"and burgres "<<loopIter->second->B.transpose();

                        if(iter->second.n.cross(loopIter->second->N).norm()>FLT_EPSILON)
                        {
                            iter->second.n.setZero();
                        }

                        // std::cout<<" After adding Normal is "<<iter->second.n.transpose();

                    }
                    
//                    if(link.sourceID<link.sinkID)
                    if(loopSourceIter->second->networkNodeID<loopSinkIter->second->networkNodeID)
                    {
                        iter->second.b+=loopIter->second->B;
                        // std::cout<<" Burgers vector "<<iter->second.b.transpose()<<std::endl;
                    }
                    else
                    {
                        iter->second.b-=loopIter->second->B;
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

        std::map<std::pair<size_t,size_t>,std::set<size_t> > segmentloopMap() const
        {
            
           std::map<std::pair<size_t,size_t>,std::set<size_t> >temp;
            
            for(const auto& link : loopLinks())
            {
                
                
                const auto loopIter=loopMap().find(link.loopID);
                assert(loopIter!=loopMap().end());
                
                const auto loopSourceIter(loopNodeMap().find(link.sourceID));
                assert(loopSourceIter!=loopNodeMap().end());

                const auto loopSinkIter(loopNodeMap().find(link.sinkID));
                assert(loopSinkIter!=loopNodeMap().end());

                
                if(link.hasNetworkLink)
                {
                    const size_t sourceID(std::min(loopSourceIter->second->networkNodeID,loopSinkIter->second->networkNodeID));
                    const size_t   sinkID(std::max(loopSourceIter->second->networkNodeID,loopSinkIter->second->networkNodeID));
                    const auto key=std::make_pair(sourceID,sinkID);

                    temp[key].insert(link.loopID);
                    
                }
            }
            
            return temp;
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
        void writeTxt(const size_t& runID)
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
                std::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
        }
        
        
        /**********************************************************************/
        void writeTxtStream(std::ostream &file)
        {
            
            // Write header
            file<<nodes().size()<<"\n";
            file<<loops().size()<<"\n";
            file<<loopLinks().size()<<"\n";
            file<<loopNodes().size()<<"\n";

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
            
            
            
            
//            std::cout<<" WRITING:  "<<nodes().size()<<" nodes "<<std::endl;
//            std::cout<<" WRITING:  "<<loops().size()<<" loops "<<std::endl;
//            std::cout<<" WRITING:  "<<loopLinks().size()<<" links "<<std::endl;
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
                const size_t nV(nodes().size());
                const size_t nL(loops().size());
                const size_t nE(loopLinks().size());
                const size_t nLN(loopNodes().size());

                binWrite(file,nV);
                binWrite(file,nL);
                binWrite(file,nE);
                binWrite(file,nLN);

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
                
                
                file.close();
                std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
                
            }
            else
            {
                std::cout<<"CANNOT OPEN "<<filename<<std::endl;
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
        
        /**********************************************************************/
        void readBin(const size_t& runID)
        {
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

                
                infile.close();
                make_maps();
                std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
                std::cout<<"  "<<nodes().size()<<" nodes "<<std::endl;
                std::cout<<"  "<<loops().size()<<" loops "<<std::endl;
                std::cout<<"  "<<loopLinks().size()<<" loopLinks "<<std::endl;
                std::cout<<"  "<<loopNodes().size()<<" loopNodes "<<std::endl;

            }
            else
            {
                std::cout<<"CANNOT OPEN "<<filename<<std::endl;
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
                std::cout<<"reading "<<filename<<std::flush;
                
                readTxtStream(infile);
                infile.close();
                
                std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
                std::cout<<"  "<<nodes().size()<<" nodes "<<std::endl;
                std::cout<<"  "<<loops().size()<<" loops "<<std::endl;
                std::cout<<"  "<<loopLinks().size()<<" loopLinks "<<std::endl;
                std::cout<<"  "<<loopNodes().size()<<" loopNodes "<<std::endl;
            }
            else
            {
                std::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
            
        }
        
        
        
        
        /**********************************************************************/
        void readTxtStream(std::istream &infile)
        {
            
            size_t sizeV;
            size_t sizeL;
            size_t sizeE;
            size_t sizeLN;

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
            
            
            make_maps();
            
            std::cout<<" READING:  "<<nodes().size()<<" nodes "<<std::endl;
            std::cout<<" READING:  "<<loops().size()<<" loops "<<std::endl;
            std::cout<<" READING:  "<<loopLinks().size()<<" links "<<std::endl;
            
        }
        
        
        
        /**********************************************************************/
        void bin2txt(const size_t& runID,const bool& writeSegments)
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
        
    };
    
}
#endif
