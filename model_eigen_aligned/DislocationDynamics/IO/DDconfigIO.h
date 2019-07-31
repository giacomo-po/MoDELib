/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDconfigIO_H_
#define model_DDconfigIO_H_

#include <vector>
#include <map>
#include <chrono>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <cfloat>      // std::ifstream
#include <DDbaseIO.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationEdgeIO.h>
#include <DislocationSegmentIO.h>
#include <MPIcout.h>


namespace model
{
    
    
    template <int dim>
    class DDconfigIO : public DDbaseIO
    /*              */,private std::vector<DislocationNodeIO<dim>>
    /*              */,private std::vector<DislocationLoopIO<dim>>
    /*              */,private std::vector<DislocationEdgeIO<dim>>
    /*              */,private std::map<size_t,const DislocationNodeIO<dim>* const>
    /*              */,private std::map<size_t, const DislocationLoopIO<dim>* const>
    {
        
        /**********************************************************************/
        void make_maps()
        {
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
        
        DDconfigIO(const std::string& suffix="") :
        /* init */ DDbaseIO("evl","evl",suffix)
        {
            
        }
        
        /**********************************************************************/
        const std::vector<DislocationNodeIO<dim>>& nodes() const
        {
            return *this;
        }
        
        std::vector<DislocationNodeIO<dim>>& nodes()
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
        const std::vector<DislocationLoopIO<dim>>& loops() const
        {
            return *this;
        }
        
        std::vector<DislocationLoopIO<dim>>& loops()
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
        const std::vector<DislocationEdgeIO<dim>>& links() const
        {
            return *this;
        }
        
        std::vector<DislocationEdgeIO<dim>>& links()
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>> segments() const
        {
            
            std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>> temp;
            
            for(const auto& link : links())
            {
            
                
                const auto loopIter=loopMap().find(link.loopID);
                assert(loopIter!=loopMap().end());
                
                const size_t sourceID(std::min(link.sourceID,link.sinkID));
                const size_t sinkID(std::max(link.sourceID,link.sinkID));
                const auto key=std::make_pair(sourceID,sinkID);
                
                const auto insertPair=temp.insert(std::make_pair(key,DislocationSegmentIO<dim>(sourceID,sinkID)));
                const auto& iter=insertPair.first;
                const bool& success=insertPair.second;
                
                if(success)
                {
                    iter->second.n=loopIter->second->N;
                }
                else
                {
                    if(iter->second.n.cross(loopIter->second->N).norm()>FLT_EPSILON)
                    {
                        iter->second.n.setZero();
                    }
                }
                
                if(link.sourceID<link.sinkID)
                {
                    iter->second.b+=loopIter->second->B;
                }
                else
                {
                    iter->second.b-=loopIter->second->B;
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
            
            return temp;
        }
        

        
        /**********************************************************************/
        template<typename DislocationNodeType>
        void writeTxt(const DislocationNodeType& dn,
                             const long int& runID)
        {
            
            
            const std::string filename(this->getTxtFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
            if(file.is_open())
            {
                // Write header
                file<<dn.nodes().size()<<"\n";
                file<<dn.loops().size()<<"\n";
                file<<dn.loopLinks().size()<<"\n";
                
                // Write Nodes
                for(const auto& node : dn.nodes())
                {
                    file<<DislocationNodeIO<dim>(*node.second)<<"\n";
                }
                
                // Write Loops
                for(const auto& loop : dn.loops())
                {
                    file<<DislocationLoopIO<dim>(*loop.second)<<"\n";
                }
                
                // Write Edges
                for(const auto& link : dn.loopLinks())
                {
                    file<<DislocationEdgeIO<dim>(link.second)<<"\n";
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
                             const std::vector<DislocationNodeIO<dim>> nodesIO,
                             const std::vector<DislocationLoopIO<dim>> loopsIO,
                             const std::vector<DislocationEdgeIO<dim>> edgesIO)
        {
            
            const std::string filename(this->getTxtFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
//            std::cout<<"Writing to "<<filename<<std::endl;
            if(file.is_open())
            {
                // Write header
                file<<nodesIO.size()<<"\n";
                file<<loopsIO.size()<<"\n";
                file<<edgesIO.size()<<"\n";
                
                // Write Nodes
                for(const auto& node : nodesIO)
                {
                    file<<node<<"\n";
                }
                // Write Loops
                for(const auto& loop : loopsIO)
                {
                    file<<loop<<"\n";
                }
                
                // Write Edges
                for(const auto& link : edgesIO)
                {
                    file<<link<<"\n";
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
        void bin2txt(const size_t& runID,const bool& writeSegments)
        {
            readBin(runID);
            writeTxt(runID,nodes(),loops(),links());
            if(writeSegments)
            {// std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>
                const std::string segmentsFilename(getTxtSegmentFilename(runID));
                const auto segmts(segments());
                std::ofstream segFile(segmentsFilename.c_str());
                for(const auto& seg : segmts)
                {
                    segFile<<seg.second<<"\n";
                }
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
                const size_t nV(dn.nodes().size());
                const size_t nL(dn.loops().size());
                const size_t nE(dn.loopLinks().size());
                binWrite(file,nV);
                binWrite(file,nL);
                binWrite(file,nE);
                
                // Write Nodes
                for(const auto& node : dn.nodes())
                {
                    binWrite(file,DislocationNodeIO<dim>(*node.second));
                }
                
                // Write Loops
                for(const auto& loop : dn.loops())
                {
                    binWrite(file,DislocationLoopIO<dim>(*loop.second));
                }
                
                // Write Edges
                for(const auto& link : dn.loopLinks())
                {
                    binWrite(file,DislocationEdgeIO<dim>(link.second));
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
                             const std::vector<DislocationNodeIO<dim>> nodesIO,
                             const std::vector<DislocationLoopIO<dim>> loopsIO,
                             const std::vector<DislocationEdgeIO<dim>> edgesIO)
        {
            
            const std::string filename(this->getBinFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
            if(file.is_open())
            {
                // Write header
                const size_t nV(nodesIO.size());
                const size_t nL(loopsIO.size());
                const size_t nE(edgesIO.size());
                binWrite(file,nV);
                binWrite(file,nL);
                binWrite(file,nE);
                
                // Write Nodes
                for(const auto& node : nodesIO)
                {
                    binWrite(file,node);
                }
                
                // Write Loops
                for(const auto& loop : loopsIO)
                {
                    binWrite(file,loop);
                }
                
                // Write Edges
                for(const auto& link : edgesIO)
                {
                    binWrite(file,link);
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
                size_t sizeV;
                infile.read (reinterpret_cast<char*>(&sizeV), 1*sizeof(sizeV));
                size_t sizeL;
                infile.read (reinterpret_cast<char*>(&sizeL), 1*sizeof(sizeL));
                size_t sizeE;
                infile.read (reinterpret_cast<char*>(&sizeE), 1*sizeof(sizeE));
                
                // Read vertices
                nodes().resize(sizeV);
                infile.read (reinterpret_cast<char*>(nodes().data()),nodes().size()*sizeof(DislocationNodeIO<dim>));
                // Read loops
                loops().resize(sizeL);
                infile.read (reinterpret_cast<char*>(loops().data()),loops().size()*sizeof(DislocationLoopIO<dim>));
                // Read links
                links().resize(sizeE);
                infile.read (reinterpret_cast<char*>(links().data()),links().size()*sizeof(DislocationEdgeIO<dim>));

                infile.close();
                make_maps();
                model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
                model::cout<<"  "<<nodes().size()<<" nodes "<<std::endl;
                model::cout<<"  "<<loops().size()<<" loops "<<std::endl;
                model::cout<<"  "<<links().size()<<" links "<<std::endl;

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
                
                size_t sizeV;
                size_t sizeL;
                size_t sizeE;

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
                
                for(size_t k=0;k<sizeV;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    nodes().emplace_back(ss);
                    ss.clear();
                }
                
                for(size_t k=0;k<sizeL;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    loops().emplace_back(ss);
                    ss.clear();
                }

                
                for(size_t k=0;k<sizeE;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    links().emplace_back(ss);
                    ss.clear();
                }

                infile.close();
                make_maps();
                model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
                model::cout<<"  "<<nodes().size()<<" nodes "<<std::endl;
                model::cout<<"  "<<loops().size()<<" loops "<<std::endl;
                model::cout<<"  "<<links().size()<<" links "<<std::endl;

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

//        /**********************************************************************/
//        template<typename T>
//        static void binWrite(std::ofstream& of,const T& o)
//        {
//            of.write((char *) &o, (sizeof o));
//        }
//
//        /**********************************************************************/
//        template<typename T>
//        static void binRead(std::ifstream& file,
//                            T*& memblock,
//                            const size_t& arraySize)
//        {
//            memblock = new T [arraySize];
//            file.read (reinterpret_cast<char*>(memblock), arraySize*sizeof (T));
//        }
//
//        /**********************************************************************/
//        static std::string getBinFilename(const size_t& runID,const std::string& suffix="")
//        {
//            return "evl"+suffix+"/evl_"+std::to_string(runID)+".bin";
//        }
//
//        /**********************************************************************/
//        static std::string getTxtFilename(const size_t& runID,const std::string& suffix="")
//        {
//            return "evl"+suffix+"/evl_"+std::to_string(runID)+".txt";
//        }
//
//        /**********************************************************************/
//        static std::string getTxtSegmentFilename(const size_t& runID,const std::string& suffix="")
//        {
//            return "evl"+suffix+"/S_"+std::to_string(runID)+".txt";
//        }
