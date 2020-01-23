
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <StaticID.h>
#include <CompareVectorsByComponent.h>
#include <Polygon2D.h>
#include <SegmentSegmentDistance.h>
#include <stdio.h>



using namespace model;

std::vector<Eigen::Matrix<double,2,1>> readPolyVector(const std::string& fileName)
{
    std::vector<Eigen::Matrix<double,2,1>> temp;
    std::ifstream file ( fileName.c_str() , std::ifstream::in );
    if(file.is_open())
    {
        
        std::string line;
        double x,y;
        
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            
            ss >> x >>y;
            temp.push_back((Eigen::Matrix<double,2,1>()<<x,y).finished());
            //std::cout<<x<<","<<y<<std::endl;
            
            //poly[0] <<ClipperLib::DoublePoint(x,y);
        }
        
    }
    else
    {
        std::cout<<"CANNOT READ "+fileName<<std::endl;
    }
    
    
    
    return temp;
}


struct LoopPathClipperNode;

struct LoopEdge
{
    
    const std::shared_ptr<LoopPathClipperNode> source;
    const std::shared_ptr<LoopPathClipperNode> sink;
    const int wn;
    
    LoopEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
             const std::shared_ptr<LoopPathClipperNode>& sink_in,
             const int& wn_in);
    
};

struct PatchEdge
{
    
    const std::shared_ptr<LoopPathClipperNode> source;
    const std::shared_ptr<LoopPathClipperNode> sink;
    const int wn;

    PatchEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
              const std::shared_ptr<LoopPathClipperNode>& sink_in,
              const int& wn_in);
    
};

struct LoopPathClipperNode : public Eigen::Matrix<double,2,1>,
                             public StaticID<LoopPathClipperNode>
{
    
    const LoopEdge* loopPrev;
    const LoopEdge* loopNext;
    const PatchEdge* patchPrev;
    const PatchEdge* patchNext;
    
    LoopPathClipperNode(const Eigen::Matrix<double,2,1>& pos) :
    /* init */ Eigen::Matrix<double,2,1>(pos)
    /* init */,loopPrev(nullptr)
    /* init */,loopNext(nullptr)
    /* init */,patchPrev(nullptr)
    /* init */,patchNext(nullptr)
    {
        
    }
    
};

LoopEdge::LoopEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
         const std::shared_ptr<LoopPathClipperNode>& sink_in,
         const int& wn_in) :
/* init */ source(source_in)
/* init */,sink(sink_in)
/* init */,wn(wn_in)
{
    assert(source->loopNext==nullptr);
    source->loopNext=this;
    assert(sink->loopPrev==nullptr);
    sink->loopPrev=this;
}

PatchEdge::PatchEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
                     const std::shared_ptr<LoopPathClipperNode>& sink_in,
                     const int& wn_in) :
/* init */ source(source_in)
/* init */,sink(sink_in)
/* init */,wn(wn_in)
{
    assert(source->patchNext==nullptr);
    source->patchNext=this;
    assert(sink->patchPrev==nullptr);
    sink->patchPrev=this;
    
    
}

struct LoopPathClipper : private std::map<Eigen::Matrix<double,2,1>,const std::weak_ptr<LoopPathClipperNode>,CompareVectorsByComponent<double,2,float>>
/*                    */,public std::map<std::pair<size_t,size_t>,LoopEdge>
/*                    */,public std::map<std::pair<size_t,size_t>,PatchEdge>

{
    typedef Eigen::Matrix<double,2,1> VectorDim;
    typedef std::map<Eigen::Matrix<double,2,1>,const std::weak_ptr<LoopPathClipperNode>,CompareVectorsByComponent<double,2,float>> NodeCointainerType;
    
//    const std::vector<VectorDim>& loop;
//    const std::vector<VectorDim>& patch;
    

    
    /**********************************************************************/
    NodeCointainerType& nodes()
    {
        return *this;
    }
    
    const NodeCointainerType& nodes() const
    {
        return *this;
    }
    
    const std::map<std::pair<size_t,size_t>,LoopEdge>& loopEdges() const
    {
        return *this;
    }

    std::map<std::pair<size_t,size_t>,LoopEdge>& loopEdges()
    {
        return *this;
    }
    
    const std::map<std::pair<size_t,size_t>,PatchEdge>& patchEdges() const
    {
        return *this;
    }
    
    std::map<std::pair<size_t,size_t>,PatchEdge>& patchEdges()
    {
        return *this;
    }

    /**********************************************************************/
    std::shared_ptr<LoopPathClipperNode> getSharedNode(const VectorDim& point)
    {
        const auto iter(nodes().find(point));
        if(iter==nodes().end())
        {// point does not exist
            return nodes().emplace(point,std::shared_ptr<LoopPathClipperNode >(new LoopPathClipperNode(point))).first->second.lock();
        }
        else
        {// point exists
            if(iter->second.expired())
            {// node deleted elsewhere
                nodes().erase(iter);
                return nodes().emplace(point,std::shared_ptr<LoopPathClipperNode>(new LoopPathClipperNode(point))).first->second.lock();
            }
            else
            {
                return iter->second.lock();
            }
        }
    }
    
    LoopPathClipper(const std::vector<VectorDim>& loop,
                    const std::vector<VectorDim>& patch)
//     loop(loop_in)
//    ,patch(patch_in)
    {
        std::vector<std::map<double,VectorDim>> loopIntersections;
        std::vector<std::map<double,VectorDim>> patchIntersections;
        loopIntersections.resize(loop.size());
        patchIntersections.resize(patch.size());

        for(size_t i=0;i<loop.size();++i)
        {
            const size_t i1(i+1<loop.size()? i+1 :0);
            loopIntersections[i].emplace(0.0,loop[i]);
            loopIntersections[i].emplace(1.0,loop[i1]);
        }
        
        for(size_t j=0;j<patch.size();++j)
        {
            const size_t j1(j+1<patch.size()? j+1 :0);
            patchIntersections[j].emplace(0.0,patch[j]);
            patchIntersections[j].emplace(1.0,patch[j1]);
        }
        
        for(size_t i=0;i<loop.size();++i)
        {
            const size_t i1(i+1<loop.size()? i+1 :0);
            for(size_t j=0;j<patch.size();++j)
            {
                const size_t j1(j+1<patch.size()? j+1 :0);
                SegmentSegmentDistance<2> ssd(loop[i],loop[i1],patch[j],patch[j1]);
                for(const auto& inter : ssd.intersectionSegment())
                {
                    loopIntersections[i].emplace(std::get<1>(inter),std::get<0>(inter));
                    patchIntersections[j].emplace(std::get<2>(inter),std::get<0>(inter));
                }

            }
        }
        
        for(const auto& map_val : loopIntersections)
        {
            
            for(auto iter=map_val.begin();iter!=map_val.end();++iter)
            {
                auto iter1(iter);
                iter1++;
                if(iter1!=map_val.end())
                {
                    const auto source(getSharedNode( iter->second));
                    const auto   sink(getSharedNode(iter1->second));
                    const VectorDim midPoint(0.5*(*source+*sink));
                    loopEdges().emplace(std::make_pair(source->sID,sink->sID),LoopEdge(source,sink,Polygon2D::windingNumber(midPoint,patch)));
                }
            }

        }
        
        for(const auto& map_val : patchIntersections)
        {
            
            for(auto iter=map_val.begin();iter!=map_val.end();++iter)
            {
                auto iter1(iter);
                iter1++;
                if(iter1!=map_val.end())
                {
                    const auto source(getSharedNode( iter->second));
                    const auto   sink(getSharedNode(iter1->second));
                    const VectorDim midPoint(0.5*(*source+*sink));
//                    patchEdges().emplace(std::make_pair(source->sID,sink->sID),PatchEdge(source,sink));
                    patchEdges().emplace(std::make_pair(source->sID,sink->sID),PatchEdge(source,sink,Polygon2D::windingNumber(midPoint,loop)));

                }
            }
            
        }
        
    }
    
    /**********************************************************************/
    template <class T>
    friend T& operator << (T& os, const LoopPathClipper& lpc)
    {
        for(const auto& edge : lpc.loopEdges())
        {
            os  << edge.second.source->transpose()<<" "<<edge.second.sink->transpose()<<" "<<edge.second.wn<<"\n";
        }
        for(const auto& edge : lpc.patchEdges())
        {
            os  << edge.second.source->transpose()<<" "<<edge.second.sink->transpose()<<" "<<edge.second.wn<<"\n";
        }
        return os;
    }
    
};



int main(int argc, char * argv[])
{
    
    std::vector<Eigen::Matrix<double,2,1>> polyPoints(readPolyVector("polyPoints.txt"));
    std::vector<Eigen::Matrix<double,2,1>> testPoints(readPolyVector("testPoints.txt"));
    std::vector<Eigen::Matrix<double,2,1>> patchPoints(readPolyVector("patchPoints.txt"));

    std::ofstream ofs("results.txt");
    for(const auto& testPt : testPoints)
    {
        //windingNumbers.push_back();
        ofs<<testPt.transpose()<<" "<<Polygon2D::windingNumber(testPt,polyPoints)<<"\n";
    }
    
    std::ofstream ofs1("edges.txt");

    
    LoopPathClipper lpc(polyPoints,patchPoints);
    ofs1<<lpc;
    
    return 0;
}
