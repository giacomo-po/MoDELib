
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <set>

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

/******************************************************************************/
struct ClipperEdge
{
    const std::shared_ptr<LoopPathClipperNode> source;
    const std::shared_ptr<LoopPathClipperNode> sink;
    const int wn;
    
    ClipperEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
              const std::shared_ptr<LoopPathClipperNode>& sink_in,
              const int& wn_in);
    std::string tag() const;
    ClipperEdge(const ClipperEdge& other) =delete;
    ClipperEdge(ClipperEdge&& other) =delete;
    ClipperEdge& operator=(const ClipperEdge& other)=delete;
    ClipperEdge& operator=(ClipperEdge&& other)=delete;
};

/******************************************************************************/
struct LoopEdge : public ClipperEdge
{
    LoopEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
             const std::shared_ptr<LoopPathClipperNode>& sink_in,
             const int& wn_in);
};

struct PatchEdge : public ClipperEdge
{
    PatchEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
              const std::shared_ptr<LoopPathClipperNode>& sink_in,
              const int& wn_in);
};

/******************************************************************************/
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
    
    LoopPathClipperNode(const LoopPathClipperNode& other) =delete;
    LoopPathClipperNode(LoopPathClipperNode&& other) =delete;
    LoopPathClipperNode& operator=(const LoopPathClipperNode& other)=delete;
    LoopPathClipperNode& operator=(LoopPathClipperNode&& other)=delete;
    
    /**********************************************************************/
    template <class T>
    friend T& operator << (T& os, const LoopPathClipperNode& node)
    {
        os<<node.sID<<" "<<node.transpose();
        return os;
    }
    
};

/******************************************************************************/
ClipperEdge::ClipperEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
                   const std::shared_ptr<LoopPathClipperNode>& sink_in,
                   const int& wn_in) :
/* init */ source(source_in)
/* init */,sink(sink_in)
/* init */,wn(wn_in)
{
}

std::string ClipperEdge::tag() const
{
    assert(source.get() && "BAD SOURCE");
        assert(sink.get() && "BAD SINK");
    return std::to_string(source->sID)+"->"+std::to_string(sink->sID);
}

/******************************************************************************/
LoopEdge::LoopEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
         const std::shared_ptr<LoopPathClipperNode>& sink_in,
         const int& wn_in) :
ClipperEdge(source_in,sink_in,wn_in)
{
    assert(source->loopNext==nullptr);
    source->loopNext=this;
    assert(sink->loopPrev==nullptr);
    sink->loopPrev=this;
}

/******************************************************************************/
PatchEdge::PatchEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
                     const std::shared_ptr<LoopPathClipperNode>& sink_in,
                     const int& wn_in) :
ClipperEdge(source_in,sink_in,wn_in)
{
    assert(source->patchNext==nullptr);
    source->patchNext=this;
    assert(sink->patchPrev==nullptr);
    sink->patchPrev=this;
}

/******************************************************************************/
struct LoopPathClipper : private std::map<std::pair<size_t,size_t>,LoopEdge>
/*                    */,private std::map<std::pair<size_t,size_t>,PatchEdge>
/*                    */,private std::map<size_t,const std::shared_ptr<LoopPathClipperNode>>
/*                    */,private std::vector<std::vector<const ClipperEdge*>>
{
    typedef Eigen::Matrix<double,2,1> VectorDim;
    typedef std::map<size_t,const std::shared_ptr<LoopPathClipperNode>> NodeCointainerType;

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

    const std::vector<std::vector<const ClipperEdge*>>& clippedPolygons() const
    {
        return *this;
    }

    
    std::vector<std::vector<const ClipperEdge*>>& clippedPolygons()
    {
        return *this;
    }
    
    
//    /**********************************************************************/
//    std::shared_ptr<LoopPathClipperNode> getSharedNode(const VectorDim& point)
//    {
//        const auto iter(nodes().find(point));
//        if(iter==nodes().end())
//        {// point does not exist
//            return nodes().emplace(point,std::shared_ptr<LoopPathClipperNode >(new LoopPathClipperNode(point))).first->second.lock();
//        }
//        else
//        {// point exists
//            if(iter->second.expired())
//            {// node deleted elsewhere
//                nodes().erase(iter);
//                return nodes().emplace(point,std::shared_ptr<LoopPathClipperNode>(new LoopPathClipperNode(point))).first->second.lock();
//            }
//            else
//            {
//                return iter->second.lock();
//            }
//        }
//    }
    
    LoopPathClipper(const std::vector<VectorDim>& loop,
                    const std::vector<VectorDim>& patch)
    {
        std::vector<std::map<double,std::shared_ptr<LoopPathClipperNode>>> loopIntersections;
        std::vector<std::map<double,std::shared_ptr<LoopPathClipperNode>>> patchIntersections;
        loopIntersections.resize(loop.size());
        patchIntersections.resize(patch.size());

        for(size_t i=0;i<loop.size();++i)
        {
            std::shared_ptr<LoopPathClipperNode> newNode( new LoopPathClipperNode(loop[i]));
            nodes().emplace(newNode->sID,newNode);
            const size_t prevI(i>0? i-1 : loop.size()-1);
            loopIntersections[i].emplace(0.0,newNode);
            loopIntersections[prevI].emplace(1.0,newNode);
        }
        
        for(size_t j=0;j<patch.size();++j)
        {
            std::shared_ptr<LoopPathClipperNode> newNode( new LoopPathClipperNode(patch[j]));
            nodes().emplace(newNode->sID,newNode);
            const size_t prevJ(j>0? j-1 : patch.size()-1);
            patchIntersections[j].emplace(0.0,newNode);
            patchIntersections[prevJ].emplace(1.0,newNode);
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
                    std::shared_ptr<LoopPathClipperNode> newNode( new LoopPathClipperNode(std::get<0>(inter)));
                    nodes().emplace(newNode->sID,newNode);
                    loopIntersections[i].emplace(std::get<1>(inter),newNode);
                    patchIntersections[j].emplace(std::get<2>(inter),newNode);
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
                    const auto source(iter->second);
                    const auto   sink(iter1->second);
                    const VectorDim midPoint(0.5*(*source+*sink));
//                    loopEdges().emplace(std::make_pair(source->sID,sink->sID),LoopEdge(source,sink,Polygon2D::windingNumber(midPoint,patch)));
                    loopEdges().emplace(std::piecewise_construct,
                                        std::forward_as_tuple(source->sID,sink->sID),
                                        std::forward_as_tuple(source,sink,Polygon2D::windingNumber(midPoint,patch))
                                        );

                
                    
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
                    const auto source(iter->second);
                    const auto   sink(iter1->second);
                    const VectorDim midPoint(0.5*(*source+*sink));
                    patchEdges().emplace(std::piecewise_construct,
                                        std::forward_as_tuple(source->sID,sink->sID),
                                        std::forward_as_tuple(source,sink,Polygon2D::windingNumber(midPoint,loop))
                                        );
                }
            }
            
        }
    }
    
    void makePaths()
    {
//        std::vector<std::vector<VectorDim>> clippedPolygons;
        std::set<const ClipperEdge*> checkedEdges;
        
        for(const auto& edg : loopEdges())
        {
            std::vector<const ClipperEdge*> path;
            walkPath(&edg.second,path,checkedEdges);
            if(path.size())
            {
                clippedPolygons().push_back(path);
            }
        }
    }
    
    
    int walkPath(const ClipperEdge* edg,
                   std::vector<const ClipperEdge*>& path,
                   std::set<const ClipperEdge*>& checkedEdges) const
    {
        
        if(checkedEdges.find(edg)==checkedEdges.end())
        {// checkedEdges does not cointain edg

            const bool isLoopEdg(edg->sink->loopPrev==edg);
            const bool isPatchEdg(edg->sink->patchPrev==edg);

            if(isLoopEdg)
            {// loop edges cannot be used twice
                checkedEdges.insert(edg);
            }

            if(edg->wn>0)
            {// add edg to path
                path.push_back(edg);
            }
            
            if(path.size())
            {
                if(edg->sink==path.front()->source)
                {// edg is same as beginning of path. Stop
                    return 0;
                }
            }
            
            if(isLoopEdg)
            {// currently on loop
                assert(edg->sink->loopNext && "loopNext must exist");
                if(edg->sink->patchNext)
                {
                    if(edg->sink->patchNext->wn>edg->sink->loopNext->wn)
                    {
                        return 1+walkPath(edg->sink->patchNext,path,checkedEdges);
                    }
                }
                return 1+walkPath(edg->sink->loopNext,path,checkedEdges);
            }
            else if(isPatchEdg)
            {// currently on patch
                assert(edg->sink->patchNext && "loopNext must exist");
                if(edg->sink->loopNext)
                {
                    if(edg->sink->loopNext->wn>edg->sink->patchNext->wn)
                    {
                        return 1+walkPath(edg->sink->loopNext,path,checkedEdges);
                    }
                }
                return 1+walkPath(edg->sink->patchNext,path,checkedEdges);
            }
            else
            {// impossible
                std::cout<<"edg->sink="<<edg->sink->sID<<std::endl;
                if(edg->sink->loopPrev)
                {
                    std::cout<<"loopPrev="<<edg->sink->loopPrev->tag()<<std::endl;
                }
                if(edg->sink->patchPrev)
                {
                    std::cout<<"patchPrev="<<edg->sink->patchPrev->tag()<<std::endl;
                }
                if(edg->sink->loopNext)
                {
                    std::cout<<"loopNext="<<edg->sink->loopNext->tag()<<std::endl;
                }
                if(edg->sink->patchNext)
                {
                    std::cout<<"patchNext="<<edg->sink->patchNext->tag()<<std::endl;
                }

                assert(0);
            }
        }
        return 0;
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
    

    
    LoopPathClipper lpc(polyPoints,patchPoints);

        std::ofstream ofs1("edges.txt");
    ofs1<<lpc;

    std::ofstream ofs2("nodes.txt");
    for(const auto& node : lpc.nodes())
    {
        ofs2<<*node.second<<std::endl;
    }
    



    std::cout<<"NODE#="<<StaticID<LoopPathClipperNode>::get_count()<<std::endl;
    std::cout<<"lpc.nodes().size()="<<lpc.nodes().size()<<std::endl;

    lpc.makePaths();

    std::cout<<"CLIPPED POLYGONS: "<<lpc.clippedPolygons().size()<<std::endl;
    for(size_t k=0;k<lpc.clippedPolygons().size();++k)
    {
        std::cout<<"Polygon "<<k<<std::endl;
        for(const auto& edg : lpc.clippedPolygons()[k])
        {
            std::cout<<edg->tag()<<std::endl;
        }
        
    }

    
    return 0;
}
