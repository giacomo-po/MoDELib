
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <math.h>       /* fmod */

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
    void next(std::vector<const ClipperEdge*>&) const;
    const ClipperEdge* twin() const;
    bool isLoop() const;
    bool isPatch() const;
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

bool ClipperEdge::isLoop() const
{
    return source->loopNext==this && sink->loopPrev==this;
}

bool ClipperEdge::isPatch() const
{
    return source->patchNext==this && sink->patchPrev==this;
}

void ClipperEdge::next(std::vector<const ClipperEdge*>& validStarts) const
{
    const bool isLoopEdg(isLoop());
    const bool isPatchEdg(isPatch());
    assert(((isLoopEdg!=isPatchEdg) && (isLoopEdg||isPatchEdg)) && "edge must be either loop or patch");

    if(isLoopEdg)
    {// currently on loop
        assert(sink->loopNext && "loopNext must exist");
        if(sink->patchNext)
        {
            if(sink->patchNext->wn>sink->loopNext->wn)
            {
                if(sink->loopNext->wn!=0)
                {
                    validStarts.push_back(sink->loopNext);
                }
                if(!sink->patchNext->twin())
                {
                    validStarts.push_back(sink->patchNext);
                }
//                break;
//                return sink->patchNext;
            }
            else
            {
                if(sink->patchNext->wn!=0 && !sink->patchNext->twin())
                {
                    validStarts.push_back(sink->patchNext);
                }
                validStarts.push_back(sink->loopNext);
            }
        }
        else
        {
            validStarts.push_back(sink->loopNext);
        }
//        return sink->loopNext;
    }
    else
    {// currently on patch
        assert(sink->patchNext && "loopNext must exist");
        if(sink->loopNext)
        {
            if(sink->loopNext->wn>sink->patchNext->wn)
            {
                if(sink->patchNext->wn!=0)
                {
                    validStarts.push_back(sink->patchNext);
                }
                if(!sink->loopNext->twin())
                {
                    validStarts.push_back(sink->loopNext);
                }
//                return sink->loopNext;
            }
            else
            {
                if(sink->loopNext->wn!=0 && !sink->loopNext->twin())
                {
                    validStarts.push_back(sink->loopNext);
                }
                validStarts.push_back(sink->patchNext);
            }
        }
        else
        {
            validStarts.push_back(sink->patchNext);
        }
//        return sink->patchNext;
    }
}

const ClipperEdge* ClipperEdge::twin() const
{
    if(isLoop())
    {
        if(sink->patchPrev)
        {
            if(sink->patchPrev->source==source)
            {// a boundary edge
                return sink->patchPrev;
            }
        }
    }
    else if(isPatch())
    {
        if(sink->loopPrev)
        {
            if(sink->loopPrev->source==source)
            {// a boundary edge
                return sink->loopPrev;
            }
        }
    }
    else
    {
        assert(false);
    }
    return nullptr;
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
        std::map<float,std::shared_ptr<LoopPathClipperNode>> loopIntersections;  // intersection points along loop. Key with float precision to round off errors
        std::map<float,std::shared_ptr<LoopPathClipperNode>> patchIntersections; // intersection points along patch. Key with float precision to round off errors

        
        for(size_t i=0;i<loop.size();++i)
        {
            const size_t i1(i+1<loop.size()? i+1 :0);
            for(size_t j=0;j<patch.size();++j)
            {
                const size_t j1(j+1<patch.size()? j+1 :0);
                SegmentSegmentDistance<2> ssd(loop[i],loop[i1],patch[j],patch[j1]);
                for(const auto& inter : ssd.intersectionSegment())
                {
                    const double t(fmod(i+std::get<1>(inter), loop.size())); // intersection parameter in loopIntersections
                    const double u(fmod(j+std::get<2>(inter),patch.size())); // intersection parameter in patchIntersections
                    
//                    std::cout<<t<<" "<<u<<std::endl;
                    const auto loopIter(  loopIntersections.find(t));
                    const auto patchIter(patchIntersections.find(u));

                    if(    loopIter!= loopIntersections.end()
                       && patchIter!=patchIntersections.end())
                    {// both u and t exist, so intersection was already found
                        assert(loopIter->second==patchIter->second && "INTERSECTION MUST BE SAME NODE");
                        
                    }
                    else if(    loopIter!= loopIntersections.end()
                            && patchIter==patchIntersections.end())
                    {// u not found and t found. Intersection point exists on loop only
//                        assert(false && "Intersection point exists on patch only. FINISH HERE");
                        patchIntersections.emplace(u,loopIter->second);
                    }
                    else if(    loopIter== loopIntersections.end()
                            && patchIter!=patchIntersections.end())
                    {// u found and t not found. Intersection point exists on patch only
//                        assert(false && "Intersection point exists on loop only. FINISH HERE");
                        loopIntersections.emplace(t,patchIter->second);
                    }
                    else
                    {// neither u nor t exist. New intersection point
                        std::shared_ptr<LoopPathClipperNode> newNode( new LoopPathClipperNode(std::get<0>(inter)));
                        nodes().emplace(newNode->sID,newNode);
                        loopIntersections.emplace(t,newNode);
                        patchIntersections.emplace(u,newNode);
                    }
                }

            }
        }
        
        for(size_t i=0;i<loop.size();++i)
        {// insert loop vertices as loopIntersections points, if not already existing
            if(loopIntersections.find(i+0.0)==loopIntersections.end())
            {
                std::shared_ptr<LoopPathClipperNode> newNode( new LoopPathClipperNode(loop[i]));
                nodes().emplace(newNode->sID,newNode);
                loopIntersections.emplace(i+0.0,newNode);
            }
        }
        
        for(size_t j=0;j<patch.size();++j)
        {// insert patch vertices as patchintersection points, if not already existing
            if(patchIntersections.find(j+0.0)==patchIntersections.end())
            {
                std::shared_ptr<LoopPathClipperNode> newNode( new LoopPathClipperNode(patch[j]));
                nodes().emplace(newNode->sID,newNode);
                patchIntersections.emplace(j+0.0,newNode);
            }
        }

        
        for(auto iter=loopIntersections.begin();iter!=loopIntersections.end();++iter)
        {
            auto iter1(iter);
            iter1++;
            if(iter1==loopIntersections.end())
            {
                iter1=loopIntersections.begin();
            }
            const auto source(iter->second);
            const auto   sink(iter1->second);
            const VectorDim midPoint(0.5*(*source+*sink));
            loopEdges().emplace(std::piecewise_construct,
                                std::forward_as_tuple(source->sID,sink->sID),
                                std::forward_as_tuple(source,sink,Polygon2D::windingNumber(midPoint,patch))
                                );
        }

        for(auto iter=patchIntersections.begin();iter!=patchIntersections.end();++iter)
        {
            auto iter1(iter);
            iter1++;
            if(iter1==patchIntersections.end())
            {
                iter1=patchIntersections.begin();
            }
            const auto source(iter->second);
            const auto   sink(iter1->second);
            const VectorDim midPoint(0.5*(*source+*sink));
            patchEdges().emplace(std::piecewise_construct,
                                 std::forward_as_tuple(source->sID,sink->sID),
                                 std::forward_as_tuple(source,sink,Polygon2D::windingNumber(midPoint,loop))
                                 );
        }

    }

    void makePaths()
    {
        //        std::vector<std::vector<VectorDim>> clippedPolygons;
        std::vector<const ClipperEdge*> validStarts;
        
        for(const auto& edg : loopEdges())
        {
            if(edg.second.wn>0)
            {
                validStarts.push_back(&edg.second);
                break;
            }
//            std::vector<const ClipperEdge*> path;
//            clipPolygon(&edg.second,path,checkedEdges);
//            if(path.size())
//            {
//                //                assert(path.back()->sink==path.front()->source && "ERROR, OPEN PATH");
//                clippedPolygons().push_back(path);
//            }
        }
        
        while(validStarts.size())
        {
            std::vector<const ClipperEdge*> path;
            clipPolygon(path,validStarts);
            if(path.size())
            {
                clippedPolygons().push_back(path);
            }

        }
        
    }
    
    
    void clipPolygon(std::vector<const ClipperEdge*>& path,
                     std::vector<const ClipperEdge*>& validStarts) const
    {
        
//        if(checkedEdges.find(edg)==checkedEdges.end())
//        {// checkedEdges does not cointain edg
        
        
        const auto edg(validStarts.back());
        std::cout<<"edg= "<<edg->tag()<<std::endl;
        validStarts.pop_back();
            
            const auto twin(edg->twin());
            
            
//            if(edg->isLoop())
//            {// loop edges cannot be used twice
//                checkedEdges.insert(edg);
//            }
//            else
//            {
//                if(twin)
//                {
//                    checkedEdges.insert(twin);
//                }
//            }
        
            // need to determine if edg goes in path
            if(edg->wn>0)
            {// add edg to path
                path.push_back(edg);
            }
            else if(edg->wn==0)
            {// check if this is a boundary edge
                if(twin)
                {
                    path.push_back(twin);
                }
            }
            
            
            if(path.size())
            {
                if(path.back()->sink!=path.front()->source)
                {// edg is same as beginning of path. Stop
                    edg->next(validStarts);
                    clipPolygon(path,validStarts);
//                    clipPolygon(edg->next(),path,checkedEdges);
                }
            }
            else
            {
                edg->next(validStarts);
                clipPolygon(path,validStarts);
//                clipPolygon(edg->next(),path,checkedEdges);
            }
//        }
    }
    
//    void makePaths()
//    {
////        std::vector<std::vector<VectorDim>> clippedPolygons;
//        std::set<const ClipperEdge*> checkedEdges;
//
//        for(const auto& edg : loopEdges())
//        {
//            std::vector<const ClipperEdge*> path;
//            clipPolygon(&edg.second,path,checkedEdges);
//            if(path.size())
//            {
////                assert(path.back()->sink==path.front()->source && "ERROR, OPEN PATH");
//                clippedPolygons().push_back(path);
//            }
//        }
//    }
//
//
//    void clipPolygon(const ClipperEdge* edg,
//                   std::vector<const ClipperEdge*>& path,
//                   std::set<const ClipperEdge*>& checkedEdges) const
//    {
//
//        if(checkedEdges.find(edg)==checkedEdges.end())
//        {// checkedEdges does not cointain edg
//
//            std::cout<<"edg= "<<edg->tag()<<std::endl;
//
//
//            const auto twin(edg->twin());
//
//
//            if(edg->isLoop())
//            {// loop edges cannot be used twice
//                checkedEdges.insert(edg);
//            }
//            else
//            {
//                if(twin)
//                {
//                    checkedEdges.insert(twin);
//                }
//            }
//
//            // neet to determine if
//            if(edg->wn>0)
//            {// add edg to path
//                path.push_back(edg);
//            }
//            else if(edg->wn==0)
//            {// check if this is a boundary edge
//                if(twin)
//                {
//                    path.push_back(twin);
//                }
//            }
//
//
//            if(path.size())
//            {
//                if(path.back()->sink!=path.front()->source)
//                {// edg is same as beginning of path. Stop
//                    clipPolygon(edg->next(),path,checkedEdges);
//                }
//            }
//            else
//            {
//                clipPolygon(edg->next(),path,checkedEdges);
//            }
//        }
//    }
    
    
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
        os<<std::flush;
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
                        assert(lpc.clippedPolygons()[k].back()->sink==lpc.clippedPolygons()[k].front()->source && "ERROR, OPEN PATH");
    }

    
    return 0;
}
