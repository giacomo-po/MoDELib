
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <math.h>       /* fmod */
#include <deque>

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


struct IntersectionNode
{
    typedef Eigen::Matrix<double,2,1> VectorDim;

    const float u;
    const VectorDim P;
    
    IntersectionNode(const float& u_in,
                     const VectorDim& P_in) :
    /* init */ u(u_in)
    /* init */,P(P_in)
    {
        
    }
};

///******************************************************************************/
struct LoopEdge //: public ClipperEdge
{
    
    typedef Eigen::Matrix<double,2,1> VectorDim;

    const IntersectionNode& source;
    const IntersectionNode& sink;
    const int wn;
    
    LoopEdge(const IntersectionNode& source_in,
             const IntersectionNode& sink_in,
             const int& wn_in) :
    /* init */ source(source_in)
    /* init */,sink(sink_in)
    /* init */,wn(wn_in)
    {
        
    }
    
};


/******************************************************************************/
struct LoopPathClipper : private std::map<float,IntersectionNode>
/*                    */,private std::deque<LoopEdge>
//: private std::map<std::pair<size_t,size_t>,LoopEdge>
///*                    */,private std::map<std::pair<size_t,size_t>,PatchEdge>
///*                    */,private std::map<size_t,const std::shared_ptr<LoopPathClipperNode>>
///*                    */,private std::vector<std::vector<const ClipperEdge*>>
{
    typedef Eigen::Matrix<double,2,1> VectorDim;
//    typedef std::map<size_t,const std::shared_ptr<LoopPathClipperNode>> NodeCointainerType;

    std::map<float,IntersectionNode>& nodes()
    {
        return *this;
    }

    const std::map<float,IntersectionNode>& nodes() const
    {
        return *this;
    }
//
    const std::deque<LoopEdge>& loopEdges() const
    {
        return *this;
    }

    std::deque<LoopEdge>& loopEdges()
    {
        return *this;
    }

    
    LoopPathClipper(const std::vector<VectorDim>& loop,
                    const std::vector<VectorDim>& patch)
    {

//        std::map<float,std::pair<float,VectorDim>> loopIntersections;  // intersection points along loop. Key with float precision to round off errors


//        std::map<float,std::shared_ptr<LoopPathClipperNode>> loopIntersections;  // intersection points along loop. Key with float precision to round off errors
//        std::map<float,std::shared_ptr<LoopPathClipperNode>> patchIntersections; // intersection points along patch. Key with float precision to round off errors

        
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
                    const auto loopIter(  nodes().find(t));
//                    const auto patchIter(patchIntersections.find(u));

                    if(    loopIter!= nodes().end())
                    {// both u and t exist, so intersection was already found
                        std::cout<<"t="<<t<<std::endl;
                        std::cout<<"u="<<u<<std::endl;
                        std::cout<<"loopIter->second.u="<<loopIter->second.u<<std::endl;
                        assert(std::fabs(loopIter->second.u-u)<FLT_EPSILON && "INTERSECTION MUST BE SAME NODE");
                    }
                    else
                    {
                        nodes().emplace(t,IntersectionNode(u,std::get<0>(inter)));
                    }
                }

            }
        }
        
        if(nodes().size())
        {// we have intersections between loop and patch
            for(size_t i=0;i<loop.size();++i)
            {// insert loop vertices as loopIntersections points, if not already existing
                if(nodes().find(i+0.0)==nodes().end())
                {
                    nodes().emplace(i+0.0,IntersectionNode(-1.0,loop[i])); // here -1.0 means not an intersection with patch
                }
            }
            
            for(auto iter=nodes().begin();iter!=nodes().end();++iter)
            {
                auto iter1(iter);
                iter1++;
                if(iter1==nodes().end())
                {
                    iter1=nodes().begin();
                }
                const auto& source(iter->second);
                const auto&   sink(iter1->second);
                const VectorDim midPoint(0.5*(source.P+sink.P));
                loopEdges().emplace_back(source,sink,Polygon2D::windingNumber(midPoint,patch));
            }
            
            makeIntersectingPaths();
            
        }
        else
        {// no intersection between loop and patch
            assert(false && "FINISH HERE CASE OF NO INTERSECTION");
        }
        


    }
//
    void makeIntersectingPaths()
    {
        //        std::vector<std::vector<VectorDim>> clippedPolygons;
//        std::vector<const ClipperEdge*> validStarts;
//
//        for(const auto& edg : loopEdges())
//        {
//            if(edg.second.wn>0)
//            {
//                validStarts.push_back(&edg.second);
//                break;
//            }
////            std::vector<const ClipperEdge*> path;
////            clipPolygon(&edg.second,path,checkedEdges);
////            if(path.size())
////            {
////                //                assert(path.back()->sink==path.front()->source && "ERROR, OPEN PATH");
////                clippedPolygons().push_back(path);
////            }
//        }

        
//        std::vector<std::pair<size_t,size_t>> inOutIDs;
        
        while(true)
        {// make sure that loopEdges.front() is outside
            if(loopEdges().front().wn!=0)
            {
                loopEdges().push_back(loopEdges().front());
                loopEdges().pop_front();
            }
            else
            {
                break;
            }
        }
        
        std::vector<std::vector<const LoopEdge*>> insideEdges;

        int insideID(-1);
        for(size_t k=0;k<loopEdges().size();++k)
        {
            const size_t k1(k==loopEdges().size()-1? 0 : k+1);
            
            const int& wnk(loopEdges()[k].wn);
            const int& wnk1(loopEdges()[k1].wn);
            
            if(wnk==0 && wnk1!=0)
            {// out->in transition
                //insideID=k1;
                assert(insideID==-1);
                insideEdges.push_back(std::vector<const LoopEdge*>());
            }
            
            if(wnk!=0)
            {
                insideEdges.back().push_back(&loopEdges()[k]);
            }
            
            if(wnk!=0 && wnk1==0 && insideID>=0)
            {// in->out transition
                insideID=-1;
            }
        }
        
        // HERE INSERT PATCH VERTICES with u in between start and end nodes
        
//        for(const auto& v : insideEdges)
//        {
//            const float& tStart(v.);
//            std::cout<<"path"<<std::endl;
//            for(const auto& e : v)
//            {
//                std::cout<<e->source.transpose()<<" "<<e->sink.transpose()<<std::endl;
//            }
//        }
        
        
        std::cout<<"insideEdges"<<std::endl;
        for(const auto& v : insideEdges)
        {
            std::cout<<"path"<<std::endl;
            for(const auto& e : v)
            {
                std::cout<<e->source.P.transpose()<<" "<<e->sink.P.transpose()<<std::endl;
            }
        }
        
        
        
//        while(validStarts.size())
//        {
//            std::vector<const ClipperEdge*> path;
//            clipPolygon(path,validStarts);
//            if(path.size())
//            {
//                clippedPolygons().push_back(path);
//            }
//
//        }

    }
//
//
//    void clipPolygon(std::vector<const ClipperEdge*>& path,
//                     std::vector<const ClipperEdge*>& validStarts) const
//    {
//
////        if(checkedEdges.find(edg)==checkedEdges.end())
////        {// checkedEdges does not cointain edg
//
//
//        const auto edg(validStarts.back());
//        std::cout<<"edg= "<<edg->tag()<<std::endl;
//        validStarts.pop_back();
//
//            const auto twin(edg->twin());
//
//
////            if(edg->isLoop())
////            {// loop edges cannot be used twice
////                checkedEdges.insert(edg);
////            }
////            else
////            {
////                if(twin)
////                {
////                    checkedEdges.insert(twin);
////                }
////            }
//
//            // need to determine if edg goes in path
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
//                    edg->next(validStarts);
//                    clipPolygon(path,validStarts);
////                    clipPolygon(edg->next(),path,checkedEdges);
//                }
//            }
//            else
//            {
//                edg->next(validStarts);
//                clipPolygon(path,validStarts);
////                clipPolygon(edg->next(),path,checkedEdges);
//            }
////        }
//    }
    
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
            os  << edge.source.P.transpose()<<" "<<edge.sink.P.transpose()<<" "<<edge.wn<<"\n";
        }
//        for(const auto& edge : lpc.patchEdges())
//        {
//            os  << edge.second.source->transpose()<<" "<<edge.second.sink->transpose()<<" "<<edge.second.wn<<"\n";
//        }
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
        ofs2<<node.first<<" "<<node.second.u<<" "<<node.second.P.transpose()<<std::endl;
    }
    



//    std::cout<<"NODE#="<<StaticID<LoopPathClipperNode>::get_count()<<std::endl;
    std::cout<<"lpc.nodes().size()="<<lpc.nodes().size()<<std::endl;

//    lpc.makePaths();

//    std::cout<<"CLIPPED POLYGONS: "<<lpc.clippedPolygons().size()<<std::endl;
//    for(size_t k=0;k<lpc.clippedPolygons().size();++k)
//    {
//        std::cout<<"Polygon "<<k<<std::endl;
//        for(const auto& edg : lpc.clippedPolygons()[k])
//        {
//            std::cout<<edg->tag()<<std::endl;
//        }
//                        assert(lpc.clippedPolygons()[k].back()->sink==lpc.clippedPolygons()[k].front()->source && "ERROR, OPEN PATH");
//    }

    
    return 0;
}
