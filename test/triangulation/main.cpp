#include <iostream>
#include <chrono>
#include <vector>
#include <deque>
#include <list>
#include <model/Geometry/earcut.hpp>
#include <model/IO/IDreader.h>
#include <model/IO/SequentialOutputFile.h>

#include <model/Geometry/PlanarPolygon.h>
#include <model/Geometry/SegmentSegmentDistance.h>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Constrained_Delaunay_triangulation_2.h>
//#include <CGAL/Delaunay_mesh_vertex_base_2.h>
//#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <vector>

#include <iterator>
// http://geomalgorithms.com/a09-_intersect-3.html
// https://stackoverflow.com/questions/27314724/how-divide-self-intersecting-polygon-into-simple-polygons?rq=1
// http://www.angusj.com/delphi/clipper.php
// https://stackoverflow.com/questions/24315355/using-cgal-to-partition-a-non-simple-polygon

using namespace model;

using Coord = double;
using Point = std::array<Coord, 2>;

//typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
//typedef K::Point_2                                              Point_2;
//typedef CGAL::Delaunay_mesh_vertex_base_2<K>                    Vb;
//typedef CGAL::Delaunay_mesh_face_base_2<K>                      Fb;
//typedef CGAL::Triangulation_data_structure_2<Vb, Fb>            Tds;
//typedef CGAL::Constrained_triangulation_2<K>      CDT;
////typedef CGAL::Constrained_triangulation_plus_2<CDT>                       CDTP;
//typedef CDT::Vertex_handle                                      Vertex_handle;
//typedef CDT::Face_iterator                                      Face_iterator;
//typedef CDT::Point          Point;


std::list<Point>::const_iterator beforeLast(const std::list<Point>& loopPoints)
{
std::list<Point>::const_iterator stopIter=loopPoints.end();
    stopIter--;
//    stopIter--;
    return stopIter;
}

Eigen::Map<const Eigen::Vector2d> i2v(std::list<Point>::iterator iter)
{

    return Eigen::Map<const Eigen::Vector2d>(iter->data());
}


std::list<Point>::iterator nextLoop(std::list<Point>& loopPoints,
                                    std::deque<std::deque<Point>>& loops,
                                       const std::list<Point>::iterator& startIter)
{

    const size_t loopID=loops.size();
    std::cout<<"Starting Loop "<<loopID<<std::endl;
    
    
//    std::deque<Point> tempLoop;
    loops.emplace_back();
    std::deque<Point>& currentLoop(*loops.rbegin());
    currentLoop.push_back(*startIter);
    std::cout<<"loop "<<loopID<<" adding ("<<currentLoop[0][0]<<","<<currentLoop[0][1]<<")"<<std::endl;

    
    std::list<Point>::iterator returnIter(startIter);
    std::advance(returnIter,1);
    
    bool closed=false;
    
    for(std::list<Point>::iterator iter0=startIter;iter0!=beforeLast(loopPoints);)
    {
        std::list<Point>::iterator iter1(iter0);
        std::advance(iter1,1);
        returnIter=iter1;
        
        std::list<Point>::iterator tempIter(iter1);
        std::advance(tempIter,1);

        
//        const Eigen::Map<const Eigen::Vector2d> i2v(iter0)(iter0->data());
//        const Eigen::Map<const Eigen::Vector2d> i2v(iter1)(iter1->data());

        //        bool closed=false;

        if((i2v(iter1)-Eigen::Map<const Eigen::Vector2d>(currentLoop[0].data())).norm()<FLT_EPSILON)
        {
            
            std::cout<<"loop "<<loopID<<" closed"<<std::endl;

//            currentLoop.push_back(*iter2);
//            returnIter=iter3;
//            closed=true;
            
//            std::cout<<"("<<i2v(iter0).transpose()<<")"
//            <<"("<<i2v(iter1).transpose()<<")"
//            <<"("<<i2v(iter2).transpose()<<")"
//            <<"("<<i2v(iter3).transpose()<<") "
//            <<closed<<std::endl;
            
            closed=true;
            break;
            
        }
        
        
        std::map<double,std::pair<std::list<Point>::iterator,Point>> intersectionMap;

        
        //bool closedIter=loopPoints.end();
        
        for(std::list<Point>::iterator iter2=tempIter;iter2!=loopPoints.end();++iter2)
        {

            
            std::list<Point>::iterator iter3(iter2);
            std::advance(iter3,1);
            
//            const Eigen::Map<const Eigen::Vector2d> i2v(iter2)(iter2->data());
//            const Eigen::Map<const Eigen::Vector2d> i2v(iter3)(iter3->data());
            
            
            

            
//            std::cout<<"intersecting "<<std::distance(loopPoints.begin(),iter0)
//            /*                 */<<"->"<<std::distance(loopPoints.begin(),iter1)
//            /*                 */<<" ## "<<std::distance(loopPoints.begin(),iter2)
//            /*                 */<<"->"<<std::distance(loopPoints.begin(),iter3)
//            /*                 */<<std::endl;
//            
//            std::cout<<i2v(iter0).transpose()
//            /*                 */<<"->"<<i2v(iter1).transpose()
//            /*                 */<<" ## "<<i2v(iter2).transpose()
//            /*                 */<<"->"<<i2v(iter3).transpose()
//            /*                 */<<std::endl;
            
            SegmentSegmentDistance<2> ssd(i2v(iter0),i2v(iter1),
                                          i2v(iter2),i2v(iter3));
            
            if( ssd.dMin<DBL_EPSILON
               && ssd.t>DBL_EPSILON
               && ssd.t<1.0-DBL_EPSILON
               && ssd.u>DBL_EPSILON
               && ssd.u<1.0-DBL_EPSILON
               )
            {// intersection found
                const Point X{{0.5*(ssd.x0(0)+ssd.x1(0)),0.5*(ssd.x0(1)+ssd.x1(1))}};
                intersectionMap.insert(std::make_pair(ssd.t,std::make_pair(iter3,X)));
            }

            

        }
        
//        if(closed)
//        {
//            break;
//            
////            return closedIter;
//            
//        }
//        else
//        {
            if(intersectionMap.size())
            {// segment i2v(iter0)->i2v(iter1) intersects another segment. The first intersection is with segment ending at intersectionMap.begin()->second.first
                
                //            std::cout<<"("<<intersectionMap.begin()->second.second[0]<<","<<intersectionMap.begin()->second.second[1]<<") with "<<std::distance(loopPoints.begin(),intersectionMap.begin()->second.first)<<std::endl;
                
                currentLoop.push_back(intersectionMap.begin()->second.second);
                std::cout<<"loop "<<loopID<<" adding ("<<intersectionMap.begin()->second.second[0]<<","<<intersectionMap.begin()->second.second[1]<<")"<<std::endl;

                
                loopPoints.insert(intersectionMap.begin()->second.first,intersectionMap.begin()->second.second);
                const std::list<Point>::iterator newStart=loopPoints.insert(iter1,intersectionMap.begin()->second.second);
                
//                std::cout<<"intersection found, loopPoints"<<std::endl;
//                for(const auto& point : loopPoints)
//                {
//                    std::cout<<point[0]<<","<<point[1]<<std::endl;
//                }
                
                
                iter0=nextLoop(loopPoints,loops,newStart);
                iter1=iter0;
                std::advance(iter1,1);
                returnIter=iter1;
                std::cout<<"loop "<<loopID<<" returning to ("<<iter0->operator[](0)<<","<<iter0->operator[](1)<<")"<<std::endl;

                if((i2v(iter0)-Eigen::Map<const Eigen::Vector2d>(currentLoop.rbegin()->data())).norm()<FLT_EPSILON)
                {
                    
                    std::cout<<"loop "<<loopID<<" closed 0"<<std::endl;
                    returnIter=iter0;

                    //            currentLoop.push_back(*iter2);
                    //            returnIter=iter3;
                    //            closed=true;
                    
                    //            std::cout<<"("<<i2v(iter0).transpose()<<")"
                    //            <<"("<<i2v(iter1).transpose()<<")"
                    //            <<"("<<i2v(iter2).transpose()<<")"
                    //            <<"("<<i2v(iter3).transpose()<<") "
                    //            <<closed<<std::endl;
                    
                    closed=true;
                    break;
                    
                }
                
                if((i2v(iter1)-Eigen::Map<const Eigen::Vector2d>(currentLoop[0].data())).norm()<FLT_EPSILON)
                {
                    
                    std::cout<<"loop "<<loopID<<" closed 1"<<std::endl;
                    returnIter=iter1;
                    
                    //            currentLoop.push_back(*iter2);
                    //            returnIter=iter3;
                    //            closed=true;
                    
                    //            std::cout<<"("<<i2v(iter0).transpose()<<")"
                    //            <<"("<<i2v(iter1).transpose()<<")"
                    //            <<"("<<i2v(iter2).transpose()<<")"
                    //            <<"("<<i2v(iter3).transpose()<<") "
                    //            <<closed<<std::endl;
                    
                    closed=true;
                    break;
                    
                }
                
            }
            else
            {
                
                currentLoop.emplace_back(*iter1);
                std::cout<<"loop "<<loopID<<" adding ("<<iter1->operator[](0)<<","<<iter1->operator[](1)<<")"<<std::endl;

                iter0++;
                
                
            }
//        }
        
        if(closed)
        {
            break;
        }

        
    }
    
    

    
    if(!closed)
    {
        if(loopID>0)
        {
            
            for(int k=1;k<currentLoop.size();++k)
            {
                std::cout<<"appending ("<<currentLoop[k][0]<<","<<currentLoop[k][1]<<")"<<std::endl;
                loops[loopID-1].push_back(currentLoop[k]);
            }
            
            std::cout<<"loop "<<loopID<<" clearing"<<std::endl;
//            returnIter=startIter;
            currentLoop.clear();
            
        }
        else
        {
            assert(0 && "LOOP 0 MUST CLOSE");
        }
        

    }
    // Current loop not closed
    std::cout<<"loop "<<loopID<<" returning "<<returnIter->operator[](0)<<","<<returnIter->operator[](1)<<std::endl;
    
    return returnIter;
    
}

std::deque<std::deque<Point>> splitSelfIntersectingPolygon(std::list<Point>& loopPoints)
{
    loopPoints.insert(loopPoints.end(),*loopPoints.begin()); // add a ghost point
    std::deque<std::deque<Point>> loops;
    nextLoop(loopPoints,loops,loopPoints.begin());
    return loops;
}


int main(int argc, char** argv)
{
    
    
    // The index type. Defaults to uint32_t, but you can also pass uint16_t if you know that your
    // data won't have more than 65536 vertices.
    using N = uint32_t;
    
    int fileID=0;
    model::IDreader<'I',1,2,double> reader;
    reader.read(fileID,true);
    
    // Create array
    //using Point = std::array<Coord, 2>;
    
    // Fill polygon structure with actual data. Any winding order works.
    // The first polyline defines the main polygon.
//    std::vector<Point> externalPoints;
    
    std::list<Point> pointList;
    //polygon[0].reserve(reader.size());
    for(const auto& pair : reader)
    {
//        polygon[0].push_back(pair.second);
        pointList.push_back(pair.second);
    }
    
    
    std::deque<std::deque<Point>> loops=splitSelfIntersectingPolygon(pointList);
    

    model::SequentialOutputFile<'P',true>::set_count(fileID);
    model::SequentialOutputFile<'P',true> pointsFile;
    for(const auto& point : pointList)
    {
        pointsFile<<point[0]<<" "<<point[1]<<std::endl;
    }
    
    model::SequentialOutputFile<'T',true>::set_count(fileID);
    model::SequentialOutputFile<'T',true> triFile;
    
    model::SequentialOutputFile<'X',true>::set_count(fileID);
    model::SequentialOutputFile<'X',true> xFile;
    
    double t0(clock());

    
    int loopID=0;
    for(const auto& loop : loops)
    {
        std::cout<<"Loop "<<loopID<<std::endl;
        for(const auto& point : loop)
        {
            std::cout<<" "<<point[0]<<" "<<point[1]<<std::endl;
            xFile<<loopID<<" "<<point[0]<<" "<<point[1]<<std::endl;
        }
        
        std::deque<std::deque<Point>> polygon;
//        polygon.resize(1); // first vector is external polyline, further vectors define holes
        polygon.push_back(loop); // first vector is external polyline, further vectors define holes

        
        const std::vector<N> indices = mapbox::earcut<N>(polygon);
        
        
        assert((indices.size()%3)==0);
        const size_t nTri=indices.size()/3;
        
        
        for(size_t k=0;k<nTri;++k)
        {
            triFile<<loopID<<" "<<indices[3*k+0]<<" "<<indices[3*k+1]<<" "<<indices[3*k+2]<<"\n";
        }

        
        
        
        loopID++;
    }
    
    std::cout<<" ("<<reader.size()<<" nodes) ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;


//    polygon.push_back({{100, 0}, {100, 100}, {0, 100}, {0, 0}});
//    // Following polylines define holes.
//    polygon.push_back({{75, 25}, {75, 75}, {25, 75}, {25, 25}});
    
    // Run tessellation
    // Returns array of indices that refer to the vertices of the input polygon.
    // e.g: the index 6 would refer to {25, 75} in this example.
    // Three subsequent indices form a triangle. Output triangles are clockwise.
//    double t0(clock());
//    const std::vector<N> indices = mapbox::earcut<N>(polygon);
//    std::cout<<" ("<<reader.size()<<" nodes) ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
//
//    
//    assert((indices.size()%3)==0);
//    const size_t nTri=indices.size()/3;
//
//    
//    for(size_t k=0;k<nTri;++k)
//    {
//        triFile<<indices[3*k+0]<<" "<<indices[3*k+1]<<" "<<indices[3*k+2]<<"\n";
//    }
    
    Eigen::Vector3d x(Eigen::Vector3d::Random());
    Eigen::Vector3d n(Eigen::Vector3d::Random().cross(x));
    model::PlanarPolygon pp(x,n);
    
    std::deque<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d>> vec3;
    for(int k=0;k<10;++k)
    {
        vec3.push_back(Eigen::Vector3d::Random());
    }

    pp.assignPoints(vec3);
    std::deque<std::array<size_t, 3>> triangles=pp.triangulate();
    
//    // Create a vector of the points
//    //
//    std::vector<Point_2> points2D ;
//    for(const auto& pair : reader)
//    {
//    points2D.push_back(Point_2(  pair.second[0],  pair.second[1]));
////        polygon[0].push_back(pair.second);
////        pointList.push_back(pair.second);
//    }
////    points2D.push_back(Point_2(  1,  1));
////    points2D.push_back(Point_2(  1, -1));
////    points2D.push_back(Point_2( -1,  1));
////    points2D.push_back(Point_2( -1, -1));
////    points2D.push_back(Point_2( 0, 0));
//    
//    size_t numTestPoints = points2D.size();
//    
//    // Create a constrained delaunay triangulation and add the points
//    //
//    CDT cdt;
////    std::vector<Vertex_handle> vhs;
//    for (unsigned int i=0; i<numTestPoints; ++i)
//    {
////        vhs.push_back(cdt.insert(points2D[i]));
//        int i1(i==numTestPoints-1? 0 : i+1);
//    cdt.insert_constraint( points2D[i], points2D[i1]);
//    }
//
//
//    
//    model::SequentialOutputFile<'C',true>::set_count(fileID);
//    model::SequentialOutputFile<'C',true> cgalFile;
//    for (Face_iterator fit = cdt.faces_begin()  ; fit != cdt.faces_end(); ++fit)
//    {
//        cgalFile<<fit->vertex(0)->point().x()<<" "<<fit->vertex(0)->point().y()<<"\n"
//        /*    */<<fit->vertex(1)->point().x()<<" "<<fit->vertex(1)->point().y()<<"\n"
//        /*    */<<fit->vertex(2)->point().x()<<" "<<fit->vertex(2)->point().y()<<"\n";
//    }
//
//    
//    int i=0;
//    for (Face_iterator fit = cdt.faces_begin()  ; fit != cdt.faces_end(); ++fit) {
//        printf("Face %d is (%f,%f) -- (%f,%f) -- (%f,%f) \n",i++,
//               fit->vertex(0)->point().x(),fit->vertex(0)->point().y(),
//               fit->vertex(1)->point().x(),fit->vertex(1)->point().y(),
//               fit->vertex(2)->point().x(),fit->vertex(2)->point().y() );
//        
//    }

    
    return 0;
}
