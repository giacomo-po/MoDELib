#include <iostream>
#include <chrono>
#include <vector>
#include <deque>
#include <model/Geometry/earcut.hpp>
#include <model/IO/IDreader.h>
#include <model/IO/SequentialOutputFile.h>

#include <model/Geometry/PlanarPolygon.h>


int main(int argc, char** argv)
{
    
    using Coord = double;
    
    // The index type. Defaults to uint32_t, but you can also pass uint16_t if you know that your
    // data won't have more than 65536 vertices.
    using N = uint32_t;
    
    int fileID=0;
    model::IDreader<'I',1,2,double> reader;
    reader.read(fileID,true);
    
    // Create array
    using Point = std::array<Coord, 2>;
    
    // Fill polygon structure with actual data. Any winding order works.
    // The first polyline defines the main polygon.
//    std::vector<Point> externalPoints;
    
    std::deque<std::deque<Point>> polygon;
    polygon.resize(1); // first vector is external polyline, further vectors define holes
    //polygon[0].reserve(reader.size());
    for(const auto& pair : reader)
    {
        polygon[0].push_back(pair.second);
    }
    
    

//    polygon.push_back({{100, 0}, {100, 100}, {0, 100}, {0, 0}});
//    // Following polylines define holes.
//    polygon.push_back({{75, 25}, {75, 75}, {25, 75}, {25, 25}});
    
    // Run tessellation
    // Returns array of indices that refer to the vertices of the input polygon.
    // e.g: the index 6 would refer to {25, 75} in this example.
    // Three subsequent indices form a triangle. Output triangles are clockwise.
    double t0(clock());
    const std::vector<N> indices = mapbox::earcut<N>(polygon);
    std::cout<<" ("<<reader.size()<<" nodes) ["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;

    
    assert((indices.size()%3)==0);
    const size_t nTri=indices.size()/3;

    model::SequentialOutputFile<'O',true>::set_count(fileID);
    model::SequentialOutputFile<'O',true> outFile;
    for(size_t k=0;k<nTri;++k)
    {
        outFile<<indices[3*k+0]<<" "<<indices[3*k+1]<<" "<<indices[3*k+2]<<"\n";
    }
    
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
    
    return 0;
}
