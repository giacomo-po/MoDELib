
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
#include <Polygon2D.h>



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




int main(int argc, char * argv[])
{
    
    std::vector<Eigen::Matrix<double,2,1>> poly1(readPolyVector("poly1.txt"));
    std::vector<Eigen::Matrix<double,2,1>> poly2(readPolyVector("poly2.txt"));

    std::ofstream ofs("wn1.txt");
    
    const int I=1;
    
    for(size_t k=0;k<poly1.size();++k)
    {
        const size_t k1(k<poly1.size()-1? k+1 : 0);
        const auto P0(poly1[k]);
        const auto P1(poly1[k1]);

        for(int i=0;i<I;++i)
        {
            const Eigen::Matrix<double,2,1> q(P0*(1.0-double(i)/I)+P1*double(i)/I);
            ofs<<q.transpose()<<" "<<Polygon2D::windingNumber(q,poly2)<<"\n";

        }
    }
    


    
    return 0;
}
