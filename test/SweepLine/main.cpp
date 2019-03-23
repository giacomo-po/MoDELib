#include <iostream>
#include <model/Geometry/SweepPlane.h>
#include <model/Geometry/SegmentSegmentDistance.h>
#include <Eigen/Dense>
#include <deque>
#include <chrono>

template<int dim>
struct SimpleLine
{
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    VectorDim P0;
    VectorDim P1;
    SimpleLine(const VectorDim& p0,const VectorDim& p1) :
    P0(p0),P1(p1)
    {
    
    }

};

template<int dim>
std::pair<size_t,double> timeN2method(const std::deque<SimpleLine<dim>>& sld,const double& collisionTol)
{
    const auto t0= std::chrono::system_clock::now();
    //std::cout<<"N^2 method... "<<std::flush;
    
    size_t nIntersections=0;
    for(size_t i=0;i<sld.size();++i)
    {
        for(size_t j=i+1;j<sld.size();++j)
        {
            model::SegmentSegmentDistance<dim> ssd(sld[i].P0,sld[i].P1,
                                   sld[j].P0,sld[j].P1);
            if(ssd.dMin<collisionTol)
            {
                nIntersections++;
            }
        }
    
    }
    
    //std::cout<<nIntersections<<" intersections"<<std::flush;
    //std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;

    //double x=;
    
    return std::make_pair(nIntersections,(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count());
}

template<int dim>
std::pair<size_t,double> timeSweepPlane(const std::deque<SimpleLine<dim>>& sld,const double& collisionTol)
{
    const auto t0= std::chrono::system_clock::now();
    //std::cout<<"SweepPlane method... "<<std::flush;
    
    model::SweepPlane<SimpleLine<dim>,dim> sw;
    for(const auto& line : sld)
    {
        sw.addSegment(line.P0(0),line.P1(0),line);
    }
    sw.computeIntersectionPairs();
    
    size_t nIntersections=0;
    for(const auto& pair : sw.potentialIntersectionPairs())
    {
        model::SegmentSegmentDistance<dim> ssd(pair.first->P0,pair.first->P1,
                                               pair.second->P0,pair.second->P1);
        if(ssd.dMin<collisionTol)
        {
            nIntersections++;
        }
    }
    
    //std::cout<<nIntersections<<" intersections"<<std::flush;
    //std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;

    return std::make_pair(nIntersections,(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count());
}


int main()
{

    const int dim=2;
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    std::deque<SimpleLine<dim>> sld;
    
    const double collisionTol=1e-2;
    
    const double Lmax=1.0;
    const double boxSize=100*Lmax;
    
    for(int n=0;n<6;++n)
    {
        const int N=std::pow(10,n);

        for(int k=0;k<N;++k)
        {
            const VectorDim P0(VectorDim::Random()*boxSize);
            const VectorDim P1(P0+VectorDim::Random()*Lmax);
            sld.emplace_back(P0,P1);
        }
        
        std::cout<<N<<std::flush;
        const std::pair<size_t,double> isp=timeSweepPlane(sld,collisionTol);
        std::cout<<" "<<isp.first<<" "<<isp.second<<std::flush;
        
        const std::pair<size_t,double> in2=timeN2method(sld,collisionTol);
        //assert(isp.first==in2.first);
        std::cout<<" "<<in2.first<<" "<<in2.second<<std::endl;

    }

    

    
    return 0;
}
