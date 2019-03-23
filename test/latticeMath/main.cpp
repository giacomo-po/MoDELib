#include <iostream>
#include <chrono>
//#include <model/DislocationDynamics/SnapToLattice.h>
//#include <model/DislocationDynamics/Materials/CrystalStructures.h>
#include <model/DislocationDynamics/Materials/FCCCrystal.h>
#include <model/DislocationDynamics/Materials/BCCCrystal.h>


#include <model/LatticeMath/LatticeMath.h>

using namespace model;




int main()
{

    
    typedef Eigen::Matrix<double,3,1> VectorDimD;
    
    const double a=1.0;
    const double b=0.5*sqrt(2.0)*a; // FCC
//    const double b=0.5*sqrt(3.0)*a; // BCC
    
    LatticeBase<3>::setLatticeBasis(FCC::getLatticeBasis<3>());
//    LatticeBase<3>::setLatticeBasis(BCC::getLatticeBasis<3>());

//    /*******************************************/
//    std::cout<<"slip normals are"<<std::endl;
//    std::vector<LatticePlaneBase> normals=FCC::template reciprocalPlaneNormals<3>();
//    for(unsigned int i=0;i<normals.size();++i)
//    {
//            std::cout<<normals[i].cartesian().transpose()<<std::endl;
//    }
    
//    ReciprocalLatticeVector<3> n(VectorDimD(3*a/b,0*a/b,-3*a/b));
//    std::cout<<n.transpose()<<std::endl;
    
//    LatticeVector<3>    v1(VectorDimD(1*a/b,0.5*a/b,0.5*a/b));
//    LatticeVector<3>    v2(VectorDimD(1*a/b,0.5*a/b,0.5*a/b));
//    LatticeDirection<3> d1(VectorDimD(0*a/b,-1*a/b,0*a/b));
//    LatticeDirection<3> d2(VectorDimD(0*a/b,0*a/b,-1*a/b));
//    
//
//    /*******************************************/
//    std::cout<<"slip directions are"<<std::endl;
//    std::vector<ReciprocalLatticeDirection<3>> normals=FCC::template reciprocalPlaneNormals<3>();
//    for(unsigned int i=0;i<normals.size();++i)
//    {
//        for(unsigned int j=i+1;j<normals.size();++j)
//        {
//            LatticeDirection<3> s(normals[i],normals[j]);
//            std::cout<<s.cartesian().transpose()<<std::endl;
//        }
//    }
//    
//    /*******************************************/
//
//    LatticePlane pl1(v1,normals[3]);
//
//    LatticePlane pl2(v2,normals[2]);
//    
//
//    LatticeLine l1(v1,d1);
//    LatticeLine l2(v2,d2);
//    
//    PlaneLineIntersection pli(pl1,l1);
//    std::cout<<pli<<std::endl;
//
//    PlanePlaneIntersection ppi(pl1,pl2);
//    std::cout<<ppi<<std::endl<<std::endl;
// 
//    std::cout<<"LineLineIntersection"<<std::endl;
//    LineLineIntersection lli(l1,l2);
//    std::cout<<lli<<std::endl<<std::endl;
    
    LatticeVector<3>    v1(VectorDimD(-0.5*sqrt(2.0),-sqrt(2.0),-0.5*sqrt(2.0)));
    std::cout<<v1<<std::endl;
    std::cout<<v1.cartesian()<<std::endl;
    
    return 0;
}