#include <iostream>
#include <chrono>
#include <math.h>       /* round, floor, ceil, trunc */
#include <random>
#include <iomanip>
#include <limits>
#include <fstream>
#define _MODEL_NON_SINGULAR_DD_ 1
#include <model/DislocationDynamics/ElasticFields/StressStraight.h>
#include <model/DislocationDynamics/ElasticFields/DislocationStress.h>
#include <model/Geometry/SegmentSegmentDistance.h>

using namespace model;

std::random_device rd;
std::mt19937 generator;
std::uniform_real_distribution<> phiDis(0.0, 2.0*M_PI);
std::uniform_real_distribution<> thetaDis(0.0,   M_PI);
std::uniform_real_distribution<> preDis(1.0,   10.0);
std::uniform_int_distribution<> expDis(1,6);


template<typename ScalarType>
Eigen::Matrix<ScalarType,3,1> randomUnitVector()
{
    const ScalarType phi=phiDis(generator); // random angle of the dislocation line in the plane from screw orientation.
    const ScalarType theta=thetaDis(generator); // random angle of the dislocation line in the plane from screw orientation.
    return (Eigen::Matrix<ScalarType,3,1>()<<sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)).finished();
    
}

template<typename ScalarType>
ScalarType randomLength()
{
    const ScalarType p=preDis(generator); // random angle of the dislocation line in the plane from screw orientation.
    const int e=expDis(generator); // random angle of the dislocation line in the plane from screw orientation.
    return p*std::pow(10,e);
    
}

template<typename ScalarType>
void stressStraighDistance(const int& N,std::ofstream& myfile)
{
    
    //    const ScalarType LMin=
    
    auto t0 = std::chrono::system_clock::now();
    // Repeat N times
    for(int n=0;n<N;++n)
    {
        
        // Random direction of source segment
        const Eigen::Matrix<ScalarType,3,1> xi(randomUnitVector<ScalarType>());
        
        // Random Burgers vector of source segment
        Eigen::Matrix<ScalarType,3,1> b(randomUnitVector<ScalarType>());
        
        // End points of source segment
        const ScalarType L0(randomLength<ScalarType>());
        const Eigen::Matrix<ScalarType,3,1> P0(-0.5*xi*L0);
        const Eigen::Matrix<ScalarType,3,1> P1(+0.5*xi*L0);
        
        StressStraight<3,ScalarType> ss(P0,P1,b);
        
        // Random direction of center of field segment
        
        const ScalarType c(randomLength<ScalarType>());
        
        Eigen::Matrix<ScalarType,3,1> C(c*randomUnitVector<ScalarType>());
        // Random direction of field segment
        const Eigen::Matrix<ScalarType,3,1> xi1(randomUnitVector<ScalarType>());
        
        const ScalarType L1(randomLength<ScalarType>());
        const Eigen::Matrix<ScalarType,3,1> P2(C-0.5*xi1*L1);
        const Eigen::Matrix<ScalarType,3,1> P3(C+0.5*xi1*L1);
        
        SegmentSegmentDistance<3> ssd(P0,P1,P2,P3);
        
        const Eigen::Matrix<ScalarType,3,3> sP1=ss.stress(C); // stress at midpoint
        const Eigen::Matrix<ScalarType,3,3> sMin=ss.stress(ssd.x1); // stress  point of minimum distance
        const Eigen::Matrix<ScalarType,3,3> s2=ss.stress(P2); // stress at P2
        const Eigen::Matrix<ScalarType,3,3> s3=ss.stress(P3); // stress at P3
        //        const Eigen::Matrix<ScalarType,3,3> si=s2+ssd.u*(s3-s2); // stress at P3
        
        //        const ScalarType e((s-sC).norm()/s.norm());
        //        const ScalarType e1((s-si).norm()/s.norm());
        
        std::set<ScalarType> uSet; // parametric values along field segments where to compute stress
        uSet.insert(ssd.u); // insert point of minumum distance
        
        const int K=100;
        for(int k=0;k<K+1;++k)
        {
            uSet.insert(k*1.0/K);
        }
        
        double eMax1=0.0;
        double u1=-1.0;

        double eMax2=0.0;
        double u2=-1.0;

        double eMax3C=0.0;
        double eMax3=0.0;
        double u3=-1.0;

        
        for(const auto& u : uSet)
        {
            const Eigen::Matrix<ScalarType,3,1> P=P2+u*(P3-P2);   // point along field segment
            const Eigen::Matrix<ScalarType,3,3> sP=ss.stress(P);  // true stress at P
            
            // 1-point interpolation
            const double temp1=(sP-sP1).norm()/sP.norm();
            if(temp1>eMax1)
            {
                eMax1=temp1;
                u1=u;
            }
            
            // 2-point interpolation
            const Eigen::Matrix<ScalarType,3,3> sP2=s2+u*(s3-s2); // stress at P interpolated with two points
            const double temp2=(sP-sP2).norm()/sP.norm();
            if(temp2>eMax2)
            {
                eMax2=temp2;
                u2=u;
            }
            
            // 3-point interpolation with uMin
            const Eigen::Matrix<ScalarType,3,3> sP3=(u<=ssd.u)? s2+u/ssd.u*(sMin-s2) : sMin+(u-ssd.u)/(1.0-ssd.u)*(s3-sMin); // stress at P interpolated with three points using point of min distance
            if(ssd.u==0.0 || ssd.u==1.0)
            {
                eMax3= std::numeric_limits<double>::quiet_NaN();  // note, std::max removes the nan from sP3
                u3=std::numeric_limits<double>::quiet_NaN();
            }
            else
            {
                const double temp3=(sP-sP3).norm()/sP.norm();
                if(temp3>eMax3)
                {
                    eMax3=temp3;
                    u3=u;
                }
//                eMax3=std::max(eMax3,);  // note, std::max removes the nan from sP3
                
            }
            
            //            Eigen::Matrix<ScalarType,3,3> sP3(sP2);
            //            if(ssd.u>0.0 && ssd.u<1.0)
            //            {
            //                sP3=(u<=ssd.u)? s2+u/ssd.u*(sMin-s2) : sMin+(u-ssd.u)/(1.0-ssd.u)*(s3-sMin); // stress at P interpolated with three points using point of min distance
            //            }
            
            
            const Eigen::Matrix<ScalarType,3,3> sP3C=(u<=0.5)? s2+u/0.5*(sP1-s2) : sP1+(u-0.5)/(1.0-0.5)*(s3-sP1); // stress at P interpolated with three points using C
            

            

            


            eMax3C=std::max(eMax3C,(sP-sP3C).norm()/sP.norm());
            
        }
        
        myfile <<c<<" "<<ssd.dMin<<" "<<L0<<" "<<L1<<" "<<ssd.u<<" "<<u1<<" "<<u2<<" "<<u3<<" "<<eMax1<<" "<<eMax2<<" "<<eMax3<<" "<<eMax3C<<"\n";
        
    }
    
    std::cout<<std::chrono::duration<double>(std::chrono::system_clock::now()-t0).count()<<std::endl;
}

//template<typename ScalarType>
//void stressStraighDistance(const int& N,std::ofstream& myfile)
//{
//    auto t0 = std::chrono::system_clock::now();
//    // Repeat N times
//    for(int n=0;n<N;++n)
//    {
//
//        // Random direction of source segment
//        const Eigen::Matrix<ScalarType,3,1> xi(randomUnitVector<ScalarType>());
//
//        // Random Burgers vector of source segment
//        Eigen::Matrix<ScalarType,3,1> b(randomUnitVector<ScalarType>());
//
//        // End points of source segment
//        const ScalarType L0(1.0);
//        const Eigen::Matrix<ScalarType,3,1> P0(-0.5*xi*L0);
//        const Eigen::Matrix<ScalarType,3,1> P1(+0.5*xi*L0);
//
//        StressStraight<3,ScalarType> ss(P0,P1,b);
//
//        // Random direction of center of field segment
//        Eigen::Matrix<ScalarType,3,1> c(randomUnitVector<ScalarType>());
//
//        // Random direction of field segment
//        const Eigen::Matrix<ScalarType,3,1> xi1(randomUnitVector<ScalarType>());
//
//
//        for(int dExp=-2;dExp<3;++dExp)
//        {
//            // Distance of center of field segment from origin
//            for (int i=1;i<10;++i)
//            {
//                const ScalarType cNorm(i*std::pow(10.0,dExp));
//
//                for(int lExp=-3;lExp<3;++lExp)
//                {
//                    // Distance of center of field segment from origin
//                    for (int j=1;j<10;++j)
//                    {
//                        const ScalarType lNorm(j*std::pow(10.0,lExp)*cNorm);
//
//                        // End points of field segment
//                        const Eigen::Matrix<ScalarType,3,1> P2(cNorm*c-0.5*xi1*lNorm);
//                        const Eigen::Matrix<ScalarType,3,1> P3(cNorm*c+0.5*xi1*lNorm);
//
//                        SegmentSegmentDistance<3> ssd(P0,P1,P2,P3);
//
//
//                        const Eigen::Matrix<ScalarType,3,3> sC=ss.stress(cNorm*c); // stress at midpoint
//                        const Eigen::Matrix<ScalarType,3,3> s=ss.stress(ssd.x1); // stress at midpoint
//
//
//                        const Eigen::Matrix<ScalarType,3,3> s2=ss.stress(P2); // stress at P2
//                        const Eigen::Matrix<ScalarType,3,3> s3=ss.stress(P3); // stress at P3
//                        const Eigen::Matrix<ScalarType,3,3> si=s2+ssd.u*(s3-s2); // stress at P3
//
//                        const ScalarType e((s-sC).norm()/s.norm());
//                        const ScalarType e1((s-si).norm()/s.norm());
//
//                        const int K=100;
//                        double eMax2=e1;
//                        double eMaxC=e;
//                        double eMax3=0.0;
//                        double eMax3C=0.0;
//
//                        for(int k=0;k<K+1;++k)
//                        {
//                            const double u=k*1.0/K;
//                            const Eigen::Matrix<ScalarType,3,1> P=P2+u*(P3-P2);
//                            const Eigen::Matrix<ScalarType,3,3> sP=ss.stress(P); // true stress at P
//                            const Eigen::Matrix<ScalarType,3,3> sPi=s2+u*(s3-s2); // stress at P interpolated linearly
//
//                            eMax2=std::max(eMax2,(sP-sPi).norm()/sP.norm());
//                            eMaxC=std::max(eMaxC,(sP-sC).norm()/sP.norm());
//
//
//
//                            const Eigen::Matrix<ScalarType,3,3> sPi3=(u<=ssd.u)? s2+u/ssd.u*(s-s2) : s+(u-ssd.u)/(1.0-ssd.u)*(s3-s);
//                            const Eigen::Matrix<ScalarType,3,3> sPi3C=(u<=0.5)? s2+u/0.5*(sC-s2) : sC+(u-0.5)/(1.0-0.5)*(s3-sC);
//
//                            eMax3=std::max(eMax3,(sP-sPi3).norm()/sP.norm());
//                            eMax3C=std::max(eMax3C,(sP-sPi3C).norm()/sP.norm());
//                        }
//
//                        myfile << cNorm<<" "<<ssd.dMin<<" "<<lNorm<<" "<<e<<" "<<e1<<" "<<eMax2<<" "<<eMaxC<<" "<<eMax3<<" "<<eMax3C<<"\n";
//                        
//                        
//                    }
//                }
//                
//                
//                //                // End points of field segment
//                //                const Eigen::Matrix<ScalarType,3,1> P2(cNorm*c-0.5*xi1*L1);
//                //                const Eigen::Matrix<ScalarType,3,1> P3(cNorm*c+0.5*xi1*L1);
//                
//            }
//        }
//    }
//    
//    std::cout<<std::chrono::duration<double>(std::chrono::system_clock::now()-t0).count()<<std::endl;
//}


int main(int argc, char** argv)
{
    
    // Take meshID as a user input
    int N(1000);
    if (argc>1)
    {
        N=atoi(argv[1]);
    }
    
    std::ofstream myfile;
    myfile.open ("output.txt");
    stressStraighDistance<double>(N,myfile);
    myfile.close();
    
    
    
    
    
    return 0;
}
