

//#include <stdlib.h> // atoi
//#include <math.h>
//#include <chrono>
#define _MODEL_BENCH_SPLINEINTERSECTION_

#include <iostream>
#include <fstream>
#include <iomanip>      // std::setprecision
#include <model/Geometry/Splines/Intersection/PlanarSplineImplicitization.h>
#include <model/Geometry/Splines/Coeff2Hermite.h>
#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>

//#include <model/DislocationDynamics/DislocationNetwork.h>
//#include <model/Utilities/SequentialOutputFile.h>
//#include <model/ParticleInteraction/FieldPoint.h>
//#include <model/DislocationDynamics/StressStraight.h>

using namespace model;


Eigen::Matrix<double,3,3> rotMatrix(Eigen::Matrix<double,3,1> v1, Eigen::Matrix<double,3,1> n)
{
    v1.normalize();
    n.normalize();
    return (Eigen::Matrix<double,3,3>()<< v1, n.cross(v1), n).finished().transpose();
}


/******************************************************************************/
int main(int argc, char * argv[])
{
    
    typedef VertexReader<'V',9,double> VertexReaderType;
    
    
    
    constexpr int dim=3;
    constexpr int polyDegree=3;
    constexpr int polyCoeff=polyDegree+1;
    
    const double tol=FLT_EPSILON;
    
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    
    
    int fileID = 376;
    
    typedef VertexReader<'V',9,double> VertexReaderType;
    VertexReaderType  vReader;	// sID,Px,Py,Pz,Tx,Ty,Tz,snID
    vReader.read(fileID,true);
    
    typedef EdgeReader  <'E',11,double>	EdgeReaderType;
    EdgeReaderType    eReader;	// sourceID,sinkID,Bx,By,Bz,Nx,Ny,Nz,sourceF,sinkF,snID
    eReader.read(fileID,true);
    
    const int sourceID1 = 214;
    const int sinkID1 = 95;
    const std::pair<int,int> key1(sourceID1,sinkID1);
    const int sourceID2 = 595;
    const int sinkID2 = 213;
    const std::pair<int,int> key2(sourceID2,sinkID2);
    
    
    //    Eigen::Matrix<double,dim,polyCoeff> H1;
    //    H1<<
    //
    //    Eigen::Matrix<double,polyCoeff,dim>()<<
    //
    const Eigen::Matrix<double,dim,polyCoeff> H1=(Eigen::Matrix<double,polyCoeff,dim>()<<
                                                  vReader[sourceID1].segment(0,dim),
                                                  vReader[sourceID1].segment(dim,dim)*eReader[key1](2*dim),
                                                  vReader[sinkID1].segment(0,dim),
                                                  -vReader[sinkID1].segment(dim,dim)*eReader[key1](2*dim+1)
                                                  ).finished().transpose();
    
    const Eigen::Matrix<double,dim,polyCoeff> H2=(Eigen::Matrix<double,polyCoeff,dim>()<<
                                                  vReader[sourceID2].segment(0,dim),
                                                  vReader[sourceID2].segment(dim,dim)*eReader[key2](2*dim),
                                                  vReader[sinkID2].segment(0,dim),
                                                  -vReader[sinkID2].segment(dim,dim)*eReader[key2](2*dim+1)
                                                  ).finished().transpose();
    
    const VectorDim N1 = eReader[key1].segment(dim,dim);	  // end point of the spline
    const VectorDim P0 = H1.col(0);	  // end point of the spline
    const VectorDim T0 = H1.col(1);	  // tangent of the spline source
    const VectorDim P1 = H1.col(2);   // end point of the spline
    const VectorDim T1 = H1.col(3);   // tangent of the spline sink
    const VectorDim chord1 = P1-P0;	  // end point of the spline
    
    const VectorDim N2 = eReader[key2].segment(dim,dim);	  // end point of the spline
    const VectorDim P2 = H2.col(0);	  // end point of the spline
    const VectorDim T2 = H2.col(1);	  // tangent of the spline source
    const VectorDim P3 = H2.col(2);   // end point of the spline
    const VectorDim T3 = H2.col(3);   // tangent of the spline sink
    const VectorDim chord2 = P3-P2;	  // end point of the spline
    
    
    const Eigen::Matrix<double,3,3> R1=rotMatrix(chord1,N1);
    
    

    
    Eigen::Matrix<double,dim-1,polyCoeff> H1L;
    H1L.col(0)=(R1*(P0-P0)).template segment<dim-1>(0);
    H1L.col(1)=(R1*T0).template segment<dim-1>(0);
    H1L.col(2)=(R1*(P1-P0)).template segment<dim-1>(0);
    H1L.col(3)=(R1*T1).template segment<dim-1>(0);
    
    Eigen::Matrix<double,dim-1,polyCoeff> H2L;
    H2L.col(0)=(R1*(P2-P0)).template segment<dim-1>(0);
    H2L.col(1)=(R1*T2).template segment<dim-1>(0);
    H2L.col(2)=(R1*(P3-P0)).template segment<dim-1>(0);
    H2L.col(3)=(R1*T3).template segment<dim-1>(0);
    
    std::cout<<H1<<std::endl<<std::endl;
    std::cout<<H2<<std::endl<<std::endl;
    
    std::cout<<H1L<<std::endl<<std::endl;
    std::cout<<H2L<<std::endl<<std::endl;

    //
    PlanarSplineImplicitization<polyDegree> sli(Coeff2Hermite<polyDegree>::template h2c<dim-1>(H1L));
    
    std::set<std::pair<double,double> > intersectionParameters = sli.template intersectWith<polyDegree>(Coeff2Hermite<polyDegree>::template h2c<dim-1>(H2L));
    
    std::cout<<"Spline intersections are:"<<std::endl;
    for (std::set<std::pair<double,double> >::const_iterator sIter= intersectionParameters.begin();
         sIter!=intersectionParameters.begin();++sIter)
    {
        std::cout<<sIter->first<<" "<<sIter->second<<std::endl;
    }
    //    Coeff2Hermite<polyDegree>::template h2c<dim-1>(H2L);
    //
    //
    
    
    return 0;
}

