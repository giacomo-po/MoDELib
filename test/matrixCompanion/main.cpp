/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <map>
#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>
#include <model/Network/Readers/EdgeReader.h>
#include <model/Geometry/Splines/Coeff2Hermite.h>
#include <model/Math/MatrixCompanion.h>


using namespace model;

template<int dim,int Ncoeff>
Eigen::Matrix<double,dim,1> get_P(const Eigen::Matrix<double,dim,Ncoeff>& coeffs,const double& u)
{
    Eigen::Matrix<double,3,1> temp(Eigen::Matrix<double,3,1>::Zero());
    for(int k=0;k<Ncoeff;++k)
    {
        temp+=coeffs.col(k)*std::pow(u,k);
    }
    return temp;
}


/******************************************************************************/
template<int dim,int Ncoeff>
std::pair<double,std::pair<double,Eigen::Matrix<double,dim,1>> > closestPoint(const Eigen::Matrix<double,dim,Ncoeff>& coeffs1,
                                                                              const Eigen::Matrix<double,dim,1>& P0)
{
    
    // solve (P-P0)*dP/du=0
    
    // The polynomial coefficients of this spine segment
    
    Eigen::Matrix<double,dim,Ncoeff> coeffs(coeffs1);
    coeffs.col(0)-=P0;
    
    
    // The derivative of the polynomial coefficients
    Eigen::Matrix<double,dim,Ncoeff-1> dcoeffs(Eigen::Matrix<double,dim,Ncoeff-1>::Zero());
    for (int i=0;i<Ncoeff-1;++i)
    {
        dcoeffs.col(i)=(i+1)*coeffs.col(i+1);
    }
    
    Eigen::Matrix<double,1,2*Ncoeff-2> pcoeffs(Eigen::Matrix<double,1,2*Ncoeff-2>::Zero()); // degree of product = pOrder+(pOrder-1)=2*pOrder-1. nCoeffs of product = 2*pOrder-1+1= 2*pOrder = 2*Ncoeff-2
    
    // The polynomial coefficients of (P-P0)*dP/du, in reverse order
    for (int i=0;i<Ncoeff;++i)
    {
        for (int j=0;j<Ncoeff-1;++j)
        {
            pcoeffs(2*Ncoeff-3-i-j) += coeffs.col(i).dot(dcoeffs.col(j));
        }
    }
    
    // Compute roots using the eigenvalue method
    MatrixCompanion<2*Ncoeff-3> mc(pcoeffs);
    
    // sort roots according to distance to P0
    std::map<double,std::pair<double,Eigen::Matrix<double,dim,1>> > rootMap;
    
    //    for (int k=0;k<2*Ncoeff-3;++k)
    for (size_t k=0;k<mc.rootSize;++k)
    {
        std::cout<<"complex root: "<<mc.root(k);
        
        if (std::fabs(mc.root(k).imag())<FLT_EPSILON && mc.root(k).real()>-FLT_EPSILON && mc.root(k).real()<1.0+FLT_EPSILON )
//            if (std::fabs(mc.root(k).imag())<FLT_EPSILON && mc.root(k).real()>0.0 && mc.root(k).real()<1.0 )

        {
            
            std::cout<<" ... ok";
            
            
            Eigen::Matrix<double,dim,1> P(get_P(coeffs1,mc.root(k).real()));
            rootMap.insert(std::make_pair((P-P0).norm(), std::make_pair(mc.root(k).real(),P) ));
            
        }
        
        std::cout<<std::endl;
        
    }
    
    Eigen::Matrix<double,3,1> X0(get_P(coeffs1,0.0));
    Eigen::Matrix<double,3,1> X1(get_P(coeffs1,1.0));
    
    
    std::cout<<"X0="<<X0.transpose()<<std::endl;
    std::cout<<"X1="<<X1.transpose()<<std::endl;
    
    // check distance to limits of interval
    rootMap.insert(std::make_pair((X0-P0).norm(), std::make_pair(0.0,X0) ));
    rootMap.insert(std::make_pair((X1-P0).norm(), std::make_pair(1.0,X1) ));
    
    for(const auto& root : rootMap)
    {
        std::cout<<"root: "<<root.first<<", "<<root.second.first<<std::endl;
    }
    
    return *rootMap.begin();
    
}


void readH(const int& fileID,
           Eigen::Matrix<double,3,4>& H1,
           Eigen::Matrix<double,3,4>& H2)
{
    
    VertexReader<'V',9,double> vReader;
    EdgeReader<'E',11,double>  eReader;
    
    vReader.read(fileID,true);
    eReader.read(fileID,true);
    
    std::pair<int,int> link1(32,26);
    std::pair<int,int> link2(26,64);
    
    
    Eigen::Matrix<double,3,1> P0=vReader[link1.first].segment<3>(0);
    Eigen::Matrix<double,3,1> P1=vReader[link1.second].segment<3>(0);
    Eigen::Matrix<double,3,1> T0=vReader[link1.first].segment<3>(3)*eReader[link1](6);
    Eigen::Matrix<double,3,1> T1=vReader[link1.second].segment<3>(3)*eReader[link1](7);
    H1<<P0,T0,P1,T1;
    
    Eigen::Matrix<double,3,1> P2=vReader[link2.first].segment<3>(0);
    Eigen::Matrix<double,3,1> P3=vReader[link2.second].segment<3>(0);
    Eigen::Matrix<double,3,1> T2=vReader[link2.first].segment<3>(3)*eReader[link2](6);
    Eigen::Matrix<double,3,1> T3=vReader[link2.second].segment<3>(3)*eReader[link2](7);
    H2<<P2,T2,P3,T3;
    
}

void fillH(Eigen::Matrix<double,3,4>& H1,
           Eigen::Matrix<double,3,4>& H2)
{
    Eigen::Matrix<double,3,1> P0;
    P0<<1.485e+02, 3.620e+02, 2.000e+03;
    
    Eigen::Matrix<double,3,1> T0;
    T0<<0.000e+00, 5.173e-14, -1.035e-13;
    
    Eigen::Matrix<double,3,1> P1;
    P1<<1.110e+02, 3.246e+02, 2.000e+03;
    
    Eigen::Matrix<double,3,1> T1;
    T1<<3.014e+00, -5.063e+00, -8.077e+00;
    
    Eigen::Matrix<double,3,1> P2;
    P2<<7.425e+01, 2.878e+02, 2.000e+03;
    
    Eigen::Matrix<double,3,1> T2;
    T2<<6.594e+01, 6.594e+01, 1.456e-13;
    
    Eigen::Matrix<double,3,1> P3;
    P3<<1.485e+02, 3.620e+02, 2.000e+03;
    
    Eigen::Matrix<double,3,1> T3;
    T3<<0.000e+00,  7.281e-14, -1.456e-13;
    
    H1<<P0,T0,P1,T1;
    H2<<P2,T2,P3,T3;
    
}


int main(int argc, char** argv)
{
    
    Eigen::Matrix<double,3,4> H1;
    Eigen::Matrix<double,3,4> H2;
    //    readH(546,H1,H2);
    fillH(H1,H2);
    
    Eigen::Matrix<double,3,4> C1=Coeff2Hermite<3>::h2c<3>(H1);
    Eigen::Matrix<double,3,4> C2=Coeff2Hermite<3>::h2c<3>(H2);
    
    Eigen::Matrix<double,3,1> P1;
    P1<<1.110e+02, 3.246e+02, 2.000e+03;

    const auto temp=closestPoint(C2,P1);
    
    std::cout<<temp.second.first<<std::endl;
    
    return 0;
}
