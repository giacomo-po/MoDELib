/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// The following line is needed for shape functions or order 5 or more

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <model/LatticeMath/LatticeMath.h>

using namespace model;

Eigen::Matrix<long int,3,1> randomOrigin()
{
    Eigen::Matrix<long int,3,1> temp;
    temp(0)=rand() % 100;
    temp(1)=rand() % 100;
    temp(2)=rand() % 100;
    return temp;
}

int main()
{
    srand (time(NULL)); // init rand
    
    static constexpr int dim=3;
    static constexpr PeriodicElement<13,Isotropic> Al=PeriodicElement<13,Isotropic>();
    typedef Eigen::Matrix<long int,dim,1> VectorDimI;
    typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
    typedef LatticeVector<dim> LatticeVectorType;

    const MatrixDimD covBasis=PeriodicElement<Al.Z,Isotropic>::CrystalStructure::template getLatticeBasis<dim>();
    const MatrixDimD contraBasis=covBasis.inverse().transpose();
    
    const LatticeVectorType a1(VectorDimI(1,0,0),covBasis,contraBasis);
    const LatticeVectorType a2(VectorDimI(0,1,0),covBasis,contraBasis);
    const LatticeVectorType a3(VectorDimI(0,0,1),covBasis,contraBasis);

    // First plane
    const LatticeVectorType P0(randomOrigin(),covBasis,contraBasis);
    LatticePlaneBase pb0(a1-a3,a2-a3);     // is (1,1,1) in cartesian
    LatticePlane pl0(P0,pb0);
    std::cout<<std::endl<<"plane0 origin="<<pl0.P.cartesian().transpose()<<std::endl;
    std::cout<<"plane0 normal="<<pl0.n.cartesian().transpose()<<std::endl;
    
    // Second plane
    const LatticeVectorType P1(randomOrigin(),covBasis,contraBasis);
    LatticePlaneBase pb1(a1-a2+a3,a1+2*a2-a3); // is (3,0,-1) in cartesian
    LatticePlane pl1(P1,pb1);
    std::cout<<std::endl<<"plane1 origin="<<pl1.P.cartesian().transpose()<<std::endl;
    std::cout<<"plane1 normal="<<pl1.n.cartesian().transpose()<<std::endl;
    
    // Plane-Plane intersection
    PlanePlaneIntersection ppi(pl0,pl1);
    std::cout<<std::endl<<"ppi origin="<<ppi.P.cartesian().transpose()<<std::endl;
    std::cout<<"ppi direction="<<ppi.d.cartesian().transpose()<<std::endl;
    
    return 0;
}




