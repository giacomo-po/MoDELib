/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanePlaneIntersection_h_
#define model_PlanePlaneIntersection_h_

#include <set>
#include <assert.h>

#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/LatticePlane.h>
#include <model/LatticeMath/LatticeDirection.h>
//#include <model/LatticeMath/LatticeLine.h>
#include <model/Math/RoundEigen.h>

namespace model
{
    
    
    struct PlanePlaneIntersection
    {
        typedef Eigen::Matrix<long int,3,1> VectorDimI;
        typedef Eigen::Matrix<double,3,1> VectorDimD;
        
        enum IntersectionType
        {
            incident=0,
            parallel=1,
            coincident=2
        };
        
        const LatticeDirection<3> d;
        const IntersectionType intersectionType;
        const LatticeVector<3> P;
        
        
        
        
        /**********************************************************************/
        LatticeVector<3> midPoint(const LatticeVector<3>& p1, const ReciprocalLatticeDirection<3>& n1,
                                  const LatticeVector<3>& p2, const ReciprocalLatticeDirection<3>& n2,
                                  const LatticeDirection<3>& dir)
        {
            std::set<int> intSet;
            intSet.insert(0);
            intSet.insert(1);
            intSet.insert(2);
            
            int d=-1;
            for(int dd=0;dd<3;++dd)
            {
                if(dir(dd)!=0)
                {
                    d=dd;
                }
            }
            
            assert(d>=0 && "all components of dir are 0.");
            intSet.erase(d);
            
            const int a=*intSet.begin();
            const int b=*intSet.rbegin();
            
// before Eigen 3.2.5
//            Eigen::Matrix<long int,2,2> N;
//            N << n1(a), n1(b),
//            /**/ n2(a), n2(b);
//            
//            
//            
//            Eigen::Matrix<long int,2,1> PN;
//            PN<<p1.dot(n1),p2.dot(n2);
//            
//            Eigen::Matrix<long int,2,1> p=N.lu().solve(PN);
            
            Eigen::Matrix<double,2,2> N;
            N << n1(a), n1(b),
            /**/ n2(a), n2(b);
            
            
            
            Eigen::Matrix<double,2,1> PN;
            PN<<p1.dot(n1),p2.dot(n2);
            
            const Eigen::Matrix<double,2,1> tempD(N.inverse()*PN);
            
            const Eigen::Matrix<long int,2,1> p=RoundEigen<long int,2>::round(tempD);
            
            
            
            Eigen::Matrix<long int,3,1> temp;
            temp(a)=p(0);
            temp(b)=p(1);
            temp(d)=0;
            
//            std::cout<<"PlanePlaneIntersection: I'm here"<<std::endl;
            
            assert((temp-p1).dot(n1)==0);
            assert((temp-p2).dot(n2)==0);

            assert(dir.dot(n1)==0);
            assert(dir.dot(n2)==0);
            
            return temp;
        }
        
//        /**********************************************************************/
//        LatticeVector<3> midPoint(const LatticePlane& plane1,
//                                  const LatticePlane& plane2,
//                                  const LatticeDirection<3>& dir)
//        {
//            
//            const LatticeVector<3>& p1=plane1.P;
//            const ReciprocalLatticeDirection<3>& n1=plane1.n;
//            const LatticeVector<3>& p2=plane2.P;
//            const ReciprocalLatticeDirection<3>& n2=plane2.n;
//            
//            // Find the Cartesian point P which minimizes (P-p1)^2+(P-p2)^2 under
//            // the constraints (P-P1)*n1=0 and (P-P2)*n2=0
//            Eigen::Matrix<double,5,5> M(Eigen::Matrix<double,5,5>::Zero());
//            M.block<3,3>(0,0)=2.0*Eigen::Matrix<double,3,3>::Identity();
//            M.block<3,1>(0,3)=n1.cartesian();
//            M.block<3,1>(0,4)=n2.cartesian();
//            M.block<1,3>(3,0)=n1.cartesian().transpose();
//            M.block<1,3>(4,0)=n2.cartesian().transpose();
//            
//            Eigen::Matrix<double,5,1> b(Eigen::Matrix<double,5,1>::Zero());
//            b.segment<3>(0)=p1.cartesian()+p2.cartesian();
//            b(3)=p1.cartesian().dot(n1.cartesian());
//            b(4)=p2.cartesian().dot(n2.cartesian());
//            
//            const Eigen::Matrix<double,3,1> C=M.lu().solve(b).segment<3>(0);
//            
//            assert(fabs((C-p1.cartesian()).dot(n1.cartesian()))<FLT_EPSILON);
//            assert(fabs((C-p2.cartesian()).dot(n2.cartesian()))<FLT_EPSILON);
//            
//            // Starting from C, move along dir and snap o lattice
//            LatticeVector<3> temp(p1);
//            const Eigen::Matrix<double,3,1> dc=dir.cartesian();
//            const int kMax=100;
//            for (int k=0;k<kMax;++k)
//            {
//                temp=LatticeVector<3>(plane1.snapToLattice(C+k*1.0/kMax*dc));
//                if(plane2.contains(temp))
//                {
//                    break;
//                }
//            }
//            
////            std::cout<<"P1="<<p1.transpose()<<std::endl;
////            std::cout<<"n1="<<n1.transpose()<<" ("<<n1.cartesian().transpose()<<")"<<std::endl;
////            std::cout<<"P2="<<p2.transpose()<<std::endl;
////            std::cout<<"n1="<<n2.transpose()<<" ("<<n2.cartesian().transpose()<<")"<<std::endl;
//            
//            assert((temp-p1).dot(n1)==0);
//            assert((temp-p2).dot(n2)==0);
//            
//            
//            return temp;
//            
//        }
        
        /**********************************************************************/
        PlanePlaneIntersection(const LatticePlane& plane1, const LatticePlane& plane2) :
        /* init */ d(plane1.n,plane2.n),
        /* init */ intersectionType(d.squaredNorm()? incident : ((plane1.P-plane2.P).dot(plane1.n)? parallel : coincident)),
        /* init */ P(intersectionType==incident? midPoint(plane1.P,plane1.n,plane2.P,plane2.n,d) : (intersectionType==coincident? plane1.P : VectorDimI::Zero()))
//        /* init */ P(intersectionType==incident? midPoint(plane1,plane2,d) : (intersectionType==coincident? plane1.P : VectorDimI::Zero()))
        {
//         if(intersectionType==incident || intersectionType==coincident)
//         {
//             assert((P-plane1.P).dot(plane1.n)==0);
//             assert(d.dot(plane1.n)==0);
//
//             assert((P-plane2.P).dot(plane2.n)==0);
//             assert(d.dot(plane2.n)==0);
//         }
        }
        
        /**********************************************************************/
        friend std::ostream& operator << (std::ostream& os, const PlanePlaneIntersection& ppi)
        {
            os << ppi.intersectionType<<"\n"
            /**/<<ppi.d.cartesian().transpose()<<"\n"
            /**/<<ppi.P.cartesian().transpose();
            return os;
        }
        
    };
    
} // end namespace
#endif
