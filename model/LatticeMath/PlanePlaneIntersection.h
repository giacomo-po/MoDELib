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
            
            Eigen::Matrix<long int,2,2> N;
            N << n1(a), n1(b),
            /**/ n2(a), n2(b);
            
            Eigen::Matrix<long int,2,1> PN;
            PN<<p1.dot(n1),p2.dot(n2);
            
            
            Eigen::Matrix<long int,2,1> p=N.lu().solve(PN);
//            Eigen::Matrix<long int,2,1> p=N.inverse()*PN;
            
            
            
            Eigen::Matrix<long int,3,1> temp;
            temp(a)=p(0);
            temp(b)=p(1);
            temp(d)=0;
            
            std::cout<<"PlanePlaneIntersection: I'm here"<<std::endl;
            
            assert((temp-p1).dot(n1)==0);
            assert((temp-p2).dot(n2)==0);

            assert(dir.dot(n1)==0);
            assert(dir.dot(n2)==0);
            
            return temp;
        }
        
        /**********************************************************************/
        PlanePlaneIntersection(const LatticePlane& plane1, const LatticePlane& plane2) :
        /* init */ d(plane1.n,plane2.n),
        /* init */ intersectionType(d.squaredNorm()? incident : ((plane1.P-plane2.P).dot(plane1.n)? parallel : coincident)),
        /* init */ P(intersectionType==incident? midPoint(plane1.P,plane1.n,plane2.P,plane2.n,d) : (intersectionType==coincident? plane1.P : VectorDimI::Zero()))
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
