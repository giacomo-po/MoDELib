/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePlaneBase_h_
#define model_LatticePlaneBase_h_

#include <map>
#include <LatticeVector.h>
#include <LatticeDirection.h>
#include <ReciprocalLatticeDirection.h>

namespace model
{
    struct LatticePlaneBase : public ReciprocalLatticeDirection<3>
    {
        typedef Eigen::Matrix<double,3,1> VectorDimD;
        typedef Eigen::Matrix<long int,3,1> VectorDimI;
        typedef LatticeVector<3>    LatticeVectorType;
        typedef LatticeDirection<3>    LatticeDirectionType;
        typedef ReciprocalLatticeDirection<3> ReciprocalLatticeDirectionType;
        
        const std::pair<LatticeDirectionType,LatticeDirectionType> primitiveVectors;
        
//        /**********************************************************************/
//        static std::pair<LatticeDirectionType,LatticeDirectionType> findPrimitiveVectors(const ReciprocalLatticeDirectionType& rd)
//        {
//            std::cout<<"LatticePlaneBase::findPrimitiveVectors WARNING, THIS FUNCTION IS TEMPORARY, IT MAY HANG for large r"<<std::endl;
//            std::cout<<"r ="<<rd.transpose()<<std::endl;
//            std::cout<<"r.cartesian()"<<rd.cartesian().transpose()<<std::endl;
//            
//            std::multimap<double,LatticeDirectionType> normMap;
//            
//            const long int N=rd.array().abs().maxCoeff()*2;
//            for(long int i=-N;i<=N;++i)
//            {
//                for(long int j=-N;j<=N;++j)
//                {
//                    for(long int k=-N;k<=N;++k) // THIS LOOP MAY BE REMOVED BY CHOOSING k SUCH THAT P.dot(rd)==0
//                    {
////                        const LatticeDirectionType P(LatticeVectorType((VectorDimI()<<i,j,k).finished(),rd.covBasis,rd.contraBasis));
//                        const LatticeDirectionType P(LatticeVectorType((VectorDimI()<<i,j,k).finished(),rd.lattice));
//
//                        const size_t latticeNorm2=P.squaredNorm();
//                        if(P.dot(rd)==0 && latticeNorm2>0)
//                        {
//                            normMap.emplace(P.cartesian().squaredNorm(),P); // can we use the latticeNorm2 instead?
//                        }
//                    }
//                }
//            }
//            
//            for (const auto& iter : normMap)
//            {
//                const ReciprocalLatticeDirectionType temp(normMap.begin()->second,iter.second);
//                if(temp.squaredNorm()!=0)
//                {
//                    if((temp-rd).squaredNorm()==0)
//                    {
//                        return std::make_pair(normMap.begin()->second,iter.second);
//                    }
//                    else
//                    {
////                        assert(temp.cross(rd).squaredNorm()==0 && "CROSS PRODUCT BETWEEN PRIMITIVE VECTORS IS NOT THE PLANE NORMAL");
//                        assert((temp+rd).squaredNorm()==0 && "CROSS PRODUCT BETWEEN PRIMITIVE VECTORS IS NOT THE PLANE NORMAL");
//                        return std::make_pair(iter.second,normMap.begin()->second);
//                    }
//                }
//            }
//            
//            std::cout<<"normal="<<rd.cartesian().transpose()<<std::endl;
//            assert(0 && "PRIMITIVE VECTORS NOT FOUND");
//        }
        
        
        /**********************************************************************/
        LatticePlaneBase(const LatticeVectorType& v1_in, const LatticeVectorType& v2_in) :
        /* init base */ ReciprocalLatticeDirectionType(v1_in,v2_in),
        /* init */ primitiveVectors(std::make_pair(v1_in,v2_in))
        {
            assert(primitiveVectors.first.squaredNorm()>0);
            assert(primitiveVectors.second.squaredNorm()>0);
            assert(primitiveVectors.first.cross(primitiveVectors.second).squaredNorm()>0);
        }
        
//        /**********************************************************************/
//        LatticePlaneBase(const ReciprocalLatticeDirectionType& rd) :
//        /* init base */ ReciprocalLatticeDirectionType(rd),
//        primitiveVectors(findPrimitiveVectors(*this))
//        {
//            assert(primitiveVectors.first.squaredNorm()>0);
//            assert(primitiveVectors.second.squaredNorm()>0);
//            assert(primitiveVectors.first.cross(primitiveVectors.second).squaredNorm()>0);
//        }
        
//        /**********************************************************************/
//        double planeSpacing() const
//        {
//            return 1.0/this->cartesian().norm();
//        }
        
        /**********************************************************************/
        LatticeVectorType snapToLattice(const Eigen::Matrix<double,3,1>& P) const
        {
            const Eigen::Matrix<double,3,2> B((Eigen::Matrix<double,2,3>()<<primitiveVectors.first.cartesian().transpose(),primitiveVectors.second.cartesian().transpose()).finished().transpose());
            const Eigen::Matrix<double,2,1> nd=(B.transpose()*B).inverse()*B.transpose()*P;
            const Eigen::Matrix<double,3,1>  p=B*RoundEigen<double,2>::round(nd);
//            return LatticeVectorType(p,primitiveVectors.first.covBasis,primitiveVectors.second.contraBasis);
            return LatticeVectorType(p,primitiveVectors.first.lattice);

        }
        
        /**********************************************************************/
        VectorDimD snapToPlane(const Eigen::Matrix<double,3,1>& P) const
        {
            //            return LatticeVectorType(p,primitiveVectors.first.covBasis,primitiveVectors.second.contraBasis);
            VectorDimD cn(this->cartesian());
            const double cnNorm(cn.norm());
            assert(cnNorm > FLT_EPSILON);
            cn/=cnNorm;
            return P-P.dot(cn)*cn;
            
        }
        
        
    };
    
} // end namespace
#endif
