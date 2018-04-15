/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePlane_h_
#define model_LatticePlane_h_

#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/LatticePlaneBase.h>
#include <model/Geometry/Plane.h>

namespace model
{
    //    template <int dim>
    struct LatticePlane //: Plane<3>
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef LatticeVector<dim>    LatticeVectorType;
        
        
        //    public:
        
        
        /**********************************************************************/
        static std::pair<bool,long int> computeHeight(const ReciprocalLatticeDirection<dim>& r,
                                                      const VectorDimD& P)
        {/*! ????
          */
            assert(r.squaredNorm()>0 && "A zero normal cannot be used as valid GlidePlane normal");
            const double hd(P.dot(r.cartesian()));
            const long int h(std::lround(hd));
            
            return  std::make_pair(fabs(hd-h)<FLT_EPSILON,h);
        }
//
//        /**********************************************************************/
//        static long int height(const ReciprocalLatticeDirection<dim>& r,
//                               const VectorDimD& P)
//        {
//            assert(r.squaredNorm()>0 && "A zero normal cannot be used as valid GlidePlane normal");
//            const double hd(P.dot(r.cartesian()));
//            const long int h(std::lround(hd));
//            if(fabs(hd-h)>FLT_EPSILON)
//            {
//                model::cout<<"P="<<P.transpose()<<std::endl;
//                model::cout<<"r="<<r.cartesian().transpose()<<std::endl;
//                model::cout<<"hd="<<std::setprecision(15)<<std::scientific<<hd<<std::endl;
//                model::cout<<"h="<<h<<std::endl;
//                assert(0 && "P in not on a lattice plane.");
//            }
//            return h;
//        }
        
        //        const LatticeVectorType P;
//        const LatticePlaneBase n;
//        const long int h;

        /**********************************************************************/
        static int sign(const long int& i)
        {
            if(i>0)
            {
                return 1;
            }
            else if(i<0)
            {
                return -1;
            }
            else
            {
                return 0;
            }
            
        }
        
        /**********************************************************************/
        static std::pair<ReciprocalLatticeDirection<dim>,size_t> get_nh(const ReciprocalLatticeDirection<dim>& r,
                                                                        const VectorDimD& P)
        {
            assert(r.squaredNorm()>0 && "A zero normal cannot be used as valid GlidePlane normal");
            const double hd(P.dot(r.cartesian()));
            const long int h(std::lround(hd));
            if(fabs(hd-h)>FLT_EPSILON)
            {
                model::cout<<"P="<<P.transpose()<<std::endl;
                model::cout<<"r="<<r.cartesian().transpose()<<std::endl;
                model::cout<<"hd="<<std::setprecision(15)<<std::scientific<<hd<<std::endl;
                model::cout<<"h="<<h<<std::endl;
                assert(0 && "P in not on a lattice plane.");
            }
//            const int signh(sign(h));
            return std::pair<ReciprocalLatticeDirection<dim>,size_t>(r*sign(h),h*sign(h));
        }
        
        const std::pair<ReciprocalLatticeDirection<dim>,size_t> nh;
        const ReciprocalLatticeDirection<dim> n;
        const size_t& h;
        //        const VectorDimD unitNormal;
//        const VectorDimD P;
        
//        /**********************************************************************/
//        LatticePlane(const VectorDimD& P_in,
//                     const LatticePlaneBase& n_in) :
////        /* init */ Plane<3>(height(n_in,P_in)*n_in.planeSpacing()*n_in.cartesian().normalized(),n_in.cartesian()),
//        /* init */ n(n_in),
//        /* init */ h(height(n,P_in))
//        ///* init */ unitNormal(n.cartesian().normalized()),
////        /* init */ P(h*n.planeSpacing()*unitNormal)
//        {
//            //            assert(&P.lattice==&n.lattice && "LatticeVectors have different bases.");
//        }
        
        /**********************************************************************/
        LatticePlane(const VectorDimD& P,
                     const ReciprocalLatticeDirection<dim>& r) :
        //        /* init */ Plane<3>(height(n_in,P_in)*n_in.planeSpacing()*n_in.cartesian().normalized(),n_in.cartesian()),
        /* init */ nh(get_nh(r,P)),
        /* init */ n(nh.first),
        /* init */ h(nh.second)
        ///* init */ unitNormal(n.cartesian().normalized()),
        //        /* init */ P(h*n.planeSpacing()*unitNormal)
        {
            //            assert(&P.lattice==&n.lattice && "LatticeVectors have different bases.");
        }
        
        
//        /**********************************************************************/
//        VectorDimD snapToPlane(const VectorDimD& P0) const
//        {
//            return P+n.snapToPlane(P0-P);
//        }
//        
//        /**********************************************************************/
//        bool contains(const Eigen::Matrix<double,3,1>& P0) const
//        {
//            const double PP0((P-P0).norm());
//            return PP0<FLT_EPSILON? true : (fabs((P0-P).dot(unitNormal))<FLT_EPSILON*PP0);
//        }
        
        //        /**********************************************************************/
        //        static long int height(const std::pair<bool,long int>& p)
        //        {
        //            assert(p.first);
        //            return p.second;
        //        }
        
        //        /**********************************************************************/
        //        LatticePlane(const LatticeVectorType& P_in,const LatticePlaneBase& n_in) :
        //        /* init */ P(P_in),
        //        /* init */ n(n_in)
        //        {
        //            assert(&P.lattice==&n.lattice && "LatticeVectors have different bases.");
        //        }
        
        //        /**********************************************************************/
        //        LatticeVectorType snapToLattice(const VectorDimD& P0) const
        //        {
        //            return P+n.snapToLattice(P0-P.cartesian());
        //        }
        //
        //        /**********************************************************************/
        //        VectorDimD snapToPlane(const VectorDimD& P0) const
        //        {
        //            return P.cartesian()+n.snapToPlane(P0-P.cartesian());
        //        }
        //
        //        /**********************************************************************/
        //        bool contains(const LatticeVectorType& L) const
        //        {
        //            assert(&P.lattice==&L.lattice && "LatticeVectors have different bases.");
        //            return (L-P).dot(n)==0;
        //        }
        //
        //        /**********************************************************************/
        //        bool contains(const Eigen::Matrix<double,3,1>& P0) const
        //        {
        //            return fabs((P0-P.cartesian()).dot(n.cartesian()))<FLT_EPSILON;
        //        }
        
    };
    
} // end namespace
#endif
