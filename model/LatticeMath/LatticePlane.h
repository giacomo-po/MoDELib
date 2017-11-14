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

namespace model
{
//    template <int dim>
    struct LatticePlane
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef LatticeVector<dim>    LatticeVectorType;
        
       
//    public:
        
        
        /**********************************************************************/
        static long int height(const ReciprocalLatticeDirection<dim>& r,
                                                 const VectorDimD& P)
        {/*! ????
          */
            assert(r.squaredNorm()>0 && "A zero normal cannot be used as valid GlidePlane normal");
//            const VectorDimD rc();
            const double hd(P.dot(r.cartesian()));
            const long int h(std::lround(hd));
            if(fabs(hd-h)>FLT_EPSILON)
            {
                std::cout<<"rc="<<r.cartesian().transpose()<<std::endl;
                std::cout<<"P="<<P.transpose()<<std::endl;
                std::cout<<"hd="<<hd<<std::endl;
                std::cout<<"h="<<h<<std::endl;
                assert(0 && "GLIDE PLANE HEIGHT NOT FOUND");
            }
            std::cout<<"LatticePlane h="<<h<<std::endl;
            return h;
            
//            REVIEW THIS FUNCTION AND THE MATLAB FILE GENERATING THE BICRYSTAL. THE SPACING SHOULD BE 1/|r| not |r|
            
        }
        
//        const LatticeVectorType P;
        const LatticePlaneBase n;
        const long int h;
        const VectorDimD unitNormal;
        const VectorDimD P;

        /**********************************************************************/
        LatticePlane(const VectorDimD& P_in,const LatticePlaneBase& n_in) :
        /* init */ n(n_in),
        /* init */ h(height(n,P_in)),
        /* init */ unitNormal(n.cartesian().normalized()),
        /* init */ P(h*n.planeSpacing()*unitNormal)
        {
//            assert(&P.lattice==&n.lattice && "LatticeVectors have different bases.");
        }
        
        
        /**********************************************************************/
        VectorDimD snapToPlane(const VectorDimD& P0) const
        {
            return P+n.snapToPlane(P0-P);
        }
        
        /**********************************************************************/
        bool contains(const Eigen::Matrix<double,3,1>& P0) const
        {
            const double PP0((P-P0).norm());
            return PP0<FLT_EPSILON? true : (fabs((P0-P).dot(unitNormal))<FLT_EPSILON*PP0);
        }


        
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
