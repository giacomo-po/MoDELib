/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePlane_h_
#define model_LatticePlane_h_

#include <iomanip>
#include <LatticeVector.h>
#include <LatticePlaneBase.h>
#include <Plane.h>

namespace model
{
    struct LatticePlane
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef LatticeVector<dim>    LatticeVectorType;
        
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
            const int signh(sign(h));
            if(signh==0)
            {
                int signr=1;
                for(int d=0;d<dim;++d)
                {
                    if(r(d)>0)
                    {
                        break;
                    }

                    if(r(d)<0)
                    {
                        signr=-1;
                        break;
                    }
                }
                return std::pair<ReciprocalLatticeDirection<dim>,size_t>(r*signr,0);

            }
            else
            {
                return std::pair<ReciprocalLatticeDirection<dim>,size_t>(r*signh,h*signh);
                
            }
        }
        
        const std::pair<ReciprocalLatticeDirection<dim>,size_t> nh;
        const ReciprocalLatticeDirection<dim>& n;
        const size_t& h;
        
        /**********************************************************************/
        LatticePlane(const VectorDimD& P,
                     const ReciprocalLatticeDirection<dim>& r) :
        /* init */ nh(get_nh(r,P)),
        /* init */ n(nh.first),
        /* init */ h(nh.second)
        {

        }
        
        /**********************************************************************/
        VectorDimD planeOrigin() const
        {
            return h*n.planeSpacing()*n.cartesian().normalized();
        }
        
    };
    
}
#endif



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
