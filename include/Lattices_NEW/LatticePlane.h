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
    
    
    template<int dim>
    struct LatticePlaneKey
    {// dim ints for direction, 1 int for heigth
        
        
        typedef long int LongIntType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<LongIntType,dim,1> VectorDimI;
        typedef std::array<LongIntType,dim+2> ArrayType;
        //
        /**************************************************************************/
        static ArrayType hr2array(const VectorDimI& r, const LongIntType& h,const LongIntType& latticeID);
        
        /**************************************************************************/
        static int sgn(const LongIntType& val);
        
        /**************************************************************************/
        static ArrayType correct_h_sign(VectorDimI r,LongIntType h,const LongIntType& latticeID);
        
        
        const ArrayType array;
        
        /**********************************************************************/
        LatticePlaneKey(const VectorDimI& r,
                        const LongIntType& h,
                        const LongIntType& latticeID) ;
        
        /**********************************************************************/
        LatticePlaneKey(const VectorDimD& P,
                        const ReciprocalLatticeDirection<dim>& r) ;
        /**********************************************************************/
        LatticePlaneKey(const LatticeVector<dim>& L,
                        const ReciprocalLatticeDirection<dim>& r);
        
        /**********************************************************************/
        LatticePlaneKey(const long int& hin,
                        const ReciprocalLatticeDirection<dim>& r);
        
        
        Eigen::Map<const VectorDimI> reciprocalDirectionComponents() const;
        
        const LongIntType& planeIndex() const;
        
        const LongIntType& latticeID() const;
        
        bool operator<(const LatticePlaneKey<dim>& other) const;
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const LatticePlaneKey<dim>& key)
        {
            for(const auto& val : key.array)
            {
                os<<val<<" ";
            }
            return os;
        }
        
    };
    
    
    struct LatticePlane
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef LatticeVector<dim>    LatticeVectorType;
        
        /**********************************************************************/
        static std::pair<bool,long int> computeHeight(const ReciprocalLatticeDirection<dim>& r,
                                                      const VectorDimD& P);
        
        /**********************************************************************/
        static std::pair<bool,long int> computeHeight(const ReciprocalLatticeDirection<dim>& r,
                                                      const LatticeVector<dim>& L);
        
        const LatticePlaneKey<dim> key;
        const ReciprocalLatticeDirection<dim> n;
        const LatticePlaneKey<dim>::LongIntType planeIndex;
        
        /**********************************************************************/
        LatticePlane(const VectorDimD& P,
                     const ReciprocalLatticeDirection<dim>& r) ;
        /**********************************************************************/
        LatticePlane(const LatticeVector<dim>& L,
                     const ReciprocalLatticeDirection<dim>& r) ;
        /**********************************************************************/
        LatticePlane(const long int& hin,
                     const ReciprocalLatticeDirection<dim>& r) ;
        
        /**********************************************************************/
        VectorDimD planeOrigin() const;
        
    };
    
}
#endif

