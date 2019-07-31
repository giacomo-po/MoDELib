/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpF@ucla.edu>.
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneKey_H
#define model_GlidePlaneKey_H

#include <array>
#include <Eigen/Dense>
#include <LatticeMath.h>


namespace model
{

    template<int dim>
    struct GlidePlaneKey : public std::array<long int,dim+3>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        
        const long int& grainID;
        const Eigen::Map<const VectorDimI> r;
        const long int& h;
        //
        GlidePlaneKey(const std::array<long int,dim+3>& array) :
        /* init */ std::array<long int,dim+3>(array)
        /* init */,grainID(this->operator[](0))
        /* init */,r(&this->data()[2])
        /* init */,h(this->operator[](5))
        {}
        
        GlidePlaneKey(const GlidePlaneKey& other) :
        /* init */ std::array<long int,dim+3>(other)
        /* init */,grainID(this->operator[](0))
        /* init */,r(&this->data()[2])
        /* init */,h(this->operator[](5))
        {}
        
        /**********************************************************************/
        GlidePlaneKey(const int& grainID_in,
                      const LatticePlane& lp) :
        /* init */ grainID(this->operator[](0))
        /* init */,r(&this->data()[2])
        /* init */,h(this->operator[](5))
        {/*!\param[in] grain the grain on which the GlidePlane is defined
          * \param[in] P a point on the plane
          * \param[in] N the normal to the plane
          * \returns the key which uniquely identifies the plane.
          * The type of the key is a tuple with entries (grainID,r,h), where r
          * is the ReciprocalLatticeDirection corresponding to N, and h=P.dot(r)
          * is an integer indicating the "heigth" of the plane from the origin,
          * in integer multiples of the interplanar distance d=1/|r|.
          */
            this->operator[](0)=grainID_in;
            this->operator[](1)=grainID_in;
            for(int d=0;d<dim;++d)
            {
                this->operator[](2+d)=lp.n(d); // reciprocal lattice direction components
            }
            this->operator[](2+dim)=lp.h; // integer height (units of reciprocal lattice direction)
        }
        
        /**********************************************************************/
        GlidePlaneKey(const int& grainID,
                      const VectorDimD& P,
                      const ReciprocalLatticeDirection<dim>& r) :
        /* delegate */ GlidePlaneKey(grainID,LatticePlane(P,r))
        {
        }
        
        /**********************************************************************/
        GlidePlaneKey(const int& grainID,
                      const LatticeVector<dim>& L,
                      const ReciprocalLatticeDirection<dim>& r) :
        /* delegate */ GlidePlaneKey(grainID,LatticePlane(L,r))
        {
            
        }
        
    };

}
#endif

