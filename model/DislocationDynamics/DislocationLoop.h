/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoop_H_
#define model_DislocationLoop_H_

#include <memory>

#include <model/LoopNetwork/Loop.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/LatticePlaneBase.h>
#include <model/LatticeMath/LatticePlane.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>


namespace model
{
    template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
    /*	   */ template <short unsigned int, size_t> class QuadratureRule>
    class DislocationLoop : public Loop<DislocationLoop<_dim,corder,InterpolationType,QuadratureRule>>
    {
    
    public:

        constexpr static int dim=_dim;
        typedef DislocationLoop<dim,corder,InterpolationType,QuadratureRule> DislocationLoopType;
        typedef Loop<DislocationLoopType> LoopType;
        typedef typename TypeTraits<DislocationLoopType>::LoopNetworkType LoopNetworkType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef GlidePlane<DislocationLoopType> GlidePlaneType;
        typedef GlidePlaneObserver<DislocationLoopType> GlidePlaneObserverType;
        typedef Eigen::Matrix<long int,dim+1,1> GlidePlaneKeyType;
        
        
        const Grain<dim>& grain;
        const std::shared_ptr<GlidePlaneType> _glidePlane;
        const GlidePlaneType& glidePlane;
        const bool isSessile;
        

        
        /**********************************************************************/
        DislocationLoop(const LoopNetworkType& dn,
                        const VectorDim& B,
                        const VectorDim& N,
                        const VectorDim& P,
                        const int& grainID) :
        /* base init */ LoopType(dn,dn.shared.poly.grain(grainID).latticeVector(B)),
        /*      init */ grain(dn.shared.poly.grain(grainID)),
        /*      init */ _glidePlane(GlidePlaneObserverType::getSharedPlane(dn.shared.poly.grain(grainID),P,N)),
        /*      init */ glidePlane(*_glidePlane.get()),
        /*      init */ isSessile(this->flow().dot(glidePlane.n)!=0)
        {
            
            _glidePlane->addLoop(this);
            
        }
        
        /**********************************************************************/
        ~DislocationLoop()
        {
            
            _glidePlane->removeLoop(this);
            
        }
    
    };
    
}

#endif
