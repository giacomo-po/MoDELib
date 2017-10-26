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
#include <model/DislocationDynamics/IO/DislocationLoopIO.h>


namespace model
{
    template <int _dim, short unsigned int corder, typename InterpolationType>
    class DislocationLoop : public Loop<DislocationLoop<_dim,corder,InterpolationType>>
    {
    
    public:

        constexpr static int dim=_dim;
        typedef DislocationLoop<dim,corder,InterpolationType> DislocationLoopType;
        typedef Loop<DislocationLoopType> BaseLoopType;
        typedef typename BaseLoopType::LoopLinkType LoopLinkType;
        typedef typename TypeTraits<DislocationLoopType>::LoopNetworkType LoopNetworkType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef GlidePlane<DislocationLoopType> GlidePlaneType;
        typedef GlidePlaneObserver<DislocationLoopType> GlidePlaneObserverType;
        typedef Eigen::Matrix<long int,dim+1,1> GlidePlaneKeyType;
        
        
        const Grain<dim>& grain;
        const std::shared_ptr<GlidePlaneType> _glidePlane;
        const GlidePlaneType& glidePlane;
        const bool isGlissile;
        
        
        /**********************************************************************/
        static bool allowedSlipSystem(const LatticeVector<dim>& b,
                                      const ReciprocalLatticeDirection<dim>& n,
                                      const Grain<dim>& gr)
        {
        
            std::cout<<"DislocationLoop::FINISH HERE. ANOTHER CONDITION FOR isGlissile SHOULD BE THAT N IS ONE OF THE ALLOWED SLIP PLANES"<<std::endl;

            return true;
        }

        
        /**********************************************************************/
        DislocationLoop(LoopNetworkType* const dn,
                        const VectorDim& B,
                        const VectorDim& N,
                        const VectorDim& P,
                        const int& grainID) :
        /* base init */ BaseLoopType(dn,dn->poly.grain(grainID).latticeVector(B)),
        /*      init */ grain(dn->poly.grain(grainID)),
        /*      init */ _glidePlane(dn->sharedGlidePlane(dn->mesh,grain,P,N)),
        /*      init */ glidePlane(*_glidePlane.get()),
        /*      init */ isGlissile(this->flow().dot(glidePlane.n)==0 && allowedSlipSystem(this->flow(),glidePlane.n,grain))
        {
            model::cout<<"Creating DislocationLoop "<<this->sID<<std::endl;
            
            _glidePlane->addLoop(this);
        }
        
        /**********************************************************************/
        DislocationLoop(const DislocationLoop& other) :
        /* base init */ BaseLoopType(other),
        /* init */ grain(other.grain),
        /* init */ _glidePlane(other._glidePlane),
        /* init */ glidePlane(*_glidePlane.get()),
        /* init */ isGlissile(other.isGlissile)
        {
            model::cout<<"Copying DislocationLoop "<<this->sID<<std::endl;
            
            _glidePlane->addLoop(this);
        }
        
        /**********************************************************************/
        ~DislocationLoop()
        {
            
            _glidePlane->removeLoop(this);
            
        }
        
        /**********************************************************************/
        void addLink(LoopLinkType* const pL)
        {
            BaseLoopType::addLink(pL);

            assert(std::fabs((pL->source()->get_P()-pL->sink()->get_P()).dot(glidePlane.n.cartesian()))<FLT_EPSILON && "Chord does not belong to plane");
            
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationLoopType& dL)
        {
            os<< DislocationLoopIO<dim>(dL);
            return os;
        }
        
//        /**********************************************************************/
//        void removeLink(LoopLinkType* const pL)
//        {
//            LoopType::removeLink(pL);
//
//        }
    
    };
    
}

#endif
