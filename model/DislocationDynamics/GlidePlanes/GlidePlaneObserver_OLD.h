/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GLIDEPLANEOBSERVER_H_
#define model_GLIDEPLANEOBSERVER_H_

#include <algorithm>
#include <map>
#include <tuple>
#include <memory> // std::shared_ptr (c++11)
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <SimplicialMesh.h>
//#include <DislocationNetworkTraits.h>
#include <CompareVectorsByComponent.h>
#include <LatticeVector.h>
#include <ReciprocalLatticeDirection.h>
#include <PlanePlaneIntersection.h>
#include <MeshPlane.h>


namespace model
{
    
    
    template<typename KeyType,typename ValueType>
    struct KeyConstructableSharedPtrMap : private std::map<KeyType,const std::shared_ptr<const ValueType>>
    {
        
        
        std::shared_ptr<const ValueType> operator[](const KeyType& key)
        {
            const auto iter(this->find(key));
            if(iter!=this->end())
            {// key exists
                
            }
            
        }
        
    };
    
    
    template <int dim>
    struct GlidePlane; // class predeclaration
    
    /**************************************************************************/
    template<int dim>
    struct GlidePlaneObserver : private std::map<std::array<long int,dim+3>,const GlidePlane<dim>* const>
    /*                       */,private std::map<std::pair<size_t,size_t>,PlanePlaneIntersection<dim>>
    {

        typedef GlidePlaneObserver<dim> GlidePlaneObserverType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef MeshPlane<dim> MeshPlaneType;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef typename GlidePlaneType::GlidePlaneKeyType GlidePlaneKeyType;
        typedef std::map<GlidePlaneKeyType,const GlidePlane<dim>* const> GlidePlaneMapType;
        typedef std::shared_ptr<GlidePlaneType> GlidePlaneSharedPtrType;
        typedef PlanePlaneIntersection<dim> PlanePlaneIntersectionType;
        typedef std::map<std::pair<size_t,size_t>,PlanePlaneIntersectionType> MeshPlaneIntersectionContainerType;
        
    public:
        
        /**********************************************************************/
        MeshPlaneIntersectionContainerType& glidePlaneIntersections()
        {
            return *this;
        }
        
        /**********************************************************************/
        const MeshPlaneIntersectionContainerType& glidePlaneIntersections() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const PlanePlaneIntersectionType& glidePlaneIntersection(const MeshPlaneType* const p1,
                                                                 const MeshPlaneType* const p2)
        {/*@param[in] p1 first plane
          *@param[in] p2 second plane
          *\returns the infinite line of interseciton between the two planes,
          * if already computed. Otherwise it computes the interseciton, stores it,
          * and returns it.
          */
            const auto key=std::make_pair(std::max(p1->sID,p2->sID),std::min(p1->sID,p2->sID));
            const auto iter=glidePlaneIntersections().find(key);
            if(iter==glidePlaneIntersections().end())
            {
                const bool success=glidePlaneIntersections().emplace(std::piecewise_construct,
                                                                     std::make_tuple(key),
                                                                     std::make_tuple(p1->P,p1->unitNormal,p2->P,p2->unitNormal)
                                                                     ).second;
                assert(success);
            }
            return glidePlaneIntersections().at(key);
        }
        
        /**********************************************************************/
        GlidePlaneMapType& glidePlanes()
        {
            return *this;
        }
        
        /**********************************************************************/
        const GlidePlaneMapType& glidePlanes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        std::shared_ptr<GlidePlaneType> sharedGlidePlane(const SimplicialMesh<dim>& mesh,
                                                         const Grain<dim>& grain,
                                                         const VectorDimD& P,
                                                         const VectorDimD& N)
        {
            const LatticePlane lp(P,grain.reciprocalLatticeDirection(N));
            const GlidePlaneKeyType key=GlidePlaneType::getGlidePlaneKey(grain.grainID,lp);
            const auto planeIter=glidePlanes().find(key);
            return (planeIter!=glidePlanes().end())? planeIter->second->sharedPlane() : std::shared_ptr<GlidePlaneType>(new GlidePlaneType(this,mesh,grain,P,N));
            //            return (planeIter!=glidePlanes().end())? planeIter->second->loops().begin()->second->_glidePlane :
            //            /*                            */ std::shared_ptr<GlidePlaneType>(new GlidePlaneType(this,mesh,grain,P,N));
        }
        
        /**********************************************************************/
        void addGlidePlane(const GlidePlaneType* const pL)
        {/*!@\param[in] pS a row pointer to a DislocationSegment
          * Adds pS to *this GLidePlane
          */
            const bool success=glidePlanes().emplace(pL->glidePlaneKey,pL).second;
            assert( success && "COULD NOT INSERT GLIDE PLANE POINTER IN GLIDE PLANE OBSERVER.");
        }
        
        /**********************************************************************/
        void removeGlidePlane(const GlidePlaneType* const pL)
        {/*!@\param[in] pS a row pointer to a DislocationSegment
          * Removes pS from *this GLidePlane
          */
            const int success=glidePlanes().erase(pL->glidePlaneKey);
            assert(success==1 && "COULD NOT ERASE GLIDE PLANE POINTER FROM GLIDE PLANE OBSERVER.");
        }

    };
}
#endif


