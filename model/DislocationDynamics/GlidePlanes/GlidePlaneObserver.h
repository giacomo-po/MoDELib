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

#include <map>
#include <memory> // std::shared_ptr (c++11)
#include <Eigen/Dense>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/Utilities/CompareVectorsByComponent.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/ReciprocalLatticeVector.h>


namespace model
{
	
	// Class Predeclaration
	template <typename LoopType>
	class GlidePlane;
		
	/**************************************************************************/
	/**************************************************************************/
	template<typename LoopType>
    class GlidePlaneObserver
    {

    public:

        static constexpr int dim=LoopType::dim;
		typedef GlidePlaneObserver<LoopType> GlidePlaneObserverType;
		typedef GlidePlane<LoopType> GlidePlaneType;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;

        //typedef Eigen::Matrix<double,TypeTraits<LoopType>::dim+2,1> VectorDimPlusOneD;

        typedef Eigen::Matrix<long int,dim+2,1> GlidePlaneKeyType;
        
        
        typedef 	std::map<Eigen::Matrix<long int,LoopType::dim+2,1>,
        /*                */ const GlidePlane<LoopType>* const,
        /*                */ CompareVectorsByComponent<long int,LoopType::dim+2,long int>,
        /*                */ Eigen::aligned_allocator<std::pair<const Eigen::Matrix<long int,LoopType::dim+2,1>,const GlidePlane<LoopType>* const> > > GlidePlaneMapType;

		typedef std::shared_ptr<GlidePlaneType> GlidePlaneSharedPtrType;

        
    private:
//        
        static  GlidePlaneMapType glidePlaneMap;

    public:
        
        /**********************************************************************/
        static GlidePlaneMapType& glidePlanes()
        {
            return glidePlaneMap;
        }

        /**********************************************************************/
        static GlidePlaneKeyType getGlidePlaneKey(const Grain<dim>& grain,
                                                  const VectorDimD& P,
                                                  const VectorDimD& N)
        {
            const LatticeVector<dim> p=grain.latticeVector(P);
            const ReciprocalLatticeDirection<dim> n=grain.reciprocalLatticeDirection(N);
            assert(n.squaredNorm()>0 && "A zero normal cannot be used as valid GlidePlane key");
            //const ReciprocalLatticeVector<dim> pn=p.dot(n)*n;
            return (GlidePlaneKeyType()<<grain.grainID,n,p.dot(n)).finished();
        }
        
        /**********************************************************************/
        static std::shared_ptr<GlidePlaneType> getSharedPlane(const Grain<dim>& grain,
                                                              const VectorDimD& P,
                                                              const VectorDimD& N)
        {
            const GlidePlaneKeyType key=getGlidePlaneKey(grain,P,N);
            const auto planeIter=glidePlaneMap.find(key);
            
            return (planeIter!=glidePlaneMap.end())? planeIter->second->loops().begin()->second->_glidePlane : std::make_shared<GlidePlaneType>(grain,P,N);
        }
        
//        /**********************************************************************/
//        static GlidePlaneKeyType getGlidePlaneKey(const size_t& grainID,
//                                                  const ReciprocalLatticeVector<dim>& pn)
//        {
//            return (GlidePlaneKeyType()<<grainID,pn.transpose()).finished();
//        }
        
//        /**********************************************************************/
//        static  const GlidePlaneMapType& glidePlanes()
//        {
//            return glidePlaneMap;
//        }
        
        /**********************************************************************/
        static void addGlidePlane(const GlidePlaneType* const pL)
        {/*!@\param[in] pS a row pointer to a DislocationSegment
          * Adds pS to *this GLidePlane
          */
            const bool success=glidePlaneMap.emplace(pL->glidePlaneKey,pL).second;
            assert( success && "COULD NOT INSERT GLIDE PLANE POINTER IN GLIDE PLANE OBSERVER.");
        }
        
        /**********************************************************************/
        static void removeGlidePlane(const GlidePlaneType* const pL)
        {/*!@\param[in] pS a row pointer to a DislocationSegment
          * Removes pS from *this GLidePlane
          */
            const int success=glidePlaneMap.erase(pL->glidePlaneKey);
            assert(success==1 && "COULD NOT ERASE GLIDE PLANE POINTER FROM GLIDE PLANE OBSERVER.");
        }
//		/* begin() ***************************************************/
//		static typename GlidePlaneMapType::const_iterator begin()
//        {/*! A const iterator to the begin of the static map of GlidePlane(s)
//		  */
//			return glidePlaneMap.begin();
//		}
//		
//		/* end() *****************************************************/		
//		static typename GlidePlaneMapType::const_iterator end()
//        {/*! A const iterator to the end of the static map of GlidePlane(s)
//		  */
//			return glidePlaneMap.end();
//		}
//		
//		/* findExistingGlidePlane() **********************************/
//		static GlidePlaneSharedPtrType findExistingGlidePlane(const VectorDimD& planeNormal, const double& height)
//        {/*! A shared pointer to an existing GlidePlane defined by planeNormal and height.
//		  *  If no GlidePlane exists, a shared pointer to a new GlidePlane is returned.
//		  */
//			typename GlidePlaneMapType::const_iterator iter(glidePlaneMap.find((VectorDimPlusOneD()<<planeNormal.normalized(),height).finished()));
//			return (iter!=glidePlaneMap.end())? iter->second->getSharedPointer() : GlidePlaneSharedPtrType(new GlidePlaneType(planeNormal,height));
//		}
//
//		/* isGlidePlane() ********************************************/
//		static std::pair<bool, const GlidePlaneType* const> isGlidePlane(const VectorDimD& planeNormal, const double& height)
//        {
//			typename GlidePlaneMapType::const_iterator iter(glidePlaneMap.find((VectorDimPlusOneD()<<planeNormal.normalized(),height).finished()));
//			return (iter!=glidePlaneMap.end())?  std::make_pair(true,iter->second) : std::make_pair(false,(GlidePlaneType* const) NULL);
//		}
//		
//		/* isGlidePlane() ********************************************/
//		static std::pair<bool, const GlidePlaneType* const> isGlidePlane(const VectorDimPlusOneD& key)
//        {
//			typename GlidePlaneMapType::const_iterator iter(glidePlaneMap.find(key));
//			return (iter!=glidePlaneMap.end())?  std::make_pair(true,iter->second) : std::make_pair(false,(GlidePlaneType* const) NULL);
//		}
//		
		/* friend T& operator << *************************************/
		template <class T>
		friend T& operator << (T& os, const GlidePlaneObserverType& gpo)
        {
            for (const auto& glidePlane : gpo.glidePlanes())
            {
				os << (*glidePlane.second);
			}
			return os;
		}
		
		
	};
	
	// Static data 
	template<typename LoopType>
	std::map<Eigen::Matrix<long int,LoopType::dim+2,1>,const GlidePlane<LoopType>* const,CompareVectorsByComponent<long int,LoopType::dim+2,long int>,Eigen::aligned_allocator<std::pair<const Eigen::Matrix<long int,LoopType::dim+2,1>,const GlidePlane<LoopType>* const> > > GlidePlaneObserver<LoopType>::glidePlaneMap;
    
    
}	// close namespace
#endif

