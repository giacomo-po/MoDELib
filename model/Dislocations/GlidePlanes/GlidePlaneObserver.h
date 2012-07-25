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
#include <boost/shared_ptr.hpp>
#include <Eigen/Dense>
#include <model/Dislocations/DislocationNetworkTraits.h>
#include <model/Utilities/CompareVectorsByComponent.h>

namespace model {
	
	// Class Predeclaration
	template <typename SegmentType>
	class GlidePlane;
		
	/********************************************************************************************/
	/********************************************************************************************/	
	template<typename SegmentType>
	struct GlidePlaneObserver{
        
 //       enum{dim=SegmentType::dim};

		
		typedef GlidePlaneObserver<SegmentType> GlidePlaneObserverType;
		typedef GlidePlane<SegmentType> GlidePlaneType;
//		typedef Eigen::Matrix<double,dim  ,1> VectorDimD;
//		typedef Eigen::Matrix<double,dim+1,1> VectorDimPlusOneD;
        typedef Eigen::Matrix<double,TypeTraits<SegmentType>::dim,1> VectorDimD;
		typedef Eigen::Matrix<double,TypeTraits<SegmentType>::dim+1,1> VectorDimPlusOneD;

		
//		typedef std::map<VectorDimPlusOneD,GlidePlaneType* const,CompareVectorsByComponent<double,dim+1,float>,
//		                 Eigen::aligned_allocator<std::pair<const VectorDimPlusOneD,GlidePlaneType* const> > >  GlidePlaneMapType;
        
        typedef std::map<VectorDimPlusOneD,GlidePlaneType* const,CompareVectorsByComponent<double,TypeTraits<SegmentType>::dim+1,float>,
        /*            */ Eigen::aligned_allocator<std::pair<const VectorDimPlusOneD,GlidePlaneType* const> > > GlidePlaneMapType;

		typedef typename GlidePlaneMapType::const_iterator const_iterator;
		typedef boost::shared_ptr<GlidePlaneType> GlidePlaneSharedPtrType;
		
		/* begin() ***************************************************/
		typename GlidePlaneMapType::const_iterator begin() const {
		/*! A const iterator to the begin of the static map of GlidePlane(s)
		 */
			return glidePlaneMap.begin();
		}
		
		/* end() *****************************************************/		
		typename GlidePlaneMapType::const_iterator end() const {
		/*! A const iterator to the end of the static map of GlidePlane(s)
		 */
			return glidePlaneMap.end();
		}
		
		/* findExistingGlidePlane() **********************************/
		static GlidePlaneSharedPtrType findExistingGlidePlane(const VectorDimD& planeNormal, const double& height){
		/*! A shared pointer to an existing GlidePlane defined by planeNormal and height. 
		 *  If no GlidePlane exists, a shared pointer to a new GlidePlane is returned.
		 */
			typename GlidePlaneMapType::const_iterator iter(glidePlaneMap.find((VectorDimPlusOneD()<<planeNormal.normalized(),height).finished()));
			return (iter!=glidePlaneMap.end())? (*(iter->second->begin()))->pGlidePlane : GlidePlaneSharedPtrType(new GlidePlaneType(planeNormal,height));
		}

		/* isGlidePlane() ********************************************/
		static std::pair<bool, const GlidePlaneType* const> isGlidePlane(const VectorDimD& planeNormal, const double& height){
			typename GlidePlaneMapType::const_iterator iter(glidePlaneMap.find((VectorDimPlusOneD()<<planeNormal.normalized(),height).finished()));
			return (iter!=glidePlaneMap.end())?  std::make_pair(true,iter->second) : std::make_pair(false,(GlidePlaneType* const) NULL);
		}
		
		/* isGlidePlane() ********************************************/
//		static std::pair<bool, const GlidePlaneType* const> isGlidePlane(const Eigen::Matrix<double,dim+1,1>& key){
		static std::pair<bool, const GlidePlaneType* const> isGlidePlane(const VectorDimPlusOneD& key){
			typename GlidePlaneMapType::const_iterator iter(glidePlaneMap.find(key));
			return (iter!=glidePlaneMap.end())?  std::make_pair(true,iter->second) : std::make_pair(false,(GlidePlaneType* const) NULL);
		}
		
		/* friend T& operator << *************************************/
		template <class T>
		friend T& operator << (T& os, const GlidePlaneObserverType& gpo){
			for (typename GlidePlaneObserverType::const_iterator gpIter=gpo.begin(); gpIter!=gpo.end();++gpIter){
				os << (*gpIter->second) << "\n";
			}
			return os;
		}
		
		
	protected:
		static  GlidePlaneMapType glidePlaneMap;		
	};
	
	/**********************************************/
	/* Declare static data member *****************/	
	template<typename SegmentType>
	//std::map<Eigen::Matrix<double,dim+1,1>,GlidePlane<dim,SegmentType>* const,CompareVectorsByComponent<double,dim+1,float> > GlidePlaneObserver<dim,SegmentType>::glidePlaneMap;
	std::map<Eigen::Matrix<double,TypeTraits<SegmentType>::dim+1,1>,GlidePlane<SegmentType>* const,CompareVectorsByComponent<double,TypeTraits<SegmentType>::dim+1,float>,Eigen::aligned_allocator<std::pair<const Eigen::Matrix<double,TypeTraits<SegmentType>::dim+1,1>,GlidePlane<SegmentType>* const> > > GlidePlaneObserver<SegmentType>::glidePlaneMap;
    
	/********************************************************************************************/
	/********************************************************************************************/	
}	// close namespace
#endif

