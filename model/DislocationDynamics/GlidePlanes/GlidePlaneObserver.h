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
        typedef Eigen::Matrix<long int,dim+2,1> GlidePlaneKeyType;
        typedef std::map<Eigen::Matrix<long int,LoopType::dim+2,1>,
        /*            */ const GlidePlane<LoopType>* const,
        /*            */ CompareVectorsByComponent<long int,LoopType::dim+2,long int>,
        /*            */ Eigen::aligned_allocator<std::pair<const Eigen::Matrix<long int,LoopType::dim+2,1>,const GlidePlane<LoopType>* const> > > GlidePlaneMapType;
        
        typedef std::shared_ptr<GlidePlaneType> GlidePlaneSharedPtrType;
        
        
    private:
        
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
            return (GlidePlaneKeyType()<<grain.grainID,n,p.dot(n)).finished();
        }
        
        /**********************************************************************/
        static std::shared_ptr<GlidePlaneType> sharedGlidePlane(const SimplicialMesh<dim>& mesh,
                                                                const Grain<dim>& grain,
                                                                const VectorDimD& P,
                                                                const VectorDimD& N)
        {
            const GlidePlaneKeyType key=getGlidePlaneKey(grain,P,N);
            const auto planeIter=glidePlaneMap.find(key);
            return (planeIter!=glidePlaneMap.end())? planeIter->second->loops().begin()->second->_glidePlane : std::make_shared<GlidePlaneType>(mesh,grain,P,N);
        }
        
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
        
//        /**********************************************************************/
//        static std::deque<std::pair<VectorDimD,VectorDimD>,Eigen::aligned_allocator<std::pair<VectorDimD,VectorDimD>>> segmentSegmentIntersection(const VectorDimD& A0,
//                                                                                                                                              const VectorDimD& B0,
//                                                                                                                                              const VectorDimD& A1,
//                                                                                                                                              const VectorDimD& B1)
//        {/*!\param[in] A0 start point of first segment
//          *\param[in] B0   end point of first segment
//          *\param[in] A1 start point of second segment
//          *\param[in] B1   end point of second segment
//          *\returns std::deque containing the starting point and the end point of the intersection segment.
//          * Start and end point are the same when the two segments intersect at one point. The container is
//          * empty if there is no intersection
//          */
//            
//            const double A0B0norm2=(B0-A0).squaredNorm();
//            const double A1B1norm2=(B1-A1).squaredNorm();
//            assert(A0B0norm2>FLT_EPSILON);
//            assert(A1B1norm2>FLT_EPSILON);
//            
//            std::deque<std::pair<VectorDimD,VectorDimD>,Eigen::aligned_allocator<std::pair<VectorDimD,VectorDimD>>> temp;
//            
//            Eigen::Matrix<double,dim,2> M;
//            M.col(0)=B0-A0;
//            M.col(1)=A1-B1;
//            
//            const Eigen::LLT<Eigen::Matrix<double,2,2>> llt(M.transpose()*M);
//            std::cout<<"DO NOT USE LLT TO SEE IF SYSTEM HAS SOLUTION. See https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html#a858dc77b65dd48248299bb6a6a758abf"<<std::endl;
//            
//            
//            if(llt.info()==Eigen::Success)
//            {
//                const Eigen::Matrix<double,2,1> u=llt.solve(M.transpose()*(A1-A0));
//                std::cout<<"u0,1="<<u(0)<<" "<<u(1)<<std::endl;
//                if(   u(0)>-FLT_EPSILON && u(0)<1.0+FLT_EPSILON
//                   && u(1)>-FLT_EPSILON && u(1)<1.0+FLT_EPSILON)
//                {
//                    const VectorDimD X0(A0+u(0)*(B0-A0));
//                    const VectorDimD X1(A1+u(1)*(B1-A1));
//                    if((X0-X1).norm()<FLT_EPSILON)
//                    {
//                        temp.emplace_back(X0,X0);
//
//                    }
//                    
//                }
//                
//            }
//            else
//            {
//                if((M.transpose()*(A1-A0)).norm()<FLT_EPSILON)
//                {// coincident segments
//                    std::multimap<double,VectorDimD> ms;
//                    
//                    ms.emplace(0.0,A0);
//                    ms.emplace(1.0,B0);
//                    const double u1=(A1-A0).dot(B0-A0)/A0B0norm2;
//                    ms.emplace(u1,A0+u1*(B0-A0));
//                    const double u2=(B1-A0).dot(B0-A0)/A0B0norm2;
//                    ms.emplace(u2,A0+u2*(B0-A0));
//                    
//                    auto iter1=ms.begin();
//                    std::advance(iter1,1);
//                    auto iter2=ms.begin();
//                    std::advance(iter2,2);
//                    temp.emplace_back(iter1->second,iter2->second);
//                }
//                
//            }
//            assert(temp.size()<=1);
//            return temp;
//        }
        
//        /**********************************************************************/
//        static std::deque<std::pair<VectorDimD,VectorDimD>,Eigen::aligned_allocator<std::pair<VectorDimD,VectorDimD>>>    lineSegmentIntersection(const VectorDimD& A0,
//                                                                                                                                              const VectorDimD& D,
//                                                                                                                                              const VectorDimD& A1,
//                                                                                                                                              const VectorDimD& B1)
//        {/*!\param[in] A0 start point line
//          *\param[in] D   direction of line
//          *\param[in] A1 start point of second segment
//          *\param[in] B1   end point of second segment
//          *\returns std::deque containing the starting point and the end point of the intersection segment.
//          * Start and end point are the same when the two segments intersect at one point. The container is
//          * empty if there is no intersection
//          */
//            
//            const double Dnorm=D.norm();
//            const double A1B1norm2=(B1-A1).squaredNorm();
//            assert(Dnorm>FLT_EPSILON);
//            assert(A1B1norm2>FLT_EPSILON);
//            
//            const VectorDimD d(D/Dnorm);
//            
//            std::deque<VectorDimD,Eigen::aligned_allocator<VectorDimD>> temp;
//            
//            Eigen::Matrix<double,dim,2> M;
//            M.col(0)=d;
//            M.col(1)=A1-B1;
//            
//            const Eigen::LLT<Eigen::Matrix<double,2,2>> llt(M.transpose()*M);
//            std::cout<<"DO NOT USE LLT TO SEE IF SYSTEM HAS SOLUTION. See https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html#a858dc77b65dd48248299bb6a6a758abf"<<std::endl;
//            
//            
//            if(llt.info()==Eigen::Success)
//            {
//                const Eigen::Matrix<double,2,1> u=llt.solve(M.transpose()*(A1-A0));
//                if( u(1)>-FLT_EPSILON && u(1)<1.0+FLT_EPSILON)
//                {
//                    const VectorDimD X(A1+u(1)*(B1-A1));
//                    temp.emplace_back(X,X);
//                    
//                }
//                
//            }
//            else
//            {
//                if((M.transpose()*(A1-A0)).norm()<FLT_EPSILON)
//                {// coincident segments
//                    temp.emplace_back(A1,B1);
//                }
//                
//            }
//            
//            assert(temp.size()<=1);
//            
//            return temp;
//        }
        
        /**********************************************************************/
        static std::deque<std::pair<VectorDimD,VectorDimD>,Eigen::aligned_allocator<std::pair<VectorDimD,VectorDimD>>> planeSegmentIntersection(const VectorDimD& P0,
                                                const VectorDimD& N,
                                                const VectorDimD& v0,
                                                const VectorDimD& v1)
        {
            //            std::deque<VectorDim> temp;
            
            const double nNorm(N.norm());
            assert(nNorm>FLT_EPSILON);
            const VectorDimD n(N/nNorm);
            
            std::deque<std::pair<VectorDimD,VectorDimD>,Eigen::aligned_allocator<std::pair<VectorDimD,VectorDimD>>> temp;
            
            
            // check intersection of v0->v1 with plane
            // x=v0+u(v1-v0)
            // (x-P0).n=0
            // (v0+u(v1-v0)-P0).n=0
            // u=(P0-v0).n/(v1-v0).n;
            const double edgeNorm=(v1-v0).norm();
            if(edgeNorm<FLT_EPSILON)
            {
            
            }
            else
            {
                const double den=(v1-v0).dot(n);
                const double num=(P0-v0).dot(n);
                const double P0v0norm=(P0-v0).norm();
                
                const double numCheck= (P0v0norm<FLT_EPSILON)? 0.0 : num/P0v0norm;
                
                if (fabs(den/edgeNorm)>FLT_EPSILON)
                {
                    // edge intersects plane
                    const double u=num/den;
                    
                    if(fabs(u)<FLT_EPSILON)
                    {
                        temp.emplace_back(v0,v0);
                    }
                    else if (u>=FLT_EPSILON && u<=1.0-FLT_EPSILON)
                    {
                        const VectorDimD x((1.0-u)*v0 + u*v1);
                        temp.emplace_back(x,x);
                    }
                    else if (fabs(1.0-u)<FLT_EPSILON)
                    {
                        temp.emplace_back(v1,v1);
                    }
                    else
                    {// no roots
                        
                    }
                    
                }
                else
                {
                    if (fabs(numCheck)>FLT_EPSILON)
                    {// edge is parallel to plane, no intersection
                        
                    }
                    else
                    {// edge is coplanar
                        temp.emplace_back(v0,v1);
//                        temp.emplace_back(v1,&edge.child(1));
                    }
                }
            }
            
            

            
            return temp;
            //            return std::make_pair(&edge,temp);
        }
        
        /**********************************************************************/
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
//			typename GlidePlaneMapType::const_iterator iter(glidePlaneMap.find((VectorDimDPlusOneD()<<planeNormal.normalized(),height).finished()));
//			return (iter!=glidePlaneMap.end())? iter->second->getSharedPointer() : GlidePlaneSharedPtrType(new GlidePlaneType(planeNormal,height));
//		}
//
//		/* isGlidePlane() ********************************************/
//		static std::pair<bool, const GlidePlaneType* const> isGlidePlane(const VectorDimD& planeNormal, const double& height)
//        {
//			typename GlidePlaneMapType::const_iterator iter(glidePlaneMap.find((VectorDimDPlusOneD()<<planeNormal.normalized(),height).finished()));
//			return (iter!=glidePlaneMap.end())?  std::make_pair(true,iter->second) : std::make_pair(false,(GlidePlaneType* const) NULL);
//		}
//
//		/* isGlidePlane() ********************************************/
//		static std::pair<bool, const GlidePlaneType* const> isGlidePlane(const VectorDimDPlusOneD& key)
//        {
//			typename GlidePlaneMapType::const_iterator iter(glidePlaneMap.find(key));
//			return (iter!=glidePlaneMap.end())?  std::make_pair(true,iter->second) : std::make_pair(false,(GlidePlaneType* const) NULL);
//		}
//

