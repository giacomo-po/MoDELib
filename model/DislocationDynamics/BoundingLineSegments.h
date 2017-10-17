/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BoundingLineSegments_H_
#define model_BoundingLineSegments_H_

#include <Eigen/Dense>
#include <deque>
#include <map>
#include <utility>
#include <model/DislocationDynamics/GlidePlanes/Glideplane.h>
#include <model/Geometry/SegmentSegmentIntersection.h>

namespace model
{
    
    template <int dim>
    struct BoundingLineSegments :
    /* */ public std::deque<std::pair<Eigen::Matrix<double,dim,1>,Eigen::Matrix<double,dim,1>>,Eigen::aligned_allocator<std::pair<Eigen::Matrix<double,dim,1>,Eigen::Matrix<double,dim,1>>>>

    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef std::pair<VectorDim,VectorDim> LineSegmentType;
        typedef std::deque<LineSegmentType,Eigen::aligned_allocator<LineSegmentType>> LineSegmentContainerType;

        
        /**********************************************************************/
        BoundingLineSegments()
        {/* Empty constructor
          */
        
        }
        
        /**********************************************************************/
        BoundingLineSegments(const BoundingLineSegments<dim>& bls1,
                             const BoundingLineSegments<dim>& bls2)
        {/* Constructor from intersectin of BoundingLineSegments
          */
            for (const auto& s1 : bls1)
            {
                for (const auto& s2 : bls2)
                {
                    SegmentSegmentIntersection<dim> ssi(s1.first,s1.second,
                                                        s2.first,s2.second);
                    
                    if(ssi.size)
                    {
                        this->emplace_back(ssi.x0,ssi.x1);
                    }
                }
            }
        }
        
        
        /**********************************************************************/
        bool contains(const VectorDim& P) const
        {
            bool temp(false);
            for(const auto& vertexPair : *this)
            {
                const VectorDim segm(vertexPair.second-vertexPair.first);
                const double segmNorm2(segm.squaredNorm());
                if(segmNorm2>FLT_EPSILON)
                {
                    double u((P-vertexPair.first).dot(segm)/segmNorm2);
                    if(u<0.0)
                    {
                        u=0.0;
                    }
                    if(u>1.0)
                    {
                        u=1.0;
                    }
                    const VectorDim x(vertexPair.first+u*segm);
                    if((P-x).squaredNorm()<FLT_EPSILON)
                    {
                        temp=true;
                        break;
                    }
                }
                else
                {
                    const VectorDim x(0.5*(vertexPair.second+vertexPair.first));
                    if((P-x).squaredNorm()<FLT_EPSILON)
                    {
                        temp=true;
                        break;
                    }
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        VectorDim snap(const VectorDim& P) const
        {
            
            assert(this->size() && "CANNOT SNAP TO EMPTY BoundingLineSegments");
            
            std::map<double,VectorDim,std::less<double>,Eigen::aligned_allocator<std::pair<double,VectorDim>>> snapMap;
            
            for(const auto& vertexPair : *this)
            {
                const VectorDim segm(vertexPair.second-vertexPair.first);
                const double segmNorm2(segm.squaredNorm());
                if(segmNorm2>FLT_EPSILON)
                {
                    double u((P-vertexPair.first).dot(segm)/segmNorm2);
                    if(u<0.0)
                    {
                        u=0.0;
                    }
                    if(u>1.0)
                    {
                        u=1.0;
                    }
                    const VectorDim x(vertexPair.first+u*segm);
                    snapMap.emplace((P-x).squaredNorm(),x);
                }
                else
                {
                    const VectorDim x(0.5*(vertexPair.second+vertexPair.first));
                    snapMap.emplace((P-x).squaredNorm(),x);
                }
            }
            
            return snapMap.begin()->second;
            
        }
        
        /**********************************************************************/
        template <typename LoopType>
        void updateWithGlidePlane(const GlidePlane<LoopType>& gp)
        {
            
            
            if(this->size())
            {
                LineSegmentContainerType temp;

                for(const auto& oldPair : *this)
                {
                    const LineSegmentContainerType psi=GlidePlaneObserver<LoopType>::planeSegmentIntersection(gp.P.cartesian(),
                                                                                                              gp.n.cartesian(),
                                                                                                              oldPair.first,
                                                                                                              oldPair.second);
                    if(psi.size())
                    {// plane and current segment intersect
                        for(size_t k=0;k<gp.meshIntersections.size();++k)
                        {
                            const size_t k1((k==gp.meshIntersections.size()-1)? 0 : k+1);
                            
                            SegmentSegmentIntersection<dim> ssi(gp.meshIntersections[k].second,
                                                                gp.meshIntersections[k1].second,
                                                                oldPair.first,
                                                                oldPair.second);
                            
                            if(ssi.size)
                            {
                                temp.emplace_back(ssi.x0,ssi.x1);
                            }
                        }
                    }
                }
                
                this->swap(temp);

            }
            else
            {
                for(size_t k=0;k<gp.meshIntersections.size();++k)
                {
                    const size_t k1((k==gp.meshIntersections.size()-1)? 0 : k+1);
                    this->emplace_back(gp.meshIntersections[k].second,gp.meshIntersections[k1].second);
                }
            }
            
        }
        
    };
    
}	// close namespace
#endif
