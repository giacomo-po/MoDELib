/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarDislocationSuperLoop_H_
#define model_PlanarDislocationSuperLoop_H_


#include <memory>

#include <Loop.h>
#include <LatticeVector.h>
#include <LatticePlaneBase.h>
#include <LatticePlane.h>
#include <Grain.h>
#include <GlidePlane.h>
//#include <GlidePlaneObserver.h>

//#include <PeriodicDislocationLoopPair.h>
//#include <BoundaryLoopLinkSequence.h>
//#include <PlanarDislocationLoopIO.h>
//#include <PlanarPolygon.h>

//#ifndef NDEBUG
//#define VerbosePlanarDislocationLoop(N,x) if(verbosePlanarDislocationLoop>=N){model::cout<<x;}
//#else
//#define VerbosePlanarDislocationLoop(N,x)
//#endif


namespace model
{

//    template<typename SuperLoopType>
//    class SuperNode : public Eigen::Matrix<double,2,1>
//    {
//        static constexpr int dim=TypeTraits<LoopType>::dim;
//        typedef Eigen::Matrix<double,2,1> Vector2d;
//        typedef Eigen::Matrix<double,dim,1> VectorDim;
//
//    };
    
    template<typename SuperLoopType>
    class SuperPlane : public Eigen::Matrix<double,2,1>
    {
        static constexpr int dim=TypeTraits<LoopType>::dim;
        typedef Eigen::Matrix<double,2,1> Vector2d;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
    };
    
    template<typename SuperLoopType>
    struct SuperPlaneNode
    {
        
    };
    
    template<typename SuperLoopType>
    class SuperPlanePatch : private std::vector<std::shared_ptr<SuperPlaneNode<SuperLoopType>>>
    {
        //static constexpr int dim=TypeTraits<LoopType>::dim;
        typedef Eigen::Matrix<double,2,1> Vector2d;
        //typedef Eigen::Matrix<double,dim,1> VectorDim;

        typedef Vector2d PatchPointType;
        typedef Vector2d PatchPointType;
        typedef std::vector<PatchPointType> PatchPointTypeContainerType;
        typedef std::vector<std::shared_ptr<SuperPatchNode<SuperLoopType>>> SuperPatchNodeContainer;
        
    public:
        
        const SuperLoopType& superLoop;
        const std::shared_ptr<GLidePlane<dim>> glidePlane;
        
        SuperPatch(const std::shared_ptr<GLidePlane<dim>>& superLoop_in,
                   const SimplicialMesh<dim>& mesh,
                   const VectorDimD& P) :
        /* init */ superLoop(superLoop_in)
        /* init */ glidePlane(superLoop.glidePlane->glidePlaneObserver->sharedGlidePlane(mesh,superLoop.glidePlane->grain,P,superLoop.glidePlane->unitNormal))
        {
            
            for(const auto& lineSeg : glidePlane.meshIntersections)
            {// compute 2d points on superPlane and create patchLoop
             
                patchNodes().push_back(superLoop.sharedPatchNode(????));
                
            }
            
        }
        
        const SuperPatchNodeContainer& patchNodes() const
        {
            return *this;
        }

        SuperPatchNodeContainer& patchNodes()
        {
            return *this;
        }

        
    };
    
    
    template<int dim>
    class SuperGlidePlane
    {
        
    };

    
    template <typename LoopType>
    class PlanarDislocationSuperLoop : public std::map<Eigen::Matrix<double,2,1>,PlanarSuperNode<LoopType>,CompareVectorsByComponent<double,2,float>>

    {
        static constexpr int dim=TypeTraits<LoopType>::dim;
        typedef Eigen::Matrix<double,2,1> Vector2d;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef TypeTraits<LoopType>::NodeType NodeType;
        typedef std::vector<sta::pair<Vector2d,const NodeType* const>> Vector2dContainer;
        typedef std::map<Eigen::Matrix<double,2,1>,PlanarSuperNode<LoopType>,CompareVectorsByComponent<double,2,float>> PatchNodeContainer;

        
        /**********************************************************************/
        static MatrixDim getL2G(const VectorDim& x,
                                    const VectorDim& z)
        {
            const double xNorm(x.norm());
            const double zNorm(z.norm());
            assert(xNorm>FLT_EPSILON);
            assert(zNorm>FLT_EPSILON);
            assert(fabs(x.dot(z)<FLT_EPSILON*xNorm*zNorm));
            MatrixDim temp(Eigen::Matrix3d::Identity());
            temp.col(2)=z/zNorm;
            temp.col(0)=x/xNorm;
            temp.col(1)=temp.col(2).cross(temp.col(0));
            return temp;
        }
        
        int pnpoly(const Vector2dContainer& polygon, const Vector2d& test)
        {
            int i, j, c = 0;
            for (i = 0, j = polygon.size()-1; i < polygon.size(); j = i++)
            {
                if ( ((polygon[i](1)>test(1)) != (polygon[j](1)>test(1))) &&
                    (test(0) < (polygon[j](0)-polygon[i](0)) * (test(1)-polygon[i](1)) / (polygon[j](1)-polygon[i](1)) + polygon[i](0)) )
                {
                    c = !c;
                }
            }
            return c;
        }
        
    public:
        const GlidePlane<dim>& glidePlane;
        const MatrixDim L2G;
        
        PlanarDislocationSuperLoop(const LoopType& loop) :
        /* init */ glidePlane(*loop.glidePlane)
        /* init */,L2G(getL2G(loop.flow().cartesian().normalized(),glidePlane.unitNormal))
        {
            
            std::cout<<"Creating PlanarDislocationSuperLoop for loop "<<loop.sID<<std::endl;
            
            std::set<std::set<const PlanarMeshFace<dim>*>> faceSet;
            
            Vector2dContainer loopPoints2d;
            loopPoints2d.reserve(loop.links().size());
            for(const auto& link : loop.linkSequence())
            {
                const VectorDim pL(L2G.transpose()*(link->source()->get_P()-glidePlane.P));
//                std::cout<<pL.transpose()<<std::endl;
                assert(fabs(pL(2))<FLT_EPSILON);
                loopPoints2d.emplace_back(pL.template segment<2>(0),link->source().get());
            }
            
            
//            for(const auto& lineSeg : glidePlane.meshIntersections)
//            {// First part of algorithm:
//                if((link->source()->get_P()-lineSeg.P0).dot(lineSeg.boundaryNormal)>0.0)
//                {
//                    faceSet.insert(lineSeg.faces());
//                }
//            }
            
            Vector2dContainer meshPoints2d;
            meshPoints2d.reserve(glidePlane.meshIntersections.size());
            for(const auto& lineSeg : glidePlane.meshIntersections)
            {
                const VectorDim pL(L2G.transpose()*(lineSeg.P0-glidePlane.P));
                //                std::cout<<pL.transpose()<<std::endl;
                assert(fabs(pL(2))<FLT_EPSILON);
                meshPoints2d.emaplce_back(pL.template segment<2>(0),);
            }
            
            
            
            
        }
        
        const PatchNodeContainer& nodes() const
        {
            return *this;
        }

        PatchNodeContainer& nodes()
        {
            return *this;
        }
        
        std::shared_ptr<SuperPatchNode<SuperLoopType>> sharedPatchNode()
        {
            
        }

    };
    
    
}
#endif
