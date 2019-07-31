/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SuperGlidePlane_H_
#define model_SuperGlidePlane_H_


#include<Eigen/StdVector>
#include <memory>


#include <GlidePlane.h>
#include <GlidePlaneObserver.h>




namespace model
{

    template<int dim>
    class SuperGlidePlane;

    template<int dim>
    class SuperPlanePatch;

    template<int dim>
    class SuperPlaneNode;

    template<int dim>
    class SuperPlaneEdge;

    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct SuperPlaneEdge
    {
        
        const SuperPlanePatch<dim>* const patch;
        const std::shared_ptr<SuperPlaneNode<dim>> source;
        const std::shared_ptr<SuperPlaneNode<dim>>   sink;
        SuperPlaneEdge<dim>* next;
        SuperPlaneEdge<dim>* prev;
        
        /**********************************************************************/
        SuperPlaneEdge(const SuperPlanePatch<dim>* const patch_in,
                       const std::shared_ptr<SuperPlaneNode<dim>>& source_in,
                       const std::shared_ptr<SuperPlaneNode<dim>>& sink_in);
        ~SuperPlaneEdge();
        
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct SuperPlaneNodeConnectivity
    {
        typedef SuperPlaneEdge<dim> SuperPlaneEdgeType;

        SuperPlaneEdgeType* inEdge;
        SuperPlaneEdgeType* outEdge;
        
        SuperPlaneNodeConnectivity() :
        /* init */ inEdge(nullptr)
        /* init */,outEdge(nullptr)
        {
            
        }

    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct SuperPlaneNode : public StaticID<SuperPlaneNode<dim>>
    /*                   */,public Eigen::Matrix<double,dim-1,1>
    /*                   */,public std::map<size_t,SuperPlaneNodeConnectivity<dim>>
    {
        
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef SuperPlaneNode<dim> SuperPlaneNodeType;
        typedef SuperPlaneEdge<dim> SuperPlaneEdgeType;
        typedef SuperPlaneNodeConnectivity<dim> NeighborType;
        typedef std::map<size_t,SuperPlaneNodeConnectivity<dim>> NeighborContainerType;

        /**********************************************************************/
        SuperPlaneNode(const VectorLowerDim& pos) :
        /* init*/ VectorLowerDim(pos)
        {
            
        }
        
        /**********************************************************************/
        void addLink(SuperPlaneEdgeType* const link)
        {
            NeighborType& neighbor(neighbors()[link->patch->sID]);
            if(link->sink->sID==this->sID)
            {// edge ends at this node
                assert(neighbor.inEdge==nullptr);
                neighbor.inEdge=link;
                if(neighbor.outEdge)
                {// and outEdge of the same patch is connected to this node
                    assert(neighbor.outEdge->prev=nullptr);
                    assert(link->next=nullptr);
                    neighbor.outEdge->prev=link;
                    link->next=neighbor.outEdge;
                }
            }
            else if(link->source->sID==this->sID)
            {// edge starts at this node
                assert(neighbor.outEdge==nullptr);
                neighbor.outEdge=link;
                if(neighbor.inEdge)
                {
                    assert(neighbor.inEdge->next==nullptr);
                    assert(link->prev==nullptr);
                    neighbor.inEdge->next=link;
                    link->prev=neighbor.inEdge;
                }
            }
            else
            {
                assert(false && "CONNECTING LINK TO NON-INCIDENT NODE");
            }
            
        }
        
        /**********************************************************************/
        void removeLink(SuperPlaneEdgeType* const link)
        {
            auto iter(neighbors().find(link->patch->sID));
            if(iter!=neighbors().end())
            {
                NeighborType& neighbor(iter->second);
                if(link->sink->sID==this->sID)
                {// edge ends at this node
                    assert(neighbor.inEdge==link);
                    neighbor.inEdge=nullptr;
                }
                else if(link->source->sID==this->sID)
                {// edge starts at this node
                    assert(neighbor.outEdge==link);
                    neighbor.outEdge=nullptr;
                }
                else
                {
                    assert(false && "DISCONNECTING LINK FROM NON-INCIDENT NODE");
                }
                
                if(neighbor.inEdge==nullptr && neighbor.outEdge==nullptr)
                {
                    neighbors().erase(iter);
                }
            }
        }

        /**********************************************************************/
        const NeighborContainerType& neighbors() const
        {
            return *this;
        }
        
        /**********************************************************************/
        NeighborContainerType& neighbors()
        {
            return *this;
        }
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    struct SuperPlanePatch : public StaticID<SuperPlanePatch<dim>>
    /*                    */,private std::vector<SuperPlaneEdge<dim>>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        const SuperGlidePlane<dim>* const superPlane;
        const VectorDim shift;
        const std::shared_ptr<GlidePlane<dim>> glidePlane;
        
        /**********************************************************************/
        SuperPlanePatch(const SuperGlidePlane<dim>* const superPlane_in,
                        const VectorDim& shift_in) :
        /* init */ superPlane(superPlane_in)
        /* init */,shift(shift_in)
        /* init */,glidePlane(superPlane->get()->glidePlaneObserver->sharedGlidePlane(superPlane->mesh,superPlane->get()->grain,superPlane->get()->P+shift,superPlane->get()->unitNormal))
        {
            
            
            for(size_t k0=0;k0<glidePlane->meshIntersections.size();++k0)
            {// compute 2d points on superPlane, and coneect them
                size_t k1(k0!=glidePlane->meshIntersections.size()-1? k0+1 : 0);
                std::shared_ptr<SuperPlaneNode<dim>> source(superPlane->getSharedNode(glidePlane->meshIntersections[k0].P0-shift));
                std::shared_ptr<SuperPlaneNode<dim>>   sink(superPlane->getSharedNode(glidePlane->meshIntersections[k1].P0-shift));
                this->emplace_back(this,source,sink);
            }
            
        }
        
        /**********************************************************************/
        const std::vector<SuperPlaneEdge<dim>>& edges() const
        {
            return *this;
        }
        
    };
    
    /**********************************************************************/
    /**********************************************************************/
    template<int dim>
    class SuperGlidePlane : public  std::shared_ptr<GlidePlane<dim>>
    /*                   */,private std::map<Eigen::Matrix<double,dim-1,1>,SuperPlaneNode<dim>* const,CompareVectorsByComponent<double,dim-1,float>>
    /*                   */,public  std::vector<SuperPlanePatch<dim>>
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef std::vector<SuperPlanePatch<dim>> PatchContainerType;
        typedef std::map<Eigen::Matrix<double,dim-1,1>,SuperPlaneNode<dim>* const,CompareVectorsByComponent<double,dim-1,float>> NodeCointainerType;
        
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

    public:
        
        const SimplicialMesh<dim>& mesh;
        const MatrixDim L2G;
        
        SuperGlidePlane(const SimplicialMesh<dim>& mesh_in,
                        const std::shared_ptr<GlidePlane<dim>>& glidePlane_in,
                        const VectorDim& globalX) :
        /* init */ std::shared_ptr<GlidePlane<dim>>(glidePlane_in)
        /* init */,mesh(mesh_in)
        /* init */,L2G(getL2G(globalX,this->get()->unitNormal))
        {
        
        }
        
        /**********************************************************************/
        const PatchContainerType& patches() const
        {
            return *this;
        }
        
        PatchContainerType& patches()
        {
            return *this;
        }
        
        /**********************************************************************/
        NodeCointainerType& nodes()
        {
            return *this;
        }
        
        const NodeCointainerType& nodes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        void addPatch(const VectorDim& shift)
        {
            SuperPlanePatch patch1(this,shift);
            SuperPlanePatch patch2(this,shift);

            std::cout<<*patch1.glidePlane<<std::endl;
                        std::cout<<*patch2.glidePlane<<std::endl;
//            return patches().push_back();
        }
        
        /**********************************************************************/
        VectorLowerDim getLocalPosition(const VectorDim& point) const
        {
            const VectorDim pointLocal(L2G.transpose()*(point-this->get()->P));
            assert(fabs(pointLocal(2))<FLT_EPSILON);
            return pointLocal.template segment<2>(0);
        }
        
//        /**********************************************************************/
//        std::shared_ptr<SuperPlaneNode<dim> > getSharedNode(const VectorDim& point) const
//        {
//            return getSharedNode(getLocalPosition(point));
//        }
        
        /**********************************************************************/
        std::shared_ptr<SuperPlaneNode<dim> > getSharedNode(const VectorDim& pointDim) const
        {
            const VectorLowerDim point(getLocalPosition(pointDim));
            const auto iter(nodes().find(point));
            if(iter==nodes().end())
            {// point does not exist
                return std::shared_ptr<SuperPlaneNode<dim>>(new SuperPlaneNode<dim>(point));
            }
            else
            {
                assert(iter->second->neighbors().size()>0 && "EXISTING NODE IS ISOLATED");
                const auto neighbor(iter->second->neighbors().begin()->second);
                if(neighbor.outEdge)
                {
                    return neighbor.outEdge->source;
                }
                else
                {
                    if(neighbor.inEdge)
                    {
                        return neighbor.inEdge->sink;
                    }
                    else
                    {
                        assert(false && "EXISTING POSITION IS NOT CONNECTED");
                        return std::shared_ptr<SuperPlaneNode<dim> >(nullptr);
                    }
                }
            }
        }
        
    };
    
    
    
    /**********************************************************************/
    template<int dim>
    SuperPlaneEdge<dim>::SuperPlaneEdge(const SuperPlanePatch<dim>* const patch_in,
                   const std::shared_ptr<SuperPlaneNode<dim>>& source_in,
                   const std::shared_ptr<SuperPlaneNode<dim>>& sink_in) :
    /* init */ patch(patch_in)
    /* init */,source(source_in)
    /* init */,sink(sink_in)
    /* init */,next(nullptr)
    /* init */,prev(nullptr)
    {
        source->addLink(this);
        sink->addLink(this);
    }
    
    /**********************************************************************/
    template<int dim>
    SuperPlaneEdge<dim>::~SuperPlaneEdge()
    {
        source->removeLink(this);
        sink->removeLink(this);
        if(next)
        {
            next->prev=nullptr;
        }
        if(prev)
        {
            prev->next=nullptr;
        }
    }
    
}
#endif
