/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarMeshFace_H_
#define model_PlanarMeshFace_H_

#include <chrono>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <utility>      // std::pair, std::make_pair
#include <set>

#include <SimplexTraits.h>
#include <Simplex.h>
#include <ConvexHull.h>
#include <Plane.h>



namespace model
{
    

    
    /**************************************************************************/
    /**************************************************************************/
    template<int dim>
    struct PlanarMeshFace : public StaticID<PlanarMeshFace<dim>>
    /*                   */,public std::set<const Simplex<dim,dim-1>*>
//    /*                   */,public std::vector<const Simplex<dim,dim-2>*>
    {
     
        
        
        typedef Simplex<dim,dim-1> SimplexType;
        typedef std::set<const SimplexType*> SimplexContainerType;
        const std::pair<int,int> regionIDs;

    private:
        
        Eigen::Matrix<double,dim,1> n;
        Eigen::Matrix<double,dim,1> c;
        std::vector<const Simplex<dim,0>*> _hull;
        
    public:
        
        
        /**********************************************************************/
        PlanarMeshFace(const Simplex<dim,dim-1>* const pS) :
        /* init */ regionIDs(std::make_pair(*pS->regionIDs().begin(),*pS->regionIDs().rbegin()))
        /* init */,n(pS->outNormal(regionIDs.first))
        /* init */,c(pS->center())
        {
//            std::cout<<"Creating PlanarMeshFace "<<this->sID<<", n="<<outNormal().transpose()<<std::endl;
            assert((n.norm()-1.0)<FLT_EPSILON);
            this->insert(pS);
        }
        
        
//        const std::vector<const Simplex<dim,dim-2>*>& boundaryEdges() const
//        {
//            return *this;
//        }
//
//        std::vector<const Simplex<dim,dim-2>*>& boundaryEdges()
//        {
//            return *this;
//        }
        
        const std::set<const Simplex<dim,dim-1>*>& internalSimplices() const
        {
            return *this;
        }
        
//        /**********************************************************************/
//        ~PlanarMeshFace()
//        {
//            std::cout<<"Destroying PlanarMeshFace "<<this->sID<<std::endl;
//        }
        
        /**********************************************************************/
        const std::vector<const Simplex<dim,0>*>& convexHull() const
        {
            return _hull;
        }
        
        
        
//        void walkEdge()
//        {
//            for(const auto& parent : boundaryEdges.back()->child(1).parents())
//            {
//                if(parent.second->isBoundarySimplex() || parent.second->isRegionBoundarySimplex())
//                {
//                    if(internalSimplices().find(parent.second)==internalSimplices().end())
//                    {// boundary parent not found in this face, then child is an edge
//                        boundaryEdges().insert(child);
//                        boundaryEdges().push_back(child);
//                        walkEdge();
//                    }
//                }
//            }
//        }
//
//        /**********************************************************************/
//        void buildEdges()
//        {
//            bool found
//            for(const auto& simplex : internalSimplices())
//            {
//                for(const auto& child : simplex->children())
//                {
//                        for(const auto& parent : child->parents())
//                        {
//                            if(parent.second->isBoundarySimplex() || parent.second->isRegionBoundarySimplex())
//                            {
//                                if(internalSimplices().find(parent.second)==internalSimplices().end())
//                                {// boundary parent not found in this face, then child is an edge
//                                    boundaryEdges().insert(child);
//                                    boundaryEdges().push_back(child);
//                                    walkEdge();
//                                }
//                            }
//                        }
//                }
//            }
//        }
        
        /**********************************************************************/
        void finalize()
        {/*!\todo This function uses ConvexHull and therefore will yield a wrong
          * face boundary for non-convex faces. In order to support non-convex faces
          * we need to have a structure with face edges.
          */
            
//            buildEdges();
            
            // Recompute c and n to minimize floating point errors
            n.setZero();
            c.setZero();
            for(const auto& simplex : *this)
            {
                n+=simplex->outNormal(regionIDs.first);
                c+=simplex->center();
            }
            n/=this->size();
            c/=this->size();
            assert(fabs(n.norm()-1.0)<FLT_EPSILON);
            n.normalize();
            
            for(const auto& simplex : *this)
            {// Make sure all simplices belong to face
                assert(fabs((c-simplex->center()).dot(n))<FLT_EPSILON);
            }
            
            // Compute ConvexHull of vertices
            std::set<const Simplex<dim,0>*> vertices;
            for(const auto& simplex : *this)
            {// Collect vertices on this face
                for(const auto& v : simplex->vertices())
                {
                    vertices.insert(v);
                }
            }
            
//            PlanarMeshFaceHull<dim>::makeHull(vertices,outNormal(),center,_hull);
            
            if(vertices.size())
            {
                assert(dim==3 && "ALGORITHM ONLY VALID IN dim=3");

                // Compute local rotation matrix
                Eigen::Matrix<double,dim,dim> R;
                R.col(0)=((*vertices.begin())->P0-center()).normalized();
                R.col(2)=n;
                R.col(1)=R.col(2).cross(R.col(0));

                ConvexHull<2,Simplex<dim,0>> hull;
                for(const auto& v : vertices)
                {// compute local position in the plane as insert in ConvexHull
                    const Eigen::Matrix<double,dim,1> x(R.transpose()*(v->P0-center()));
                    hull.emplace(std::array<double,2>{x[0],x[1]},v);
                }


                const auto hullPts=hull.getPoints();
                _hull.clear();
                for(const auto& hp : hullPts)
                {
                    _hull.push_back(hp.t);
                }

            }
            
            
        }
        
        /**********************************************************************/
        const SimplexContainerType& simplices() const
        {
            return *this;
        }
        
        /**********************************************************************/
        SimplexContainerType& simplices()
        {
            return *this;
        }
        
        /**********************************************************************/
        bool isExternal() const
        {
            return regionIDs.first==regionIDs.second;
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,1> outNormal() const
        {
            return isExternal()? n : Eigen::Matrix<double,dim,1>::Zero();
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,1> outNormal(const int& k) const
        {
            return k==regionIDs.first? n : (k==regionIDs.second? -n : Eigen::Matrix<double,dim,1>::Zero());
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& center() const
        {
            return c;
        }
        
        /**********************************************************************/
        Plane<dim> asPlane() const
        {
            return Plane<dim>(center(),n);
        }
        
    };
    
}
#endif
