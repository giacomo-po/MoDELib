/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarMeshFace_cpp_
#define model_PlanarMeshFace_cpp_

#include <chrono>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <utility>      // std::pair, std::make_pair
#include <set>

#include <MeshModule.h>
#include <Plane.h>

namespace model
{

    /**********************************************************************/
    template<int dim>
    PlanarMeshFace<dim>::PlanarMeshFace(const Simplex<dim,dim-1>* const pS) :
    /* init */ regionIDs(std::make_pair(*pS->regionIDs().begin(),*pS->regionIDs().rbegin()))
    /* init */,n(pS->outNormal(regionIDs.first))
    /* init */,c(pS->center())
    /* init */,periodicFacePair(std::make_pair(VectorDim::Zero(),nullptr))
    {
        //            std::cout<<"Creating PlanarMeshFace "<<this->sID<<", n="<<outNormal().transpose()<<std::endl;
        assert((n.norm()-1.0)<FLT_EPSILON);
        this->insert(pS);
    }
    
    /**********************************************************************/
    template<int dim>
    const std::set<const Simplex<dim,dim-1>*>& PlanarMeshFace<dim>::internalSimplices() const
    {
        return *this;
    }
    
    /**********************************************************************/
    template<int dim>
    const std::vector<const Simplex<dim,0>*>& PlanarMeshFace<dim>::convexHull() const
    {
        return _hull;
    }
    
    /**********************************************************************/
    template<int dim>
    void PlanarMeshFace<dim>::finalize()
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
            if(fabs((c-simplex->center()).dot(n))>FLT_EPSILON)
            {
                throw std::runtime_error("simplex does not belong to PlanarMeshFace.");
            }
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
    template<int dim>
    const typename PlanarMeshFace<dim>:: SimplexContainerType& PlanarMeshFace<dim>::simplices() const
    {
        return *this;
    }
    
    /**********************************************************************/
    template<int dim>
    typename PlanarMeshFace<dim>:: SimplexContainerType& PlanarMeshFace<dim>::simplices()
    {
        return *this;
    }
    
    /**********************************************************************/
    template<int dim>
    bool PlanarMeshFace<dim>::isExternal() const
    {
        return regionIDs.first==regionIDs.second;
    }
    
    /**********************************************************************/
    template<int dim>
    Eigen::Matrix<double,dim,1> PlanarMeshFace<dim>::outNormal() const
    {
        return isExternal()? n : Eigen::Matrix<double,dim,1>::Zero();
    }
    
    /**********************************************************************/
    template<int dim>
    Eigen::Matrix<double,dim,1> PlanarMeshFace<dim>::outNormal(const int& k) const
    {
        return k==regionIDs.first? n : (k==regionIDs.second? (-n).eval() : Eigen::Matrix<double,dim,1>::Zero());
    }
    
    /**********************************************************************/
    template<int dim>
    const Eigen::Matrix<double,dim,1>& PlanarMeshFace<dim>::center() const
    {
        return c;
    }
    
    /**********************************************************************/
    template<int dim>
    Plane<dim> PlanarMeshFace<dim>::asPlane() const
    {
        return Plane<dim>(center(),n);
    }
    
//    template class PlanarMeshFace<1>;
    // template class PlanarMeshFace<1>;
    // template class PlanarMeshFace<2>;
    template struct PlanarMeshFace<3>;

}
#endif
