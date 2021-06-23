/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Simplex_H_
#define model_Simplex_H_

#include <memory> // shared_ptr
#include <set>
#include <map>
#include <vector>
#include <float.h>
#include <Eigen/Dense>
#include <SimplexTraits.h>
#include <SimplexReader.h>
#include <SimplexObserver.h>
#include <SimplexBase.h>
#include <SimplexParent.h>
#include <SimplexChild.h>
#include <MeshRegion.h>
//#include <VertexReader.h>
#include <BarycentricTraits.h>
#include <SimplexVolume.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim>
    class Simplex<dim,0> :
    /* inheritance      */ public SimplexBase<dim,0>,
    /* inheritance      */ public SimplexChild <dim,0>
    {
        
    public:
        
        
        enum{order=0};
        enum{nVertices=SimplexTraits<dim,order>::nVertices};
        enum{nFaces=SimplexTraits<dim,order>::nFaces};
        typedef typename SimplexTraits<dim,order>::SimplexIDType SimplexIDType;
        
        const Eigen::Matrix<double,dim,1> P0;
        
        Simplex(SimplicialMesh<dim>* const m,
                const SimplexIDType& vIN);
        ~Simplex();
        std::set<const Simplex<dim,0>*> boundaryNeighbors() const;
        Eigen::Matrix<double,dim,1> outNormal() const;
        Eigen::Matrix<double,dim,1> outNormal(const size_t& rID) const;
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim, short int order>
    class Simplex : public SimplexParent<dim,order>
    /*          */, public SimplexBase<dim,order>
    /*          */, public SimplexChild <dim,order>
    {
        
        typedef typename SimplexTraits<dim,order>::BaseArrayType BaseArrayType;
        
        
    public:
        
        
        //        enum{nVertices=SimplexTraits<dim,order>::nVertices};
        static constexpr int nVertices=SimplexTraits<dim,order>::nVertices;
        static constexpr int nFaces=SimplexTraits<dim,order>::nFaces;
        //        enum{nFaces=SimplexTraits<dim,order>::nFaces};
        typedef typename SimplexTraits<dim,order>::SimplexIDType SimplexIDType;
        
        typedef Simplex<dim,order-1> ChildSimplexType;
        typedef typename SimplexTraits<dim,order-1>::SimplexIDType ChildIDType;
        
        double vol0; // SHOULD BE CONST, SEE BELOW
        
        /**********************************************************************/
        Simplex(SimplicialMesh<dim>* const m,
                const SimplexIDType& vIN);
        ~Simplex();
        //        BaseArrayType& children();
        //        const BaseArrayType& children() const;
        //        ChildSimplexType& child(const int& n);
        //        const ChildSimplexType& child(const int& n) const;
        //        const std::shared_ptr<ChildSimplexType>& child(const ChildIDType& xID) const;
        size_t childOrder(const ChildIDType& childID) const;
        std::vector<int> boundaryFaces() const;
        std::set<const Simplex<dim,order>*> spouses() const;
        std::set<const Simplex<dim,order>*> boundaryNeighbors() const;
        std::set<const Simplex<dim,order>*> regionBoundaryNeighbors() const;
        Eigen::Matrix<double,dim,1> outNormal() const;
        Eigen::Matrix<double,dim,1> outNormal(const size_t& rID) const;
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim>
    class Simplex<dim,dim> : public SimplexParent<dim,dim>
    /*                   */, public SimplexBase<dim,dim>
    {
        
        typedef typename SimplexTraits<dim,dim>::BaseArrayType BaseArrayType;
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,dim+1> get_nda() const;
        
    public:
        
        
        enum{order=dim};
        enum{nVertices=SimplexTraits<dim,order>::nVertices};
        enum{nFaces=SimplexTraits<dim,order>::nFaces};
        typedef typename SimplexTraits<dim,order>::SimplexIDType SimplexIDType;
        
        typedef Simplex<dim,order-1> ChildSimplexType;
        typedef typename SimplexTraits<dim,order-1>::SimplexIDType ChildIDType;
        
        typedef MeshRegion<dim> MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        
        //! Shared pointer to MeshRegioin containing *this
        const std::shared_ptr<MeshRegionType> region;
        
        //! The position to barycentric-coordinate transformation matrix
        const Eigen::Matrix<double,dim+1,dim+1> p2b;
        
        //! The column matrix of face normals
        const Eigen::Matrix<double,dim,dim+1> nda;
        
        const double vol0;
        
        Simplex(SimplicialMesh<dim>* const m,
                const SimplexIDType& vIN, const int regionID=0);
        ~Simplex();
        //        BaseArrayType& children();
        //        const BaseArrayType& children() const;
        //        ChildSimplexType& child(const int& n);
        //        const ChildSimplexType& child(const int& n) const;
        //        const std::shared_ptr<ChildSimplexType>& child(const ChildIDType& xID) const;
        size_t childOrder(const ChildIDType& childID) const;
        std::vector<int> boundaryFaces() const;
        bool isBoundarySimplex() const;
        bool isRegionBoundarySimplex() const;
        Eigen::Matrix<double,dim+1,1> pos2bary(const Eigen::Matrix<double,dim,1>& P) const;
        void convexDelaunaynSearch(const bool& searchAllRegions,
                                   const Eigen::Matrix<double,dim,1>& P,
                                   std::pair<bool,const Simplex<dim,dim>*>& lastSearched,
                                   std::set<const Simplex<dim,dim>*>& searchSet) const; // TO DO: searchSet is not necessary, because baryMin changes sign in next Simplex
        Eigen::Matrix<double,dim+1,1> faceLineIntersection(const Eigen::Matrix<double,dim+1,1>& bary0,
                                                           const Eigen::Matrix<double,dim+1,1>& bary1,
                                                           const int& faceID) const;
        std::set<size_t> regionIDs() const;
        
    };
    
    
}	// close namespace
#endif
