/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Simplex_cpp_
#define model_Simplex_cpp_

#include <memory> // shared_ptr
#include <set>
#include <map>
#include <vector>
#include <float.h>
#include <Eigen/Dense>

#include <MeshModule.h>

namespace model
{






template<short int dim>
Simplex<dim,0>::Simplex(SimplicialMesh<dim>* const m,
                        const SimplexIDType& vIN) :
/* init list */ SimplexBase<dim,order>(m,vIN),
//        /* init list */ P0(get_P0())
/* init list */ P0(m->simplexReader().get_P0(this->xID))
{/*!@param[in] vIN the (possibly unsorted) ID of this Simplex
  *
  * Constructur performs the following operations:
  */
    //! -1 Adds this to the SimplexObserver
    this->observer().insertSimplex(*this);
}


template<short int dim>
Simplex<dim,0>::~Simplex()
{/*!
  */
    this->observer().removeSimplex(*this);
}


template<short int dim>
std::set<const Simplex<dim,0>*> Simplex<dim,0>::boundaryNeighbors() const // boundary spouse
{
    return std::set<const Simplex<dim,0>*>();
}


template<short int dim>
Eigen::Matrix<double,dim,1> Simplex<dim,0>::outNormal() const
{
    return BoundarySimplex<dim,dim-order>::outNormal(*this);
}


template<short int dim>
Eigen::Matrix<double,dim,1> Simplex<dim,0>::outNormal(const size_t& rID) const
{
    return BoundarySimplex<dim,dim-order>::outNormal(*this,rID);
}




template<short int dim, short int order>
Simplex<dim,order>::Simplex(SimplicialMesh<dim>* const m,
                            const SimplexIDType& vIN) :
/* init */ SimplexParent<dim,order>(m,vIN)
/* init */,SimplexBase<dim,order>(m,vIN)
//        /* init */ vol0(SimplexVolume<dim,order>::volume(this->vertexPositionMatrix())) // THIS GIVES SEGMENTATION FAULT, WHY?
{/*!@param[in] vIN the (possibly unsorted) ID of this
  *
  * Constructur performs the following operations:
  */
    //! -1 inserts *this into the SimplexObserver
    this->observer().insertSimplex(*this);
    
    //! -2 inserts this into children Simplices
    for (int k=0;k<nFaces;++k)
    {
        this->child(k).addToParents(this);
    }
    
    vol0=SimplexVolume<dim,order>::volume(this->vertexPositionMatrix());
}


template<short int dim, short int order>
Simplex<dim,order>::~Simplex()
{/*! Destructor performs the following operations:
  */
    //! -1 removes this in SimplexObserver
    this->observer().removeSimplex(*this);
    
    //! -2 remove this fomr children parentContainers
    for (int k=0;k<nFaces;++k)
    {
        this->child(k).removeFromParents(this);
    }
}

//
//    template<short int dim, short int order>
//    BaseArrayType& Simplex<dim,order>::children()
//    {
//        return *this;
//    }
//
//
//    template<short int dim, short int order>
//    const BaseArrayType& Simplex<dim,order>::children() const
//    {
//        return *this;
//    }

//
//    template<short int dim, short int order>
//    ChildSimplexType& Simplex<dim,order>::child(const int& n)
//    {
//        return *(children()[n].get());
//    }
//
//
//    template<short int dim, short int order>
//    const ChildSimplexType& Simplex<dim,order>::child(const int& n) const
//    {
//        return *(children()[n].get());
//    }

//
//    template<short int dim, short int order>
//    const std::shared_ptr<ChildSimplexType>& Simplex<dim,order>::child(const ChildIDType& xID) const
//    {
//        size_t n=nFaces;
//        for (size_t k=0;k<nFaces;++k)
//        {
//            if(children()[k]->xID==xID)
//            {
//                n=k;
//                break;
//            }
//        }
//        assert(n!=nFaces && "CHILD NOT FOUND");
//        return children()[n];
//    }


template<short int dim, short int order>
size_t Simplex<dim,order>::childOrder(const ChildIDType& childID) const
{
    return SimplexTraits<dim,order>::faceOrder(this->xID,childID);
}


template<short int dim, short int order>
std::vector<int> Simplex<dim,order>::boundaryFaces() const
{
    std::vector<int> temp;
    temp.reserve(nFaces);
    for (int n=0;n<nFaces;++n)
    {
        if(this->child(n).isBoundarySimplex())
        {
            temp.push_back(n);
        }
    }
    return temp;
}



template<short int dim, short int order>
std::set<const Simplex<dim,order>*> Simplex<dim,order>::spouses() const
{/*!\return the set of Simplices sharing children with this
  */
    std::set<const Simplex<dim,order>*> temp;
    for(const auto& child : this->children())
    {
        for(const auto& parent : child->parents())
        {
            if(parent.second!=this)
            {
                temp.insert(parent.second);
            }
        }
    }
    return temp;
}


template<short int dim, short int order>
std::set<const Simplex<dim,order>*> Simplex<dim,order>::boundaryNeighbors() const // boundary spouse
{
    std::set<const Simplex<dim,order>*> temp;
    for(const auto& child : this->children())
    {
        if(child->isBoundarySimplex())
        {
            for(const auto& parent : child->parents())
            {
                if(parent.second->isBoundarySimplex())
                {
                    temp.insert(parent.second);
                }
            }
        }
    }
    return temp;
}


template<short int dim, short int order>
std::set<const Simplex<dim,order>*> Simplex<dim,order>::regionBoundaryNeighbors() const // region boundary spouse
{
    std::set<const Simplex<dim,order>*> temp;
    for(const auto& child : this->children())
    {
        if(child->isRegionBoundarySimplex())
        {
            for(const auto& parent : child->parents())
            {
                if(parent.second->isRegionBoundarySimplex())
                {
                    temp.insert(parent.second);
                }
            }
        }
    }
    return temp;
}



template<short int dim, short int order>
Eigen::Matrix<double,dim,1> Simplex<dim,order>::outNormal() const
{
    return BoundarySimplex<dim,dim-order>::outNormal(*this);
}


template<short int dim, short int order>
Eigen::Matrix<double,dim,1> Simplex<dim,order>::outNormal(const size_t& rID) const
{
    return BoundarySimplex<dim,dim-order>::outNormal(*this,rID);
}




template<short int dim>
const Eigen::Matrix<double,dim,dim+1> Simplex<dim,dim>::get_nda() const
{
    Eigen::Matrix<double,dim,dim+1> vP(this->vertexPositionMatrix());
    Eigen::Matrix<double,dim,dim> F(vP.template block<dim,dim>(0,0));
    F.colwise() -= vP.col(dim);
    const double jFabs(std::fabs(F.determinant()));
    if(jFabs<DBL_EPSILON)
    {
        std::cout<<this->xID<<", volume="<<jFabs<<std::endl;
        std::cout<<F<<std::endl;
        assert(0 && "SIMPLEX HAS ZERO VOLUME");
    }
    return jFabs*F.inverse().transpose()*BarycentricTraits<dim>::NdA;
}



template<short int dim>
Simplex<dim,dim>::Simplex(SimplicialMesh<dim>* const m,
                          const SimplexIDType& vIN, const int regionID) :
/* init */ SimplexParent<dim,order>(m,vIN)
/* init */,SimplexBase<dim,order>(m,vIN)
/* init */,region(m->getSharedRegion(regionID))
/* init */,p2b(this->b2p.fullPivLu().solve(Eigen::Matrix<double,dim+1,dim+1>::Identity()))
/* init */,nda(get_nda())
/* init */,vol0(SimplexVolume<dim,order>::volume(this->vertexPositionMatrix()))
{/*!
  */
    
    this->observer().insertSimplex(*this);
    
    for (int k=0;k<nFaces;++k)
    {
        this->child(k).addToParents(this);
    }
    
    region->simplices().emplace(this);
    
}


template<short int dim>
Simplex<dim,dim>::~Simplex()
{/*! Destructor performs the following operations:
  */
    
    //! -1 removes this in SimplexObserver
    this->observer().removeSimplex(*this);
    
    //! -2 remove this fomr children parentContainers
    for (int k=0;k<nFaces;++k)
    {
        this->child(k).removeFromParents(this);
    }
    
    region->simplices().erase(this);
    
}

//
//    template<short int dim>
//    BaseArrayType& Simplex<dim,dim>::children()
//    {
//        return *this;
//    }
//
//
//    template<short int dim>
//    const BaseArrayType& Simplex<dim,dim>::children() const
//    {
//        return *this;
//    }

//
//    template<short int dim>
//    ChildSimplexType& Simplex<dim,dim>::child(const int& n)
//    {
//        return *(children()[n].get());
//    }
//
//
//    template<short int dim>
//    const ChildSimplexType& Simplex<dim,dim>::child(const int& n) const
//    {
//        return *(children()[n].get());
//    }
//
//
//    template<short int dim>
//    const std::shared_ptr<ChildSimplexType>& Simplex<dim,dim>::child(const ChildIDType& xID) const
//    {
//        size_t n=nFaces;
//        for (size_t k=0;k<nFaces;++k)
//        {
//            if(this->operator[](k)->xID==xID)
//            {
//                n=k;
//                break;
//            }
//        }
//        assert(n!=nFaces && "CHILD NOT FOUND");
//        return this->operator[](n);
//    }


template<short int dim>
size_t Simplex<dim,dim>::childOrder(const ChildIDType& childID) const
{
    return SimplexTraits<dim,dim>::faceOrder(this->xID,childID);
}


template<short int dim>
std::vector<int> Simplex<dim,dim>::boundaryFaces() const
{
    std::vector<int> temp;
    temp.reserve(nFaces);
    for (int n=0;n<nFaces;++n)
    {
        if(this->child(n).isBoundarySimplex())
        {
            temp.push_back(n);
        }
    }
    return temp;
}


template<short int dim>
bool Simplex<dim,dim>::isBoundarySimplex() const
{
    return BoundarySimplex<dim,dim-order>::isBoundarySimplex(*this);
}


template<short int dim>
bool Simplex<dim,dim>::isRegionBoundarySimplex() const
{
    return BoundarySimplex<dim,dim-order>::isRegionBoundarySimplex(*this);
}


template<short int dim>
Eigen::Matrix<double,dim+1,1> Simplex<dim,dim>::pos2bary(const Eigen::Matrix<double,dim,1>& P) const
{
    return p2b*(Eigen::Matrix<double,dim+1,1>()<<P,1.0).finished();
}


template<short int dim>
void Simplex<dim,dim>::convexDelaunaynSearch(const bool& searchAllRegions,
                                             const Eigen::Matrix<double,dim,1>& P,
                                             std::pair<bool,const Simplex<dim,dim>*>& lastSearched,
                                             std::set<const Simplex<dim,dim>*>& searchSet) const // TO DO: searchSet is not necessary, because baryMin changes sign in next Simplex
{
    if(searchSet.find(this)==searchSet.end())
    {// this simplex has not been searched yet
        searchSet.insert(this);
        lastSearched.second=this;
        
        int kMin;
        const double baryMin=pos2bary(P).minCoeff(&kMin);
        
#ifdef _MODEL_BENCH_BARYSEARCH_
        std::cout<<"Searching "<<this->xID<<std::endl;
        std::cout<<"bary= "<<pos2bary(P)<<std::endl;
        searchFile<<this->bary2pos(Eigen::Matrix<double,dim+1,1>::Ones()/(dim+1)).transpose()<<" "
        /*      */<<this->xID<<"\n";
#endif
        
        if (baryMin>=-FLT_EPSILON)
        {
            lastSearched.first=true;
        }
        else
        {
            
            for(auto& pParent : this->child(kMin).parents())
            {
                if(pParent.second->region->regionID==region->regionID || searchAllRegions)
                {
                    pParent.second->convexDelaunaynSearch(searchAllRegions,P,lastSearched,searchSet);
                    if (lastSearched.first)
                    {
                        break;
                    }
                }
            }
        }
    }
}

template<short int dim>
Eigen::Matrix<double,dim+1,1> Simplex<dim,dim>::faceLineIntersection(const Eigen::Matrix<double,dim+1,1>& bary0,
                                                                     const Eigen::Matrix<double,dim+1,1>& bary1,
                                                                     const int& faceID) const
{/*!@param[in] bary0 barycentric coordinate of the initial point on the line
  * @param[in] P1 barycentric coordinate of the final point on the line
  * @param[in] faceID the ID of the face
  * \returns The barycentric cooridinate of the point on the line bary0->bary1
  * that intersects the facedID-face. If the line bary0->bary1 in on the face,
  * the mean point between bary0 and bary1 is returned.
  */
    
    assert((faceID>=0) && (faceID<=dim) && "0 <= faceID <= dim");
    
    //
    // Check that baricentric coordinates sum to 1
    assert(std::fabs(bary0.sum()-1.0)<=FLT_EPSILON && "bary0 must sum to 1");
    assert(std::fabs(bary1.sum()-1.0)<=FLT_EPSILON && "bary1 must sum to 1");
    
    Eigen::Matrix<double,dim+1,1> temp(bary0*0.5+bary1*0.5);
    const double den(bary1(faceID)-bary0(faceID));
    if (std::fabs(den)>DBL_EPSILON) // non-parallel points
    {
        double u(-bary0(faceID)/den); // interpolate linearly
        
        temp=bary0*(1.0-u)+bary1*u;
        temp(faceID)=0.0; // make sure
    }
    else
    {
        std::cout<<"Parallel faceLineIntersection"<<std::endl;
    }
    return temp;
    
}

template<short int dim>
std::set<size_t> Simplex<dim,dim>::regionIDs() const
{
    std::set<size_t> temp;
    temp.insert(region->regionID);
    return temp;
}

template class Simplex<3,0>;
template class Simplex<3,1>;
template class Simplex<3,2>;
template class Simplex<3,3>;

}	// close namespace
#endif
