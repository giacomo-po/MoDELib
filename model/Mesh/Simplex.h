/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Simplex_H_
#define model_Simplex_H_

#include <set>
#include <map>
#include <Eigen/Dense>
#include <model/Mesh/SimplexTraits.h>
#include <model/Mesh/SimplexObserver.h>
#include <model/Mesh/SimplexBase.h>
#include <model/Mesh/SimplexChild.h>

namespace model {
    
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
        
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,1> get_P0() const
        {
            const typename VertexReader<'N',dim+1,double>::const_iterator nIter(nodeReader.find((this->xID)(0)));
            assert(nIter!=nodeReader.end() && "MESH VERTEX NOT FOUND IN N_0.txt.");
            return nIter->second;
        }
        
    public:
        
        static VertexReader<'N',dim+1,double> nodeReader;
        
		/**********************************************************************/
        Simplex(const SimplexIDType& vIN) :
        /* init list */ SimplexBase<dim,order>(vIN),
        /* init list */ P0(get_P0())
        {/*!@param[in] vIN the (possibly unsorted) ID of this Simplex
          *
          * Constructur performs the following operations:
          */
            //! -1 Adds this to the SimplexObserver
            SimplexObserver<dim,order>::insertSimplex(*this);
        }
        
        /**********************************************************************/
        ~Simplex()
        {/*!
          */
            SimplexObserver<dim,order>::removeSimplex(*this);
        }
        
	};
    
    template<short int dim>
    VertexReader<'N',dim+1,double> Simplex<dim,0>::nodeReader;
    
    
    /**************************************************************************/
	/**************************************************************************/
	template<short int dim, short int order>
	class Simplex :
    /* inheritance */ public SimplexBase<dim,order>,
    /* inheritance */ public SimplexTraits<dim,order>::BaseArrayType,
    /* inheritance */ public SimplexChild <dim,order>
    {
        
        typedef typename SimplexTraits<dim,order>::BaseArrayType BaseArrayType;
        
        
    public:
        enum{nVertices=SimplexTraits<dim,order>::nVertices};
        enum{nFaces=SimplexTraits<dim,order>::nFaces};
        typedef typename SimplexTraits<dim,order>::SimplexIDType SimplexIDType;
        
        typedef Simplex<dim,order-1> ChildSimplexType;
        typedef typename SimplexTraits<dim,order-1>::SimplexIDType ChildIDType;
        
        
    public:
		/**********************************************************************/
        Simplex(const SimplexIDType& vIN) :
        /* init list */ SimplexBase<dim,order>(vIN),
        /* init list */ BaseArrayType(SimplexObserver<dim,order>::faces(vIN))
        {/*!@param[in] vIN the (possibly unsorted) ID of this
          *
          * Constructur performs the following operations:
          */
            //! -1 inserts *this into the SimplexObserver
            SimplexObserver<dim,order>::insertSimplex(*this);
            
            //! -2 inserts this into children Simplices
            for (int k=0;k<nFaces;++k)
            {
                this->child(k).addToParents(this);
            }
        }
        
        /**********************************************************************/
        ~Simplex()
        {/*! Destructor performs the following operations:
          */
            //! -1 removes this in SimplexObserver
            SimplexObserver<dim,order>::removeSimplex(*this);
            
            //! -2 remove this fomr children parentContainers
            for (int k=0;k<nFaces;++k)
            {
                this->child(k).removeFromParents(this);
            }
        }
        
        /**********************************************************************/
        ChildSimplexType& child(const int& n)
        {
            return *(this->operator[](n).get());
        }
        
        /**********************************************************************/
        const ChildSimplexType& child(const int& n) const
        {
            return *(this->operator[](n).get());
        }
        
        /**********************************************************************/
        const std::shared_ptr<ChildSimplexType>& child(const ChildIDType& xID) const
        {
            size_t n=nFaces;
            for (size_t k=0;k<nFaces;++k)
            {
                if(this->operator[](k)->xID==xID)
                {
                    n=k;
                    break;
                }
            }
            assert(n!=nFaces && "CHILD NOT FOUND");
            return this->operator[](n);
        }
        
        /**********************************************************************/
        std::vector<int> boundaryFaces() const
        {
            std::vector<int> temp;
            temp.reserve(nFaces);
            for (int n=0;n<nFaces;++n)
            {
                if(child(n).isBoundarySimplex())
                {
                    temp.push_back(n);
                }
            }
            return temp;
        }
        
	};
    
    /**************************************************************************/
	/**************************************************************************/
	template<short int dim>
	class Simplex<dim,dim> :
    /* inheritance */ public SimplexBase<dim,dim>,
    /* inheritance */ public SimplexTraits<dim,dim>::BaseArrayType
    {
        
        typedef typename SimplexTraits<dim,dim>::BaseArrayType BaseArrayType;
        
    public:
        enum{order=dim};
        enum{nVertices=SimplexTraits<dim,order>::nVertices};
        enum{nFaces=SimplexTraits<dim,order>::nFaces};
        typedef typename SimplexTraits<dim,order>::SimplexIDType SimplexIDType;
        
        typedef Simplex<dim,order-1> ChildSimplexType;
        typedef typename SimplexTraits<dim,order-1>::SimplexIDType ChildIDType;
        
        
		/**********************************************************************/
        Simplex(const SimplexIDType& vIN) :
        /* init list */ SimplexBase<dim,order>(vIN),
        /* init list */ BaseArrayType(SimplexObserver<dim,dim>::faces(vIN))
        {/*!
          */
            SimplexObserver<dim,order>::insertSimplex(*this);
            
            for (int k=0;k<nFaces;++k)
            {
                this->child(k).addToParents(this);
            }
        }
        
        /**********************************************************************/
        ~Simplex()
        {/*! Destructor performs the following operations:
          */
            //! -1 removes this in SimplexObserver
            SimplexObserver<dim,order>::removeSimplex(*this);
            
            //! -2 remove this fomr children parentContainers
            for (int k=0;k<nFaces;++k)
            {
                this->child(k).removeFromParents(this);
            }
        }
        
        /**********************************************************************/
        ChildSimplexType& child(const int& n)
        {
            return *(this->operator[](n).get());
        }
        
        /**********************************************************************/
        const ChildSimplexType& child(const int& n) const
        {
            return *(this->operator[](n).get());
        }
        
        /**********************************************************************/
        const std::shared_ptr<ChildSimplexType>& child(const ChildIDType& xID) const
        {
            size_t n=nFaces;
            for (size_t k=0;k<nFaces;++k)
            {
                if(this->operator[](k)->xID==xID)
                {
                    n=k;
                    break;
                }
            }
            assert(n!=nFaces && "CHILD NOT FOUND");
            return this->operator[](n);
        }
        
        /**********************************************************************/
        std::vector<int> boundaryFaces() const
        {
            std::vector<int> temp;
            temp.reserve(nFaces);
            for (int n=0;n<nFaces;++n)
            {
                if(child(n).isBoundarySimplex())
                {
                    temp.push_back(n);
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        bool isBoundarySimplex() const
        {
            return BoundarySimplex<dim,dim-order>::isBoundarySimplex(*this);
        }
        
	};
    
    
}	// close namespace
#endif