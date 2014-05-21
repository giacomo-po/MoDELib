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
#include <model/Network/Readers/VertexReader.h>

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
//            if (isFirstSimplex)
//            {
//                nodeReader.read(meshID,true);
//                isFirstSimplex=false;
//            }
            const typename VertexReader<'N',dim+1,double>::const_iterator nIter(nodeReader.find((this->xID)(0)));
            assert((nIter!=nodeReader.end()) && "MESH VERTEX NOT FOUND IN N/N_x.txt.");
            return nIter->second;
        }
        
//        static bool isFirstSimplex;
//        static VertexReader<'N',dim+1,double> nodeReader;
    
    public:
        
//        static int meshID;
//        static bool isFirstSimplex;
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

//    template<short int dim>
//    int Simplex<dim,0>::meshID=0;
//
//    template<short int dim>
//    bool Simplex<dim,0>::isFirstSimplex=true;

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
        
//        bool beingSearched;
        
        /**********************************************************************/
        Eigen::Matrix<double,dim+1,dim+1> get_b2p() const
        {
            Eigen::Matrix<double,dim+1,dim+1> temp(Eigen::Matrix<double,dim+1,dim+1>::Ones());
            temp.template block<dim,dim+1>(0,0)=this->vertexPositionMatrix();
            return temp;
        }
        
    public:
        enum{order=dim};
        enum{nVertices=SimplexTraits<dim,order>::nVertices};
        enum{nFaces=SimplexTraits<dim,order>::nFaces};
        typedef typename SimplexTraits<dim,order>::SimplexIDType SimplexIDType;
        
        typedef Simplex<dim,order-1> ChildSimplexType;
        typedef typename SimplexTraits<dim,order-1>::SimplexIDType ChildIDType;
        
        const Eigen::Matrix<double,dim+1,dim+1> b2p;
        const Eigen::Matrix<double,dim+1,dim+1> p2b;
        
		/**********************************************************************/
        Simplex(const SimplexIDType& vIN) :
        /* init list */ SimplexBase<dim,order>(vIN),
        /* init list */ BaseArrayType(SimplexObserver<dim,dim>::faces(vIN)),
//        /* init list */ beingSearched(false),
        /* init list */ b2p(get_b2p()),
        /* init list */ p2b(b2p.inverse())
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
        
        
        /**********************************************************************/
        Eigen::Matrix<double,dim+1,1> pos2bary(const Eigen::Matrix<double,dim,1>& P) const
        {
            return p2b*(Eigen::Matrix<double,dim+1,1>()<<P,1.0).finished();
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,1> bary2pos(const Eigen::Matrix<double,dim+1,1>& bary) const
        {
            return (b2p*bary).template segment<dim>(0);
        }
        
        /**********************************************************************/
        std::pair<bool,const Simplex<dim,dim>*> search(const Eigen::Matrix<double,dim,1>& P, std::set<int>& searchSet) const
        {
            std::pair<bool,const Simplex<dim,dim>*> temp(false,NULL);
            if(searchSet.find(this->sID)==searchSet.end())
            {// this simplex has not been searched yet
                searchSet.insert(this->sID);
                const Eigen::Matrix<double,dim+1,1> bary(pos2bary(P));
#ifdef _MODEL_BENCH_SEARCH_
                searchFile<<bary2pos(Eigen::Matrix<double,dim+1,1>::Ones()/(dim+1)).transpose()<<"\n"; // REMOVE THIS
#endif
                int kMin;
                const double baryMin(bary.minCoeff(&kMin));
                int kMax;
                const double baryMax(bary.maxCoeff(&kMax));
                if (baryMin>=0.0 && baryMax<=1.0)
                {
                    temp = make_pair(true,this);
                }
                else
                {
                    for(typename Simplex<dim,dim-1>::ParentContainerType::const_iterator pIter=this->child(kMin).parentBegin();
                        /*                                                            */ pIter!=this->child(kMin).parentEnd();++pIter)
                    {
                        temp=(*pIter)->search(P,searchSet);
                        if (temp.first)
                        {
                            break;
                        }
                    }
                }
            }
            return temp;
        }

        
	};
    
    
}	// close namespace
#endif