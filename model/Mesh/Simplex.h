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
#include <model/Mesh/SimplexReader.h>
#include <model/Mesh/SimplexTraits.h>
#include <model/Mesh/SimplexObserver.h>
#include <model/Mesh/SimplexBase.h>
#include <model/Mesh/SimplexChild.h>
#include <model/Mesh/MeshRegion.h>
#include <model/Network/Readers/VertexReader.h>
#include <model/FEM/BarycentricTraits.h>
#include <model/Mesh/SimplexVolume.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim>
    class Simplex<dim,0> :
//    /* inheritance      */ public SimplexBase<dim,0>,
    /* inheritance      */ public SimplexChild <dim,0>
    {
        
    public:
        
        enum{order=0};
        enum{nVertices=SimplexTraits<dim,order>::nVertices};
        enum{nFaces=SimplexTraits<dim,order>::nFaces};
        typedef typename SimplexTraits<dim,order>::SimplexIDType SimplexIDType;
//        typedef SimplexChild<dim,order> SimplexChildType;

        
        const Eigen::Matrix<double,dim,1> P0;
        
        /**********************************************************************/
        Simplex(SimplicialMesh<dim>* const m,
                const SimplexIDType& vIN) :
//        /* init list */ SimplexBase<dim,order>(vIN),
        /* init list */ SimplexChild<dim,order>(m,vIN),
        //        /* init list */ P0(get_P0())
        /* init list */ P0(SimplexReader<dim>::get_P0(this->xID))
        {/*!@param[in] vIN the (possibly unsorted) ID of this Simplex
          *
          * Constructur performs the following operations:
          */
            
//                        //std::cout<<"Creating Simplex<3,"<<order<<std::endl;
            
            //! -1 Adds this to the SimplexObserver
            this->observer().insertSimplex(*this);
            
//            //std::cout<<"done"<<std::endl;

        }
        
        /**********************************************************************/
        ~Simplex()
        {/*!
          */
            //std::cout<<"Destroying Simplex<3,"<<order<<"> ..."<<this->xID<<std::flush;

            this->observer().removeSimplex(*this);
            //std::cout<<" done observer"<<std::endl;

        }
        
        /**********************************************************************/
        void updateBndNormal()
        {
            // Update normal of this SImplex
            this->outN=BoundarySimplex<dim,dim-order>::outNormal(*this);
        }
        
//        /**********************************************************************/
//        Eigen::Matrix<double,dim,1> outNormal() const
//        {
//            return BoundarySimplex<dim,dim-order>::outNormal(*this);
//        }
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim, short int order>
    class Simplex :
//    /* inheritance */ public SimplexBase<dim,order>,
    /* inheritance */ public SimplexChild <dim,order>,
    /* inheritance */ public SimplexTraits<dim,order>::BaseArrayType
    {
        
        typedef typename SimplexTraits<dim,order>::BaseArrayType BaseArrayType;
        //typedef SimplexChild<dim,order> SimplexChildType;

        
    public:
        enum{nVertices=SimplexTraits<dim,order>::nVertices};
        enum{nFaces=SimplexTraits<dim,order>::nFaces};
        typedef typename SimplexTraits<dim,order>::SimplexIDType SimplexIDType;
        
        typedef Simplex<dim,order-1> ChildSimplexType;
        typedef typename SimplexTraits<dim,order-1>::SimplexIDType ChildIDType;
        
        double vol0; // SHOULD BE CONST, SEE BELOW
        
        /**********************************************************************/
        Simplex(SimplicialMesh<dim>* const m,
                const SimplexIDType& vIN) :
//        /* init list */ SimplexBase<dim,order>(vIN),
        /* init list */ SimplexChild<dim,order>(m,vIN),
        /* init list */ BaseArrayType(this->observer().faces(vIN))
        //        /* init */ vol0(SimplexVolume<dim,order>::volume(this->vertexPositionMatrix())) // THIS GIVES SEGMENTATION FAULT, WHY?
        {/*!@param[in] vIN the (possibly unsorted) ID of this
          *
          * Constructur performs the following operations:
          */
            
//                        //std::cout<<"Creating Simplex<3,"<<order<<std::endl;
            
            //! -1 inserts *this into the SimplexObserver
            this->observer().insertSimplex(*this);
            
            //! -2 inserts this into children Simplices
            for (int k=0;k<nFaces;++k)
            {
                child(k).addToParents(this);
            }
            
            vol0=SimplexVolume<dim,order>::volume(this->vertexPositionMatrix());
            
//                        //std::cout<<"done"<<std::endl;
        }
        
        /**********************************************************************/
        ~Simplex()
        {/*! Destructor performs the following operations:
          */
            //std::cout<<"Destroying Simplex<3,"<<order<<"> ..."<<this->xID<<std::flush;
            
            //! -1 removes this in SimplexObserver
            this->observer().removeSimplex(*this);
            //std::cout<<" done observer"<<std::flush;

            
            //! -2 remove this fomr children parentContainers
            for (int k=0;k<nFaces;++k)
            {
                child(k).removeFromParents(this);
            }
            
            //std::cout<<" done parents"<<std::endl;

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
        size_t childOrder(const ChildIDType& childID) const
        {
            return SimplexTraits<dim,dim>::faceOrder(this->xID,childID);
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
        void updateBndNormal()
        {
            // Update normal of this SImplex
            this->outN=BoundarySimplex<dim,dim-order>::outNormal(*this);
            
            // Update normals of children
            for (int k=0;k<nFaces;++k)
            {
                child(k).updateBndNormal();
            }
        }
        
//        /**********************************************************************/
//        Eigen::Matrix<double,dim,1> outNormal() const
//        {
//            return BoundarySimplex<dim,dim-order>::outNormal(*this);
//        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim>
    class Simplex<dim,dim> :
    /* inheritance */ public SimplexBase<dim,dim>,
    /* inheritance */ public SimplexTraits<dim,dim>::BaseArrayType
    {
        
        typedef typename SimplexTraits<dim,dim>::BaseArrayType BaseArrayType;
        
        /**********************************************************************/
        Eigen::Matrix<double,dim+1,dim+1> get_b2p() const
        {
            Eigen::Matrix<double,dim+1,dim+1> temp(Eigen::Matrix<double,dim+1,dim+1>::Ones());
            temp.template block<dim,dim+1>(0,0)=this->vertexPositionMatrix();
            return temp;
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,dim+1> get_nda() const
        {
            Eigen::Matrix<double,dim,dim+1> vP(this->vertexPositionMatrix());
            Eigen::Matrix<double,dim,dim> F(vP.template block<dim,dim>(0,0));
            F.colwise() -= vP.col(dim);
            const double jFabs(std::fabs(F.determinant()));
            if(jFabs<DBL_EPSILON)
            {
                model::cout<<this->xID<<", volume="<<jFabs<<std::endl;
                model::cout<<F<<std::endl;
                assert(0 && "SIMPLEX HAS ZERO VOLUME");
            }
            return jFabs*F.inverse().transpose()*BarycentricTraits<dim>::NdA;
        }
        
        /**********************************************************************/
        void updateBndNormal()
        {
            for (int k=0;k<nFaces;++k)
            {
                child(k).updateBndNormal();
            }
        }
        
    public:
        enum{order=dim};
        enum{nVertices=SimplexTraits<dim,order>::nVertices};
        enum{nFaces=SimplexTraits<dim,order>::nFaces};
        typedef typename SimplexTraits<dim,order>::SimplexIDType SimplexIDType;
        
        typedef Simplex<dim,order-1> ChildSimplexType;
        typedef typename SimplexTraits<dim,order-1>::SimplexIDType ChildIDType;
        
        typedef MeshRegion<Simplex<dim,dim> > MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        
        //! Shared pointer to MeshRegioin containing *this
        const std::shared_ptr<MeshRegionType> region;
        
        //! The barycentric-coordinate to position transformation matrix
        const Eigen::Matrix<double,dim+1,dim+1> b2p;
        
        //! The position to barycentric-coordinate transformation matrix
        const Eigen::Matrix<double,dim+1,dim+1> p2b;
        
        //! The column matrix of face normals
        const Eigen::Matrix<double,dim,dim+1> nda;
        
        const double vol0;
        
        /**********************************************************************/
        Simplex(SimplicialMesh<dim>* const m,
                const SimplexIDType& vIN, const int regionID=0) :
        /* init base */ SimplexBase<dim,order>(m,vIN),
        /* init base */ BaseArrayType(this->observer().faces(vIN)),
        /* init list */ region(m->getSharedRegion(regionID)),
        /* init base */ b2p(get_b2p()),
        /* init list */ p2b(b2p.fullPivLu().solve(Eigen::Matrix<double,dim+1,dim+1>::Identity())),
        /* init list */ nda(get_nda()),
        /* init */ vol0(SimplexVolume<dim,order>::volume(this->vertexPositionMatrix()))
        {/*!
          */
            
//            //std::cout<<"Creating Simplex<3,"<<order<<std::endl;
            
            this->observer().insertSimplex(*this);
            
            
            for (int k=0;k<nFaces;++k)
            {
                child(k).addToParents(this);
            }
            
            region->simplices().emplace(this);
            updateBndNormal();
            
//                        //std::cout<<"done"<<std::endl;
        }
        
        /**********************************************************************/
        ~Simplex()
        {/*! Destructor performs the following operations:
          */
            
            //std::cout<<"Destroying Simplex<3,"<<order<<"> ..."<<this->xID<<std::flush;

            //! -1 removes this in SimplexObserver
//            SimplexObserver<dim,order>::removeSimplex(*this);
            this->observer().removeSimplex(*this);
            //std::cout<<" done observer"<<std::flush;

            
            //! -2 remove this fomr children parentContainers
            for (int k=0;k<nFaces;++k)
            {
                child(k).removeFromParents(this);
            }
            
            //std::cout<<" done parents"<<std::flush;

            
            region->simplices().erase(this);
            updateBndNormal();
            
            //std::cout<<" done region"<<std::endl;

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
        size_t childOrder(const ChildIDType& childID) const
        {
            return SimplexTraits<dim,dim>::faceOrder(this->xID,childID);
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
        bool isRegionBoundarySimplex() const
        {
            return BoundarySimplex<dim,dim-order>::isRegionBoundarySimplex(*this);
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
        void convexDelaunaynSearch(const bool& searchAllRegions,
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
                //                double tol=0.0;
                //                if(child(kMin).isBoundarySimplex())
                //                {
                //                    tol=-FLT_EPSILON;
                //                }
                
#ifdef _MODEL_BENCH_BARYSEARCH_
                model::cout<<"Searching "<<this->xID<<std::endl;
                model::cout<<"bary= "<<pos2bary(P)<<std::endl;
                searchFile<<bary2pos(Eigen::Matrix<double,dim+1,1>::Ones()/(dim+1)).transpose()<<" "
                /*      */<<this->xID<<"\n";
#endif
                
                
                
                
                if (baryMin>=-FLT_EPSILON)
                {
                    lastSearched.first=true;
                }
                else
                {
                    
                    for(auto& pParent : child(kMin).parents())
                    {
                        if(pParent->region->regionID==region->regionID || searchAllRegions)
                        {
                            pParent->convexDelaunaynSearch(searchAllRegions,P,lastSearched,searchSet);
                            if (lastSearched.first)
                            {
                                break;
                            }
                        }
                    }
                }
                
            }
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim+1,1> faceLineIntersection(const Eigen::Matrix<double,dim+1,1>& bary0,
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
                model::cout<<"Parallel faceLineIntersection"<<std::endl;
            }
            return temp;
            
        }
        
        /**********************************************************************/
        std::set<int> regionIDs() const
        {
            std::set<int> temp;
            temp.insert(region->regionID);
            return temp;
        }
        
        
    };
    
    
}	// close namespace
#endif



//        /**********************************************************************/
//        void convexDelaunaynSearch(const Eigen::Matrix<double,dim,1>& P,
//                                   std::pair<bool,const Simplex<dim,dim>*>& lastSearched,
//                                   std::set<int>& searchSet) const // TO DO: searchSet is not necessary, because baryMin changes sign in next Simplex
//        {
//            if(searchSet.find(this->sID)==searchSet.end())
//            {// this simplex has not been searched yet
//                searchSet.insert(this->sID);
//                lastSearched.second=this;
//#ifdef _MODEL_BENCH_BARYSEARCH_
//                const Eigen::Matrix<double,dim+1,1> bary(pos2bary(P));
//                searchFile<<bary2pos(Eigen::Matrix<double,dim+1,1>::Ones()/(dim+1)).transpose()<<"\n";
//#endif
//                int kMin;
//                if (pos2bary(P).minCoeff(&kMin)>=0.0)
//                {
//                    lastSearched.first=true;
//                }
//                else
//                {
////                    for(typename Simplex<dim,dim-1>::ParentContainerType::const_iterator pIter=child(kMin).parentBegin();
////                        /*                                                            */ pIter!=child(kMin).parentEnd();++pIter)
////                    {
//                    for(auto& pParent : child(kMin).parents())
//                    {
////                        (*pIter)->convexDelaunaynSearch(P,lastSearched,searchSet);
//                    pParent->convexDelaunaynSearch(P,lastSearched,searchSet);
//                        if (lastSearched.first)
//                        {
//                            break;
//                        }
//                    }
//                }
//            }
//        }
