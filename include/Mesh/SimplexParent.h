/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexParent_H_
#define model_SimplexParent_H_

#include <memory>
#include <map>
#include <SimplexTraits.h>

namespace model
{
    
    template<short int dim,short int order>
    struct SimplexParent : public SimplexTraits<dim,order>::BaseArrayType
    {
       
        typedef typename SimplexTraits<dim,order>::SimplexIDType SimplexIDType;
        typedef typename SimplexTraits<dim,order>::BaseArrayType BaseArrayType;
        typedef Simplex<dim,order-1> ChildSimplexType;
        typedef typename SimplexTraits<dim,order-1>::SimplexIDType ChildIDType;
        static constexpr int nFaces=SimplexTraits<dim,order>::nFaces;

        
        //! The barycentric-coordinate to position transformation matrix
        const Eigen::Matrix<double,dim+1,SimplexTraits<dim,order>::nVertices> b2p;

        
        /**********************************************************************/
        SimplexParent(SimplicialMesh<dim>* const mesh,
                      const SimplexIDType& xID) :
        /* init */ BaseArrayType(mesh->template observer<order>().faces(mesh,xID))
        /* init */,b2p(get_b2p(mesh,xID))
        {
            
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,SimplexTraits<dim,order>::nVertices>  vertexPositionMatrix() const
        {
            return b2p.template block<dim,SimplexTraits<dim,order>::nVertices>(0,0);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim+1,SimplexTraits<dim,order>::nVertices> get_b2p(SimplicialMesh<dim>* const mesh,
                                                                                const SimplexIDType& xID) const
        {
            Eigen::Matrix<double,dim+1,SimplexTraits<dim,order>::nVertices> temp(Eigen::Matrix<double,dim+1,SimplexTraits<dim,order>::nVertices>::Ones());
            //temp.template block<dim,SimplexTraits<dim,order>::nVertices>(0,0)=this->vertexPositionMatrix();
            for (int v=0;v<SimplexTraits<dim,order>::nVertices;++v)
            {
                const typename SimplexTraits<dim,0>::SimplexIDType vID(std::set<size_t>{xID[v]});
                temp.col(v).template segment<dim>(0)=mesh->template observer<0>().simplex(vID).P0;
            }
            return temp;
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,1> bary2pos(const Eigen::Matrix<double,SimplexTraits<dim,order>::nVertices,1>& bary) const
        {
            return (b2p*bary).template segment<dim>(0);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,1>  center() const
        {
            return bary2pos(Eigen::Matrix<double,SimplexTraits<dim,order>::nVertices,1>::Ones()/SimplexTraits<dim,order>::nVertices);
        }
        
        /**********************************************************************/
        BaseArrayType& children()
        {
            return *this;
        }
        
        /**********************************************************************/
        const BaseArrayType& children() const
        {
            return *this;
        }

        
        /**********************************************************************/
        ChildSimplexType& child(const int& n)
        {
            return *(children()[n].get());
        }
        
        /**********************************************************************/
        const ChildSimplexType& child(const int& n) const
        {
            return *(children()[n].get());
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
        
        
        
    };
    
}
#endif



//        /**********************************************************************/
//        Eigen::Matrix<double,dim,SimplexTraits<dim,order>::nVertices>  vertexPositionMatrix() const
//        {
//            Eigen::Matrix<double,dim,SimplexTraits<dim,order>::nVertices> temp;
//            for (int v=0;v<SimplexTraits<dim,order>::nVertices;++v)
//            {
//                const typename SimplexTraits<dim,0>::SimplexIDType vID(std::set<size_t>{xID[v]});
//                temp.col(v)=mesh->template observer<0>().simplex(vID).P0;
//            }
//            return temp;
//        }

//        /**********************************************************************/
//        Eigen::Matrix<double,dim,1>  center() const
//        {
//            Eigen::Matrix<double,dim,SimplexTraits<dim,order>::nVertices> P(vertexPositionMatrix());
//            Eigen::Matrix<double,dim,1> temp(Eigen::Matrix<double,dim,1>::Zero());
//            for (int v=0;v<SimplexTraits<dim,order>::nVertices;++v)
//            {
//                temp+=P.col(v);
//            }
//            return temp/SimplexTraits<dim,order>::nVertices;
//        }
