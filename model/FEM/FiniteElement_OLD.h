/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FiniteElement_H_
#define model_FiniteElement_H_

#include <deque>
#include <map>

#include <Eigen/Dense>

#include <model/Mesh/SimplicialMesh.h>
#include <model/Utilities/TerminalColors.h>
#include <model/Utilities/CompareVectorsByComponent.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/FEM/Elements/LagrangeElement.h>
#include <model/FEM/TrialFunction.h>

#include <model/FEM/WeakForms/LinearWeakForm.h>
#include <model/FEM/WeakForms/BilinearWeakForm.h>
#include <model/FEM/Domains/IntegrationDomain.h>
#include <model/FEM/Domains/ExternalBoundary.h>

//#include <model/FEM/WeakForms/LinearDomainAssemble.h>
//#include <model/FEM/WeakForms/LinearBoundaryAssemble.h>


namespace model {
    
    template<typename _ElementType>
    class FiniteElement :
    /* inherit         */ public std::deque<_ElementType>,
    /* inherit         */ public std::deque<typename _ElementType::NodeType>,
    /* inherit         */ public std::map<Eigen::Matrix<double,_ElementType::dim,1>, // key
    /*                                  */ const typename _ElementType::NodeType* const, // value
    /*                                  */ CompareVectorsByComponent<double,_ElementType::dim,float> // key compare
    /*                                  */ >,
    /*                 */ private std::map<size_t,std::list<const typename _ElementType::NodeType* > >
    {
        
    public:
        typedef std::list<const typename _ElementType::NodeType* > NodeListType;
        typedef std::map<size_t,NodeListType> NodeMapType;
        
    private:
        size_t incrementalNodeListID;
        
        
        
        Eigen::Matrix<double,_ElementType::dim,1> _xMin;
        Eigen::Matrix<double,_ElementType::dim,1> _xMax;
        
        
    public:
        
        typedef _ElementType ElementType;
        typedef FiniteElement<ElementType> FiniteElementType;
        typedef typename ElementType::NodeType NodeType;
        
        constexpr static int dim=ElementType::dim;
        constexpr static int nodesPerElement=ElementType::nodesPerElement;
        
        
        typedef SimplicialMesh<dim> MeshType;
        typedef std::deque<ElementType> ElementContainerType;
        typedef std::deque<typename _ElementType::NodeType> NodeContainerType;
        typedef std::map<Eigen::Matrix<double,_ElementType::dim,1>, // key
        /*                                  */ const typename _ElementType::NodeType* const, // value
        /*                                  */ CompareVectorsByComponent<double,_ElementType::dim,float> // key compare
        /*                                  */ > NodeFinderType;
        
        
        const MeshType& mesh;
        
        
        /**********************************************************************/
        FiniteElement(const SimplicialMesh<dim>& m) :
        /* init list */ incrementalNodeListID(0),
        /* init list */ _xMin(Eigen::Matrix<double,ElementType::dim,1>::Constant( DBL_MAX)),
        /* init list */ _xMax(Eigen::Matrix<double,ElementType::dim,1>::Constant(-DBL_MAX)),
        /* init list */ mesh(m)
        {/*!@param[in] s A const reference to a Simplex<dim,dim>
          */
            
            std::cout<<greenColor<<"Creating FiniteElement:\n"<<defaultColor<<std::flush;
            
            // THIS IS NECESSARY TO AVOID "STATIC INITIALIZATION FIASCO"
            std::cout<<"Element barycentric coordinates:\n"<<ElementType::baryNodalCoordinates<<std::endl;
            
            
            // Insert elements
            for (typename SimplicialMesh<dim>::const_iterator eIter=mesh.begin();eIter!=mesh.end();++eIter)
            {
                ElementContainerType::emplace_back(eIter->second,*this,*this);
            }
            
            // Compute _xMin and _xMax
            for (int n=0;n<nodeSize();++n)
            {
                for(int d=0;d<dim;++d)
                {
                    if (node(n).p0(d)<_xMin(d))
                    {
                        _xMin(d)=node(n).p0(d);
                    }
                    if (node(n).p0(d)>_xMax(d))
                    {
                        _xMax(d)=node(n).p0(d);
                    }
                }
                
            }
            
            std::cout<<"   # elements: "<<elementSize()    <<"\n";
            std::cout<<"   # nodes: "   <<nodeSize()       <<"\n";
            std::cout<<"   # xMin: "    <<_xMin.transpose()<<"\n";
            std::cout<<"   # xMax: "    <<_xMax.transpose()<<std::endl;
            
        }
        
        /**********************************************************************/
        template <int nComponents>
        TrialFunction<nComponents,FiniteElementType> trial() const
        {
            return TrialFunction<nComponents,FiniteElementType>(*this);
        }
        
        /**********************************************************************/
        typename ElementContainerType::const_iterator elementBegin() const
        {
            return ElementContainerType::begin();
        }
        
        /**********************************************************************/
        typename ElementContainerType::const_iterator elementEnd() const
        {
            return ElementContainerType::end();
        }
        
        /**********************************************************************/
        size_t elementSize() const
        {/*!\returns the number of elements in the FiniteElement
          */
            return ElementContainerType::size();
        }
        
        /**********************************************************************/
        const ElementType& element(const size_t& n) const
        {/*!@param[in] n the node ID
          * \returns a const reference to the n-th node in the FiniteElement
          */
            return ElementContainerType::operator[](n);
        }
        
        /**********************************************************************/
        ElementType& element(const size_t& n)
        {/*!@param[in] n the node ID
          * \returns a reference to the n-th node in the FiniteElement
          */
            return ElementContainerType::operator[](n);
        }
        
        /**********************************************************************/
        typename NodeContainerType::const_iterator nodeBegin() const
        {
            return NodeContainerType::begin();
        }
        
        /**********************************************************************/
        typename NodeContainerType::const_iterator nodeEnd() const
        {
            return NodeContainerType::end();
        }
        
        /**********************************************************************/
        size_t nodeSize() const
        {
            return NodeContainerType::size();
        }
        
        /**********************************************************************/
        const NodeType& node(const size_t& n) const
        {/*!@param[in] n the node ID
          * \returns a reference to the n-th node in the FiniteElement
          */
            return NodeContainerType::operator[](n);
        }
        
        /**********************************************************************/
        template <typename BndType, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        IntegrationDomain<dim,1,qOrder,QuadratureRule> boundary() const
        {
            return BndType::template domain<FiniteElementType,qOrder,QuadratureRule>(*this);
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& xMin() const
        {
            return _xMin;
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& xMax() const
        {
            return _xMax;
        }
        
        
        /**********************************************************************/
        template <typename ListSelectorType>
        size_t createNodeList()
        {
            std::cout<<"Adding node list... "<<std::flush;
            const ListSelectorType listSelector(*this);
            const std::pair<typename NodeMapType::iterator,bool> ip(NodeMapType::emplace(incrementalNodeListID,NodeListType()));
            assert(ip.second && "NODE LIST COULD NOT BE INSERTED");
            for(typename NodeContainerType::const_iterator nIter=nodeBegin(); nIter!=nodeEnd();++nIter)
            {
                if(listSelector(*nIter))
                {
                    ip.first->second.emplace_back(&(*nIter));
                }
            }
            std::cout<<"("<<ip.first->second.size()<<" nodes)."<<std::endl;
            ++incrementalNodeListID; // increment key for next call
            return ip.first->first;
        }
        
        /**********************************************************************/
        const NodeListType& nodeList(const size_t& n) const
        {
            return NodeMapType::at(n);
        }
        
    };
    
}	// close namespace
#endif


//        /**********************************************************************/
//        template <typename BndType, int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
//        ExternalBoundary<BndType,dim,qOrder,QuadratureRule> externalBoundary() const
//        {
//            return ExternalBoundary<BndType,dim,qOrder,QuadratureRule>(*this);
//        }