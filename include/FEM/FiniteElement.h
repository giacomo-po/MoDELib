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
#include <list>
#include <stdexcept>      // std::out_of_range
#include <type_traits> // std::is_same


#include <Eigen/Dense>

#include <TrialFunctionTraits.h>
//#include <AreSameType.h>
#include <SimplicialMesh.h>
#include <TerminalColors.h>
#include <CompareVectorsByComponent.h>
#include <SimplicialMesh.h>
#include <LagrangeElement.h>
#include <DiscontinuousLagrangeElement.h>
#include <TrialFunction.h>
#include <LinearWeakForm.h>
#include <BilinearWeakForm.h>
#include <IntegrationDomain.h>
#include <EntireDomain.h>
#include <ExternalBoundary.h>
//#include <NodeList.h>



namespace model
{

    template<typename _ElementType>
    class FiniteElement :
//    /* inherits        */ public std::map<Eigen::Matrix<size_t,_ElementType::dim+1,1>, // key
//    /*                                  */ _ElementType, // value
//    /*                                  */ CompareVectorsByComponent<size_t,_ElementType::dim+1> // key compare
//    /*                                  */ >, // element container
    /* inherits        */ public std::map<typename SimplexTraits<_ElementType::dim,_ElementType::dim>::SimplexIDType, // key
    /*                                  */ _ElementType>, // value
//    /*                                  */ CompareVectorsByComponent<size_t,_ElementType::dim+1> // key compare
//    /*                                  */ >, // element container
    /* inherits        */ public std::deque<typename _ElementType::NodeType>, // node container
    /* inherits        */ public std::map<Eigen::Matrix<double,_ElementType::dim,1>, // key
    /*                                  */ typename _ElementType::NodeType* const, // value
    /*                                  */ CompareVectorsByComponent<double,_ElementType::dim,float> // key compare
    /*                                  */ >, // nodefinder
    /* inherits        */ private std::map<size_t,std::deque<const typename _ElementType::NodeType*>>,// node-list container
    /* inherits        */ private std::map<size_t,const typename _ElementType::NodeType* const>
    {
        
    private:
        Eigen::Matrix<double,_ElementType::dim,1> _xMin;
        Eigen::Matrix<double,_ElementType::dim,1> _xMax;
        Eigen::Matrix<double,_ElementType::dim,1> _xMean;

        size_t nodeListID;

        
    public:
        
        typedef _ElementType ElementType;
        typedef FiniteElement<ElementType> FiniteElementType;
        typedef typename ElementType::NodeType NodeType;
//        typedef std::list<const NodeType*> NodeListType;
//        typedef std::map<size_t,NodeListType> NodeMapType;
        typedef std::deque<const NodeType*> NodeListType;
        typedef std::map<size_t,NodeListType> NodeListContainerType;

        constexpr static int dim=ElementType::dim;
        constexpr static int nodesPerElement=ElementType::nodesPerElement;
        typedef SimplicialMesh<dim> MeshType;
//        typedef std::map<Eigen::Matrix<size_t,_ElementType::dim+1,1>, // key
//        /*                                  */ _ElementType, // value
//        /*                                  */ CompareVectorsByComponent<size_t,_ElementType::dim+1> // key compare
//        /*                                  */ > ElementContainerType;
        typedef std::map<typename SimplexTraits<_ElementType::dim,_ElementType::dim>::SimplexIDType, // key
        /*                                  */ _ElementType> ElementContainerType;
        typedef std::deque<typename _ElementType::NodeType> NodeContainerType;
        typedef std::map<Eigen::Matrix<double,_ElementType::dim,1>, // key
        /*                                  */ typename _ElementType::NodeType* const, // value
        /*                                  */ CompareVectorsByComponent<double,_ElementType::dim,float> // key compare
        /*                                  */ > NodeFinderType;
        
        const MeshType& mesh;
        
        /**********************************************************************/
        FiniteElement(const SimplicialMesh<dim>& m) :
        /* init list */ _xMin(Eigen::Matrix<double,dim,1>::Zero()),
        /* init list */ _xMax(Eigen::Matrix<double,dim,1>::Zero()),
        /* init list */ nodeListID(0),
        /* init list */ mesh(m)
        {/*!@param[in] s A const reference to a SimplicialMesh on which *this
          * FiniteElement is constructed.
          */
            
             std::cout<<greenBoldColor<<"Creating FiniteElement:\n"<<defaultColor<<std::flush;
            
            // THIS IS NECESSARY TO AVOID "STATIC INITIALIZATION FIASCO"
             std::cout<<"Element barycentric coordinates:\n"<<ElementType::baryNodalCoordinates<<std::endl;
            
            // Insert elements
            for (const auto& simpl : mesh.simplices())
            {
                auto temp=ElementContainerType::emplace(simpl.first,ElementType(simpl.second,*this,*this,mesh2femIDmap()));
                assert(temp.second && "UNABLE TO INSERT ELEMENT IN ELEMENT CONTAINER.");
                
                // Add element pointer to each of its nodes
                for(int n=0;n<ElementType::nodesPerElement;++n)
                {
                    auto temp1=temp.first->second.node(n).emplace(&temp.first->second);
                    if(!temp1.second)
                    {
                        throw std::runtime_error("UNABLE TO INSERT ELEMENT IN NODE.");
                    }
//                    assert(temp1.second && "UNABLE TO INSERT ELEMENT IN NODE.");
                }
            }
            
            if(std::is_same<ElementType,LagrangeElement<ElementType::dim,ElementType::order>>::value)
            {
                if(mesh2femIDmap().size()!=mesh.template observer<0>().size())
                {
                    throw std::runtime_error("mesh2femIDmap has wrong size.");
                }
//                assert((mesh2femIDmap().size()==mesh.template observer<0>().size()) && "mesh2femIDmap has wrong size.");
            
            }
            
            // Check that node[k].gID==k;
            for(size_t n=0;n<nodes().size();++n)
            {
                assert(node(n).gID==n);
            }
            
            // Compute _xMin and _xMax
            
            if(nodeSize())
            {
                _xMin=node(0).P0;
                _xMax=node(0).P0;
                
                for (size_t n=0;n<nodeSize();++n)
                {
                    for(int d=0;d<dim;++d)
                    {
                        if (node(n).P0(d)<_xMin(d))
                        {
                            _xMin(d)=node(n).P0(d);
                        }
                        if (node(n).P0(d)>_xMax(d))
                        {
                            _xMax(d)=node(n).P0(d);
                        }
                    }
                }
                _xMean=0.5*(_xMin+_xMax);

            }
            
             std::cout<<"   # elements: "<<elementSize()    <<"\n";
             std::cout<<"   # nodes: "   <<nodeSize()       <<"\n";
             std::cout<<"   xMin= "    <<_xMin.transpose()<<"\n";
             std::cout<<"   xMax= "    <<_xMax.transpose()<<std::endl;
        }
        
        /**********************************************************************/
        template <char name,int nComponents>
        TrialFunction<name,nComponents,FiniteElementType> trial()  // made non-const only to allow fe.createNodeList
        {
            return TrialFunction<name,nComponents,FiniteElementType>(*this);
        }
        
        /**********************************************************************/
        const ElementContainerType& elements() const
        {
            return *this;
        }
        
        /**********************************************************************/
        size_t elementSize() const
        {/*!\returns the number of elements in the FiniteElement
          */
            return ElementContainerType::size();
        }
        
        /**********************************************************************/
        const NodeContainerType& nodes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        size_t nodeSize() const
        {
            return NodeContainerType::size();
        }
        
        /**********************************************************************/
        const NodeType& node(const size_t& n) const
        {/*!@param[in] n the n-th node stored in *this FiniteElement
          * \returns a reference to the n-th node in *this FiniteElement
          */
            return NodeContainerType::operator[](n);
        }
        
        /**********************************************************************/
        const NodeFinderType& nodeFinder() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& xMin() const
        {
            return _xMin;
        }
        
        /**********************************************************************/
        const double& xMin(const int& k) const
        {
            return _xMin(k);
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& xMax() const
        {
            return _xMax;
        }
        
        /**********************************************************************/
        const double& xMax(const int& k) const
        {
            return _xMax(k);
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& xMean() const
        {
            return _xMean;
        }
        
        /**********************************************************************/
        const double& xMean(const int& k) const
        {
            return _xMean(k);
        }
        
        /**********************************************************************/
        template <typename BndType, int qOrder, template <short unsigned int, size_t> class QuadratureRule>
        IntegrationDomain<FiniteElementType,1,qOrder,QuadratureRule> boundary() const
        {
            return BndType::template boundary<FiniteElementType,qOrder,QuadratureRule>(*this);
        }
        
        /**********************************************************************/
        template <typename DomainType, int qOrder, template <short unsigned int, size_t> class QuadratureRule>
        IntegrationDomain<FiniteElementType,0,qOrder,QuadratureRule> domain() const
        {
            return DomainType::template domain<FiniteElementType,qOrder,QuadratureRule>(*this);
        }
        
        /**********************************************************************/
        const std::map<size_t,const typename _ElementType::NodeType* const>& mesh2femIDmap() const
        {
            return *this;
        }

        /**********************************************************************/
        std::map<size_t,const typename _ElementType::NodeType* const>& mesh2femIDmap()
        {
            return *this;
        }
        
        /**********************************************************************/
        template <typename NodeSelectorType, typename... NodeSelectorArgs>
        size_t createNodeList(const NodeSelectorArgs&... args)
        {
            //const size_t nodeListID_old(nodeListID);
            auto pair=NodeListContainerType::emplace(nodeListID,NodeListType());
            assert(pair.second);
            const NodeSelectorType nodeSelector(*this,args...);
            for(auto& node : nodes())
            {
                if(nodeSelector(node))
                {
                    pair.first->second.emplace_back(&node);
                }
            }
            nodeListID++;
            return pair.first->first;
        }
        
        /**********************************************************************/
        void  clearNodeLists()
        {
            return NodeListContainerType::clear();
        }
        
        /**********************************************************************/
        const NodeListType& nodeList(const size_t& k) const
        {
            try
            {
                return NodeListContainerType::at(k);
            }
            catch (const std::out_of_range& oor)
            {
                std::cerr << "FiniteElement::nodeList(), out_of_range error: " << oor.what() << "\n";
                return NodeListContainerType::at(k);
            }
        }
        
        /**********************************************************************/
        std::pair<bool,const ElementType*> search(const Eigen::Matrix<double,dim,1>& P) const
        {/*!@param[in] P position to search for
          *\returns a pair, where:
          * -pair.first is a boolean indicating whether the
          * search succesfully found a Simplex<dim,dim> which includes P.
          * -pair.second is a pointer to the last Simplex<dim,dim> searched.
          *
          * By default the search starts at this->begin()->second
          */
            return searchWithGuess(P,&(elements().begin()->second.simplex));
        }
        
        /**********************************************************************/
        std::pair<bool,const ElementType*> searchWithGuess(const Eigen::Matrix<double,dim,1>& P, const Simplex<dim,dim>* const guess) const
        {/*!@param[in] P position to search for
          * @param[in] guess Simplex* where the search starts
          *\returns a pair, where:
          * -pair.first is a boolean indicating whether the
          * search succesfully found a Simplex<dim,dim> which includes P.
          * -pair.second is a pointer to the last Simplex<dim,dim> searched.
          */
            const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(P,guess));
            const typename ElementContainerType::const_iterator eIter(ElementContainerType::find(temp.second->xID));
//            const typename ElementContainerType::const_iterator eIter(ElementContainerType::find(Eigen::Map<const Eigen::Matrix<size_t,_ElementType::dim+1,1>>(temp.second->xID.data())));
            assert(eIter!=elements().end() && "ELEMENT NOT FOUND");
            return std::pair<bool,const ElementType*>(temp.first,&(eIter->second));
        }
        
    };
    
}	// close namespace
#endif
