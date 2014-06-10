/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NodeList_H_
#define model_NodeList_H_

#include <list>
//#include <Eigen/Dense>

namespace model
{
    template <typename FiniteElementType>
    struct NodeList : public std::list<const typename FiniteElementType::NodeType*>
    {
        
        typedef NodeList<FiniteElementType> NodeListType;
        
        const FiniteElementType& fe;
        
        /**************************************/
        NodeList(const FiniteElementType& fe_in) :
        fe(fe_in)
        {
            
        }

        /**************************************/
        NodeList(NodeListType&&) = default; // force a move assignment

        
        /**************************************/
        NodeListType& operator=(NodeListType&&) = default; // force a move assignment

    private:
        
        
        /**************************************/
        NodeList(const NodeListType&) = default; //Copy constructor is private to avoid copies.
//        {/*!
//          */
//        }
        
        /**************************************/
        NodeListType& operator=(const NodeListType&) = default; // Assignment operator is private to avoid copies.
//        {/*! Assignment operator is private to avoid copies.
//          */
//        }

        
    };
    
    
}	// close namespace
#endif

