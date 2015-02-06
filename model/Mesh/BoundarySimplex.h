/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BoundarySimplex_H_
#define model_BoundarySimplex_H_

namespace model {
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim,short int dmo>
    struct BoundarySimplex
    {
        /**********************************************************************/
        template<template <short int,short int> class SimplexChildType>
        static bool isBoundarySimplex(const SimplexChildType<dim,dim-dmo>& simplexChild)
        {
            bool temp(false);
            for (typename SimplexChildType<dim,dim-dmo>::ParentContainerType::const_iterator pIter=simplexChild.parentBegin();pIter!=simplexChild.parentEnd();++pIter)
            {
                temp=(*pIter)->isBoundarySimplex();
                if (temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        template<template <short int,short int> class SimplexChildType>
        static bool isRegionBoundarySimplex(const SimplexChildType<dim,dim-dmo>& simplexChild)
        {
            bool temp(false);
            for (typename SimplexChildType<dim,dim-dmo>::ParentContainerType::const_iterator pIter=simplexChild.parentBegin();pIter!=simplexChild.parentEnd();++pIter)
            {
                temp=(*pIter)->isRegionBoundarySimplex();
                if (temp)
                {
                    break;
                }
            }
            return temp;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim>
    struct BoundarySimplex<dim,1>
    {
        /**********************************************************************/
        template<template <short int,short int> class SimplexChildType>
        static bool isBoundarySimplex(const SimplexChildType<dim,dim-1>& simplexChild)
        {
            return simplexChild.parents().size()==1;
        }
        
        /**********************************************************************/
        template<template <short int,short int> class SimplexChildType>
        static bool isRegionBoundarySimplex(const SimplexChildType<dim,dim-1>& simplexChild)
        {
            bool temp(false);
            if(simplexChild.parents().size()==2)
            {
                if((*simplexChild.parents().begin())->region->regionID != (*simplexChild.parents().rbegin())->region->regionID)
                {
                    temp=true;
                }
                
            }
            
            return temp;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim>
    struct BoundarySimplex<dim,0>
    {
        
        /**********************************************************************/
        static bool isBoundarySimplex(const Simplex<dim,dim>& simplex)
        {
            bool temp(false);
            
            for (int n=0;n<simplex.nFaces;++n)
            {
                temp=simplex.child(n).isBoundarySimplex();
                if (temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        static bool isRegionBoundarySimplex(const Simplex<dim,dim>& simplex)
        {
            bool temp(false);
            
            for (int n=0;n<simplex.nFaces;++n)
            {
                temp=simplex.child(n).isRegionBoundarySimplex();
                if (temp)
                {
                    break;
                }
            }
            return temp;
        }
        
    };
    
}	// close namespace
#endif