/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BoundarySimplex_H_
#define model_BoundarySimplex_H_

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim,short int dmo>
    struct BoundarySimplex
    {
        /**********************************************************************/
        static bool isBoundarySimplex(const SimplexChild<dim,dim-dmo>& simplexChild)
        {
            bool temp(false);
            for (auto& pParent : simplexChild.parents())
            {
                temp=pParent->isBoundarySimplex();
                if (temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        static bool isRegionBoundarySimplex(const SimplexChild<dim,dim-dmo>& simplexChild)
        {
            bool temp(false);
            for (auto& pParent : simplexChild.parents())
            {
                temp=pParent->isRegionBoundarySimplex();
                if (temp)
                {
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim,1> outNormal(const Simplex<dim,dim-dmo>& simplexChild)
        {
            Eigen::Matrix<double,dim,1> temp(Eigen::Matrix<double,dim,1>::Zero());
            for(const auto& parent : simplexChild.parents())
            {
                if(parent->isBoundarySimplex())
                {
                    temp+=parent->outNormal();
                }
            }
            const double tempNorm(temp.norm());
            return (tempNorm>0.0)? (temp/tempNorm).eval() : temp;
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim,1> outNormal(const Simplex<dim,dim-dmo>& simplexChild,
                                                     const int& rID)
        {
            Eigen::Matrix<double,dim,1> temp(Eigen::Matrix<double,dim,1>::Zero());
            for(const auto& parent : simplexChild.parents())
            {
                    temp+=parent->outNormal(rID);
            }
            const double tempNorm(temp.norm());
            return (tempNorm>0.0)? (temp/tempNorm).eval() : temp;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim>
    struct BoundarySimplex<dim,1>     // this is a Simplex<dim,dim-1>
    {
        /**********************************************************************/
        static bool isBoundarySimplex(const SimplexChild<dim,dim-1>& simplexChild)
        {
            return simplexChild.parents().size()==1;
        }
        
        /**********************************************************************/
        static bool isRegionBoundarySimplex(const SimplexChild<dim,dim-1>& simplexChild)
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
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim,1> outNormal(const Simplex<dim,dim-1>& simplexChild)
        {
            Eigen::Matrix<double,dim,1> temp(Eigen::Matrix<double,dim,1>::Zero());
            if(simplexChild.isBoundarySimplex())
            {
                const size_t faceID=(*simplexChild.parents().begin())->childOrder(simplexChild.xID);
                temp=(*simplexChild.parents().begin())->nda.col(faceID).normalized();
            }
            return temp;
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim,1> outNormal(const Simplex<dim,dim-1>& simplexChild,
                                                     const int& rID)
        {
            Eigen::Matrix<double,dim,1> temp(Eigen::Matrix<double,dim,1>::Zero());
            for(const auto& parent : simplexChild.parents())
            {
                if(  (parent->isBoundarySimplex() || parent->isRegionBoundarySimplex())
                   && parent->region->regionID==rID)
                {
                    const size_t faceID=parent->childOrder(simplexChild.xID);
                    temp+=parent->nda.col(faceID).normalized();
                }
            }
            const double tempNorm(temp.norm());
            return (tempNorm>0.0)? (temp/tempNorm).eval() : temp;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim>
    struct BoundarySimplex<dim,0> // this is a Simplex<dim,dim>
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
