/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_FEMdomainEvaluation_H_
#define model_FEMdomainEvaluation_H_

#include <Eigen/Dense>
#include <FEMbaseEvaluation.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename ElementType,int rows,int cols>
    struct FEMdomainEvaluation : public FEMbaseEvaluation<ElementType,rows,cols>
    {
        constexpr static int dim=ElementType::dim;
        const ElementType& ele;
        const Eigen::Matrix<double,dim+1,1> domainBary;
//        const int boundarydomain;
        const double& weight;
//        const Eigen::Matrix<double,dim,1> P;
        
        /**********************************************************************/
        FEMdomainEvaluation(const ElementType& _ele,
                          const Eigen::Matrix<double,dim+1,1>& _domainBary,
                          const int& _boundarydomain,
                          const double& _weight) :
        /* init */ FEMbaseEvaluation(_ele.position(_domainBary))
        //        /* init */ Eigen::Matrix<double,rows,cols>(Eigen::Matrix<double,rows,cols>::Zero()),
        /* init */,ele(_ele)
        /* init */,domainBary(_domainBary)
//        /* init */,boundarydomain(_boundarydomain)
        /* init */,weight(_weight)
        //        /* init */ P(ele.position(domainBary))
        {
            
        }
        
        /**********************************************************************/
        const FEMdomainEvaluation& operator=(const Eigen::Matrix<double,rows,cols>& val)
        {
            static_cast<Eigen::Matrix<double,rows,cols>*>(this)->operator=(val);
            return *this;
        }
        
    };
    
}    // close namespace
#endif
