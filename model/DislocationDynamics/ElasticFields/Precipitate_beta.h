/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EshelbyInclusion_H_
#define model_EshelbyInclusion_H_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <Eigen/Dense>
//#include <Material.h>
#include <StaticID.h>

#include <EshelbyInclusion.h>
#include <PeriodicLatticeInterpolant.h>

namespace model
{
    
    
    template <int dim>
    class Precipitate : public EshelbyInclusion<dim>
    {
        typedef EshelbyInclusion<dim> BaseType;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
    
        const std::shared_ptr<PeriodicLatticeInterpolant<2>> gammaSurface;

        
        
    public:
        
        /**********************************************************************/
        Precipitate(const VectorDim& _C,
                         const double& _a,
                         const MatrixDim& _eT,
                         const double& _nu,
                         const double& _mu,
                         const double& _mobilityReduction,
                         const int& _type,
                    const std::shared_ptr<PeriodicLatticeInterpolant<2>>& gammaSurface_in) :
        /* init */ BaseType(_C,_a,_eT,_nu,_mu,_mobilityReduction,_type)
        /* init */,gammaSurface(gammaSurface_in)
        {
            
        }
        
        /**********************************************************************/
        MatrixDim latticeStress(const VectorDim& b) const
        {
            
        }

        
        /**********************************************************************/
        MatrixDim stress(const VectorDim& x,const VectorDim& b, slip system) const
        {
            if(this->contains(x) && gammaSurface)
            {
                return BaseType::stres(x)+latticeStress(b);
            }
            else
            {
                return BaseType::stres(x);
            }
        }
        
        bool isCoherent() const
        {
            return gammaSurface;
        }

    };
}
#endif
