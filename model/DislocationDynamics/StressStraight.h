/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_StressStraight_H_
#define model_StressStraight_H_

#include <Eigen/Core>
#include <model/DislocationDynamics/Materials/Material.h>
//#include <model/DislocationDynamics/NearestNeighbor/DislocationStress.h>

namespace model
{
    
    template<short unsigned int _dim>
	struct DislocationStress; // class predeclaration
    
    template <int dim>
    class StressStraight
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
        
        const VectorDim P0;
        const VectorDim P1;
        const VectorDim b;
        const VectorDim t;

        /**********************************************************************/
        MatrixDim stress_kernel(const VectorDim& r) const
        {
#if _MODEL_NON_SINGULAR_DD_ == 0
            /* Devincre, B., & Condat, M. (1992). Model validation of a 3D
             * simulation of dislocation dynamics: Discretization and line tension
             * effects. Acta Metallurgica Et Materialia, 40(10), 2629–2637.
             */
            const double L(r.dot(t));
            const double R(r.norm()+DislocationStress<dim>::a);
            const VectorDim rho(r-L*t);
            const VectorDim Y((L+R)*t+rho);
            const double Y2(Y.squaredNorm()+DislocationStress<dim>::a2);
            return (Material<Isotropic>::C1* b.cross(Y)*t.transpose()
            /*                 */ -b.cross(t)*Y.transpose()
            /*                 */ -b.dot(Y.cross(t))*(2.0/Y2*rho*Y.transpose()+0.5*(MatrixDim::Identity()+t*t.transpose()+2.0/Y2*L/R*Y*Y.transpose()))
                    )*2.0/Y2;
            
#elif _MODEL_NON_SINGULAR_DD_ == 1 /* Cai's non-singular theory */
            assert(0 && "NOT IMPLEMENTED YET");
#elif _MODEL_NON_SINGULAR_DD_ == 2 /* Lazar's non-singular theory */
            assert(0 && "NOT IMPLEMENTED YET");
#else
#error Unsupported choice of field regularization
#endif
        }
        
    public:
        /**********************************************************************/
        StressStraight(const VectorDim& _P0,const VectorDim& _P1, const VectorDim& _b) :
        /* init list */ P0(_P0),
        /* init list */ P1(_P1),
        /* init list */ b(_b),
        /* init list */ t((P1-P0).normalized())
        {
        
        }

        /**********************************************************************/
        MatrixDim nonSymmStress(const VectorDim& x) const
        {
            return stress_kernel(x-P1)-stress_kernel(x-P0);
            
        }
        
        /**********************************************************************/
        MatrixDim stress(const VectorDim& x) const
        {
            
            const MatrixDim temp = nonSymmStress(x);

            return Material<Isotropic>::C2*(temp+temp.transpose());
            
        }
        
	};	
	
	/*********************************************************************/
	/*********************************************************************/
} // end namespace
#endif

