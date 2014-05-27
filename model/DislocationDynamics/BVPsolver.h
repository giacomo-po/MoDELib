/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BVPsolver_H_
#define model_BVPsolver_H_


#include <model/FEM/FiniteElement.h>


namespace model
{
	
	template <int dim, int order>
	struct BVPsolver
    {

        typedef LagrangeElement<dim,order> ElementType;
        typedef FiniteElement<ElementType> FiniteElementType;
        typedef TrialFunction<dim,FiniteElementType> TrialFunctionType;
        typedef TrialGrad<TrialFunctionType> TrialGradType;
        typedef TrialDef<TrialFunctionType> TrialDefType;

//        SimplicialMesh<dim> mesh;
        FiniteElementType fe;
        TrialFunctionType u;      // displacement field u=[u1 u2 u3]'
        TrialGradType b;          // displacement gradient b=[u11 u12 u13 u21 u22 u23 u31 u32 u33]'
        TrialDefType e;           // strain e=[e11 e22 e23 e12 e23 e13]
//        auto s=C*e;               // stress field s=[s11; s12; s21; s22]
        
        /**********************************************************************/
        BVPsolver(const SimplicialMesh<dim>& mesh) :
        /* init  */ fe(mesh),
        /* init  */ u(fe.template trial<dim>()),
        /* init  */ b(grad(u)),
        /* init  */ e(def(u))
        {
        
        }
        
        /**********************************************************************/
        void solve()
        {
         
//            addPenaltyConstraints();
//            ...
//            removePenaltyConstraints();
            
            
        }
        
//        /**********************************************************************/
//        Eigen::Matrix<double,dim,dim> stress(P) const
//        {
//
//        Eigen::Matrix<double,6,1> tempV(s(P));
//        Eigen::Matrix<double,dim,dim> tempM;
//        tempM(0,0)=tempV(0);
//        tempM(1,1)=tempV(1);
//        tempM(2,2)=tempV(2);
//        tempM(1,0)=tempV(3);
//        tempM(2,1)=tempV(4);
//        tempM(2,0)=tempV(5);
//        tempM(0,1)=tempM(1,0);
//        tempM(1,2)=tempM(2,1);
//        tempM(0,2)=tempM(2,0);
//        }
        
        
        
	};
	

} // namespace model
#endif

