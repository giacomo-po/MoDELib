/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BVPsolver_H_
#define model_BVPsolver_H_

#include <memory> // unique_ptr
#include <model/FEM/FiniteElement.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/DislocationDynamics/BVP/DislocationNegativeStress.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/FEM/WeakForms/JGNselector.h>

namespace model
{
	
	template <int dim, int sfOrder>
	class BVPsolver
    {
        
        typedef LagrangeElement<dim,sfOrder> ElementType;
        typedef FiniteElement<ElementType> FiniteElementType;
        typedef TrialFunction<dim,FiniteElementType> TrialFunctionType;
        typedef TrialGrad<TrialFunctionType> TrialGradType;
        typedef TrialDef<TrialFunctionType> TrialDefType;
        typedef Eigen::Matrix<double,6,6> CmatrixType;
        typedef Constant<CmatrixType,6,6> CconstantType;
        typedef TrialProd<CconstantType,TrialDefType> TrialStressType;
//        typedef BilinearWeakForm<TrialDefType,TrialStressType> BilinearWeakFormType;
        typedef IntegrationDomain<FiniteElementType,0,4,GaussLegendre> IntegrationDomainType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        //        SimplicialMesh<dim> mesh;
        
        const SimplicialMesh<dim>& mesh;
        Eigen::Matrix<double,6,6> C;
        
        
        
        std::unique_ptr<FiniteElementType> fe;
        std::unique_ptr<TrialFunctionType>  u;  // displacement field *u=[u1 u2 u3]'
        std::unique_ptr<TrialGradType>  b;      // displacement gradient *b=[u11 u12 u13 u21 u22 u23 u31 u32 u33]'
        std::unique_ptr<TrialDefType>  e;       // strain *e=[e11 e22 e33 e12 e23 e13]'
        std::unique_ptr<TrialStressType> s;     // stress *s=[s11 s22 s33 s12 s23 s13]'
        IntegrationDomainType dV;
        
        
        /**********************************************************************/
        Eigen::Matrix<double,dim+1,1> face2domainBary(const Eigen::Matrix<double,dim,1>& b1,
                                                      /*                                         */ const int& boundaryFace) const
        {
            // Transform to barycentric coordinate on the volume, adding a zero on the boundaryFace-face
            Eigen::Matrix<double,dim+1,1> bary;
            for (int k=0;k<dim;++k)
            {
                bary((k<boundaryFace)? k : k+1)=b1(k);
            }
            bary(boundaryFace)=0.0;
            return bary;
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,6,6> get_C() const
        {
            const double  mu=1.0;  // dimensionless
            const double  nu=Material<Isotropic>::nu;
            const double lam=2.0*mu*nu/(1.0-2.0*nu);
            const double C11(lam+2.0*mu);
            const double C12(lam);
            const double C44(2.0*mu); // C multiplies true strain (not engineering), so 2 is necessary
            
            Eigen::Matrix<double,6,6> temp;
            temp<<C11, C12, C12, 0.0, 0.0, 0.0,
            /***/ C12, C11, C12, 0.0, 0.0, 0.0,
            /***/ C12, C12, C11, 0.0, 0.0, 0.0,
            /***/ 0.0, 0.0, 0.0, C44, 0.0, 0.0,
            /***/ 0.0, 0.0, 0.0, 0.0, C44, 0.0,
            /***/ 0.0, 0.0, 0.0, 0.0, 0.0, C44;
            return temp;
        }
        
    public:
        /**********************************************************************/
        BVPsolver(const SimplicialMesh<dim>& mesh_in) :
        /* init  */ mesh(mesh_in),
        /* init  */ C(get_C())
        {
            
        }
        
        
        /**********************************************************************/
        const FiniteElementType& finiteElement() const
        {
            return *fe;
        }
        
        /**********************************************************************/
        FiniteElementType& finiteElement()
        {
            return *fe;
        }
        
        /**********************************************************************/
        const TrialFunctionType& displacement() const
        {
            return *u;
        }
        
        /**********************************************************************/
        TrialFunctionType& displacement()
        {
            return *u;
        }
        
//        /**********************************************************************/
//        const TrialStressType& stress() const
//        {
//            return *s;
//        }
        
        /**********************************************************************/
        void init()
        {
            fe = std::move(std::unique_ptr<FiniteElementType>(new FiniteElementType(mesh)));
            u  = std::move(std::unique_ptr<TrialFunctionType>(new TrialFunctionType(fe->template trial<dim>())));
            b  = std::move(std::unique_ptr<TrialGradType>(new TrialGradType(grad(*u))));
            e  = std::move(std::unique_ptr<TrialDefType>(new TrialDefType(def(*u))));
            C=get_C(); // Material<Isotropic>  may have changed
            s  = std::move(std::unique_ptr<TrialStressType>(new TrialStressType(C**e)));
//            _bWF = std::move(std::unique_ptr<BilinearWeakFormType>(new BilinearWeakFormType(e->test(),*s)));
            
            dV = fe->template domain<EntireDomain,4,GaussLegendre>();
            
//            auto dV=fe->domain<EntireDomain,4,GaussLegendre>();
            
//            (*_bWF)*dV;
            
        }
        
        
        /**********************************************************************/
        template <typename DislocationNetworkType,int qOrder>
        void assembleAndSolve(const DislocationNetworkType& DN)
        {
            typedef DislocationNegativeStress<DislocationNetworkType> DislocationNegativeStressType;

            const DislocationNegativeStressType ds(DN);
            
            //typedef LinearWeakForm<TrialFunctionType,DislocationNegativeStressType> LinearWeakFormType;

            auto bWF=(e->test(),*s)*dV;

            
            //LinearWeakFormType lwf(u->test(),ds);
            auto ndA=fe->template boundary<ExternalBoundary,qOrder,GaussLegendre>();
            
            // Create the LinearWeakForm lWF_u=int(test(u)^T*f)ndS
            auto lWF=(u->test(),ds)*ndA;
            
            // Create the WeakProblem
            auto wp(bWF=lWF); //  weak problem
            
            //wp.assembleWithLagrangeConstraints();
            wp.assembleWithPenaltyConstraints(1000.0);
            
            wp.solve(0.0001);
            //wp_u.output();
            u->dofContainer=wp.x.segment(0,u->dofContainer.size());
            
            
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,dim> stress(const Eigen::Matrix<double,dim,1> P,
                                             const Simplex<dim,dim>* guess) const
        {
            
            Eigen::Matrix<double,6,1> tempV((*s)(P,guess));
            Eigen::Matrix<double,dim,dim> tempM;
            tempM(0,0)=tempV(0); // s11
            tempM(1,1)=tempV(1); // s22
            tempM(2,2)=tempV(2); // s33
            tempM(1,0)=tempV(3); // s21
            tempM(2,1)=tempV(4); // s32
            tempM(2,0)=tempV(5); // s31
            tempM(0,1)=tempM(1,0); //symm
            tempM(1,2)=tempM(2,1); //symm
            tempM(0,2)=tempM(2,0); //symm
            
            return tempM;
        }

        /**********************************************************************/
        Eigen::Matrix<double,dim,dim> stress(const ElementType& ele,
                                             const Eigen::Matrix<double,dim+1,1>& bary) const
        {
            
            Eigen::Matrix<double,6,1> tempV((*s)(ele,bary));
            Eigen::Matrix<double,dim,dim> tempM;
            tempM(0,0)=tempV(0); // s11
            tempM(1,1)=tempV(1); // s22
            tempM(2,2)=tempV(2); // s33
            tempM(1,0)=tempV(3); // s21
            tempM(2,1)=tempV(4); // s32
            tempM(2,0)=tempV(5); // s31
            tempM(0,1)=tempM(1,0); //symm
            tempM(1,2)=tempM(2,1); //symm
            tempM(0,2)=tempM(2,0); //symm
            
            return tempM;
        }
        
        /**********************************************************************/
		template <typename DislocationNetworkType>
        VectorDim traction(const Eigen::Matrix<double,dim-1,1>& a1, const ElementType& ele, const int& boundaryFace, const DislocationNetworkType& DN) const
        {
            const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(a1));
            const Eigen::Matrix<double,dim+1,1> bary(face2domainBary(b1,boundaryFace));
            return (stress(ele,bary)+DN.stress(ele.position(bary)))*JGNselector<dim>::jGN(ele.jGN(bary,boundaryFace));
		}
        

        
	};
	
    
} // namespace model
#endif


