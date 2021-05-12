/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BVPsolverBase_H_
#define model_BVPsolverBase_H_

//#define _MODEL_TEST_DD_DISPLACEMENT_
//#define _MODEL_TEST_DD_STRESS_

//#include <Eigen/Dense>
//#include <Eigen/Sparse>
#include <array>
#include <deque>
#include <chrono>
#include <memory> // unique_ptr
#include <FiniteElement.h>
#include <SimplicialMesh.h>
#include <JGNselector.h>
#include <SingleFieldPoint.h>
#include <FEMnodeEvaluation.h>
#include <FEMfaceEvaluation.h>
#include <OnExternalBoundary.h>
#include <LinearWeakList.h>
#include <TextFileParser.h>



namespace model
{
    
    
    
    template <int dim>
    class BVPsolverBase
    {
        
    public:
        typedef LagrangeElement<dim,2> ElementType;
        typedef FiniteElement<ElementType> FiniteElementType;
        typedef TrialFunction<'u',dim,FiniteElementType> TrialFunctionType;
        
        typedef LagrangeElement<dim, 1> cElementType;   //for concentration, YUE
        typedef FiniteElement<cElementType> cFiniteElementType; //for concentration, YUE
        typedef TrialFunction<'c',1,cFiniteElementType> ConcentrationTrialFunctionType;
        typedef TrialGrad<ConcentrationTrialFunctionType> ConcentrationTrialGradType;
        typedef TrialProd<Constant<double,1,1>,ConcentrationTrialGradType> TrialFluxType;
        typedef BilinearForm<ConcentrationTrialGradType,TrialProd<Constant<double,1,1>,TrialFluxType>> ConcentrationBilinearFormType;
        typedef IntegrationDomain<cFiniteElementType,0,4,GaussLegendre> ConcentrationIntegrationDomainType;
        typedef BilinearWeakForm<ConcentrationBilinearFormType,ConcentrationIntegrationDomainType> ConcentrationBilinearWeakFormType;

        constexpr static int dofPerNode=TrialFunctionType::dofPerNode;
        typedef TrialGrad<TrialFunctionType> TrialGradType;
        typedef TrialDef<TrialFunctionType> TrialDefType;
        typedef Eigen::Matrix<double,6,6> CmatrixType;
        typedef Constant<CmatrixType,6,6> CconstantType;
        typedef TrialProd<CconstantType,TrialDefType> TrialStressType;
        //        typedef BilinearWeakForm<TrialDefType,TrialStressType> BilinearWeakFormType;
        typedef BilinearForm<TrialDefType,TrialStressType> BilinearFormType;
        typedef IntegrationDomain<FiniteElementType,0,4,GaussLegendre> IntegrationDomainType;
        typedef BilinearWeakForm<BilinearFormType,IntegrationDomainType> BilinearWeakFormType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::SparseMatrix<double> SparseMatrixType;
        
        /**********************************************************************/
        static Eigen::Matrix<double,6,6> get_C(const double& mu, const double& nu)
        {
            const double lam=2.0*mu*nu/(1.0-2.0*nu);
            const double C11(lam+2.0*mu);
            const double C12(lam);
            const double C44(mu); // C multiplies engineering strain
            
            Eigen::Matrix<double,6,6> temp;
            temp<<C11, C12, C12, 0.0, 0.0, 0.0,
            /***/ C12, C11, C12, 0.0, 0.0, 0.0,
            /***/ C12, C12, C11, 0.0, 0.0, 0.0,
            /***/ 0.0, 0.0, 0.0, C44, 0.0, 0.0,
            /***/ 0.0, 0.0, 0.0, 0.0, C44, 0.0,
            /***/ 0.0, 0.0, 0.0, 0.0, 0.0, C44;
            return temp;
        }
        
        const Polycrystal<dim>& poly;
        Eigen::Matrix<double,6,6> C; // matrix of elastic moduli
        FiniteElementType fe;
        cFiniteElementType cfe;// for concentration, YUE
        TrialFunctionType  u;  // displacement field u=[u1 u2 u3]'
        TrialGradType  b;      // displacement gradient b=[u11 u12 u13 u21 u22 u23 u31 u32 u33]'
        TrialDefType  e;       // strain e=[e11 e22 e33 e12 e23 e13]'
        TrialStressType s;     // stress s=[s11 s22 s33 s12 s23 s13]'
        
        ConcentrationTrialFunctionType c;
        ConcentrationTrialGradType gradc;
        TrialFluxType j;
        
        
    public:
        
        
        /**********************************************************************/
        BVPsolverBase(const Polycrystal<dim>& poly_in) :
        /* init  */ poly(poly_in)
        /* init  */,C(get_C(poly.mu,poly.nu))
        /* init  */,fe(poly.mesh)
        /* init  */,cfe(poly.mesh)
        /* init  */,u(fe.template trial<'u',dim>())
        /* init  */,b(grad(u))
        /* init  */,e(def(u))
        /* init  */,s(C*e)
        /* init  */,c(cfe.template trial<'c',1>())
        /* init  */,gradc(grad(c))
        /* init  */,j(-poly.Dv/poly.Omega*gradc)
        {
            
            
            
            
        }
        
        
        /**********************************************************************/
        const FiniteElementType& finiteElement() const
        {
            return fe;
        }
        
        /**********************************************************************/
        FiniteElementType& finiteElement()
        {
            return fe;
        }
        
        /**********************************************************************/
        const TrialFunctionType& displacement() const
        {
            return u;
        }
        
        /**********************************************************************/
        TrialFunctionType& displacement()
        {
            return u;
        }
        
        /**********************************************************************/
        const TrialStressType& stress() const
        {
            return s;
        }
        
        ConcentrationTrialFunctionType& vacancyConcentration()
        {
            return c;
        }
        
        const ConcentrationTrialFunctionType& vacancyConcentration() const
        {
            return c;
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,1> displacement(const Eigen::Matrix<double,dim,1> P,
                                                 const Simplex<dim,dim>* guess) const
        {
            return eval(u)(P,guess);
        }
        
        /**********************************************************************/
        double vacancyConcentration(const Eigen::Matrix<double,dim,1> P,
                                    const Simplex<dim,dim>* guess) const
        {
            return eval(c)(P,guess)(0,0);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,dim> stress(const Eigen::Matrix<double,dim,1> P,
                                             const Simplex<dim,dim>* guess) const
        {
            Eigen::Matrix<double,6,1> tempV(eval(s)(P,guess));
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
            Eigen::Matrix<double,6,1> tempV(eval(s)(ele,bary));
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

    };
    
}
#endif


//        /**********************************************************************/
//        Eigen::VectorXd solve(const Eigen::VectorXd& b,const TrialFunctionType& guess)
//        {
//            return solve(b,guess.dofVector());
//        }
//
//        /**********************************************************************/
//        Eigen::VectorXd solve(const Eigen::VectorXd& b,const Eigen::VectorXd& y)
//        {/*!@param[in] b the rhs of A*x=b
//          * @param[in] y the guess for x
//          */
//
//            model::cout<<"Setting up null-space solver..."<<std::flush;
//            const auto t0= std::chrono::system_clock::now();
//
//            const size_t cSize(displacement().dirichletConditions().size());
//            const size_t tSize = gSize-cSize;
//
//            std::vector<Eigen::Triplet<double> > tTriplets;
//            tTriplets.reserve(tSize);
//
//            Eigen::VectorXd g(Eigen::VectorXd::Zero(gSize));
//            Eigen::VectorXd guess(Eigen::VectorXd::Zero(tSize));
//
//            size_t startRow=0;
//            size_t col=0;
//
//            for (const auto& cIter : displacement().dirichletConditions())
//            {
//                const size_t& endRow = cIter.first;
//                g(endRow)= cIter.second;
//                for (size_t row=startRow;row!=endRow;++row)
//                {
//                    tTriplets.emplace_back(row,col,1.0);
//                    guess(col)=y(row);
//                    col++;
//                }
//                startRow=endRow+1;
//            }
//            for (size_t row=startRow;row!=gSize;++row)
//            {
//                tTriplets.emplace_back(row,col,1.0);
//                guess(col)=y(row);
//                col++;
//            }
//
//            SparseMatrixType T1(gSize,tSize);
//            T1.setFromTriplets(tTriplets.begin(),tTriplets.end());
//            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
//
//            // Check if null-space has changed
//            const auto t4= std::chrono::system_clock::now();
//            model::cout<<"Checking if null-space has changed... "<<std::flush;
//            bool sameT=false;
//            if(T.cols()==T1.cols() && T.rows()==T1.rows())
//            {
//                const double normT2=(T1-T).squaredNorm();
//                if(normT2==0.0)
//                {
//                    sameT=true;
//                }
//
//            }
//            model::cout<<!sameT<<std::flush;
//            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t4)).count()<<" sec]"<<std::endl;
//
//            // Update A1 if null-space has changed
//            if(!sameT)
//            { // need to re-factorized the LLT decomposition
//                const auto t1= std::chrono::system_clock::now();
//                if(use_directSolver)
//                {
//#ifdef _MODEL_PARDISO_SOLVER_
//                    model::cout<<"PardisoLLT: factorizing..."<<std::flush;
//#else
//                    model::cout<<"SimplicialLLT: factorizing..."<<std::flush;
//#endif
//                }
//                else
//                {
//                    model::cout<<"ConjugateGradient: factorizing..."<<std::flush;
//                }
//
//                T=T1; // store new T
//                A1=T.transpose()*A*T; // store new A1
//                if(use_directSolver)
//                {
//                    directSolver.compute(A1);
//                    assert(directSolver.info()==Eigen::Success && "SOLVER  FAILED");
//                }
//                else
//                {
//                    iterativeSolver.compute(A1);
//                    assert(iterativeSolver.info()==Eigen::Success && "SOLVER  FAILED");
//                }
//
//                model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<std::endl;
//            }
//
//            // Solve
//            const auto t2= std::chrono::system_clock::now();
//            Eigen::VectorXd b1(T.transpose()*(b-A*g));
//            //            Eigen::VectorXd x(Eigen::VectorXd::Zero(gSize));
//            Eigen::VectorXd x(Eigen::VectorXd::Zero(b1.rows()));
//
//
//            if(use_directSolver)
//            {
//#ifdef _MODEL_PARDISO_SOLVER_
//                model::cout<<"PardisoLLT: solving..."<<std::flush;
//#else
//                model::cout<<"SimplicialLLT: solving..."<<std::flush;
//#endif
//                //                x=T*directSolver.solve(b1)+g;
//                x=directSolver.solve(b1);
//                assert(directSolver.info()==Eigen::Success && "SOLVER  FAILED");
//                const double b1Norm(b1.norm());
//                if(b1Norm>0.0)
//                {
//                    const double axbNorm((A1*x-b1).norm());
//                    if(axbNorm/b1Norm>tolerance)
//                    {
//                        model::cout<<"norm(A*x-b)/norm(b)="<<axbNorm/b1Norm<<std::endl;
//                        model::cout<<"tolerance="<<tolerance<<std::endl;
//                        assert(0 && "SOLVER FAILED");
//                    }
//                }
//            }
//            else
//            {
//                model::cout<<"ConjugateGradient: solving..."<<std::flush;
//                iterativeSolver.setTolerance(tolerance);
//                //                x=T*iterativeSolver.solveWithGuess(b1,guess)+g;
//                x=iterativeSolver.solveWithGuess(b1,guess);
//                model::cout<<" (relative error ="<<iterativeSolver.error()<<", tolerance="<<iterativeSolver.tolerance()<<")";
//                assert(iterativeSolver.info()==Eigen::Success && "SOLVER  FAILED");
//            }
//
//
//
//
//            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]"<<std::endl;
//            return T*x+g;
//        }

//        /**********************************************************************/
//        Eigen::VectorXd solveRemovingRigidBodyMotion(const Eigen::VectorXd& b)
//        {
//
//            SparseMatrixType C;
//            // fill C
//
//            NullSpaceSolver<SparseMatrixType> nss(A,C);
//            return nss.solve(b);
//
//        }
