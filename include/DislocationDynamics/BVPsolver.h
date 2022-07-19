/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BVPsolver_H_
#define model_BVPsolver_H_

//#define _MODEL_TEST_DD_DISPLACEMENT_
//#define _MODEL_TEST_DD_STRESS_

#include <array>
#include <deque>
#include <chrono>
#include <memory> // unique_ptr
#ifdef _MODEL_PARDISO_SOLVER_
#include <Eigen/PardisoSupport>
#endif
#include<Eigen/SparseCholesky>
#include <FiniteElement.h>
//#include <Material.h>
//#include <DislocationNegativeFields.h>
#include <SimplicialMesh.h>
#include <JGNselector.h>
//#include <SingleFieldPoint.h>
#include <FEMnodeEvaluation.h>
#include <FEMfaceEvaluation.h>

//#include <BoundaryDisplacementPoint.h>
//#include <DisplacementPoint.h>
//#include <BoundaryStressPoint.h>
#include <LinearWeakList.h>
//#include <SparseNullSpace.h>
#include <NullSpaceSolver.h>

//#include <BoundaryQuadraturePoint.h>
//
//#include <RuntimeError.h>
#include <TextFileParser.h>

#ifndef userLoadController
#define userLoadController "DummyLoadController.h"
#endif
#include userLoadController


namespace model
{
    
    template <int dim, int sfOrder>
    class BVPsolver
    {
        
        
        
    public:
        
        typedef LagrangeElement<dim,sfOrder> ElementType;
        typedef FiniteElement<ElementType> FiniteElementType;
        typedef TrialFunction<'u',dim,FiniteElementType> TrialFunctionType;
        
        typedef LagrangeElement<dim, 1> cElementType;   //for concentration, YUE
        typedef FiniteElement<cElementType> cFiniteElementType; //for concentration, YUE
        typedef TrialFunction<'c',1,cFiniteElementType> ConcentrationTrialFunctionType;

        //        typedef typename TrialFunctionType::DirichletConditionContainerType DirichletConditionContainerType;
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
        
        typedef LoadController<TrialFunctionType> LoadControllerType;
        
        const SimplicialMesh<dim>& mesh;
        size_t gSize;
        
        //        static bool apply_DD_displacement;
        
    private:
        Eigen::Matrix<double,6,6> C; // matrix of elastic moduli
        FiniteElementType fe;
        cFiniteElementType cfe;// for concentration, YUE
        TrialFunctionType  u;  // displacement field u=[u1 u2 u3]'
        TrialGradType  b;      // displacement gradient b=[u11 u12 u13 u21 u22 u23 u31 u32 u33]'
        TrialDefType  e;       // strain e=[e11 e22 e33 e12 e23 e13]'
        TrialStressType s;     // stress s=[s11 s22 s33 s12 s23 s13]'
        
        ConcentrationTrialFunctionType c;
        
        IntegrationDomainType dV;
        BilinearWeakFormType* bWF;
        
        LoadControllerType* lc;
        
        SparseMatrixType A;
        SparseMatrixType T;
        SparseMatrixType A1;
        
#ifdef _MODEL_PARDISO_SOLVER_
        Eigen::PardisoLLT<SparseMatrixType> directSolver;
#else
        Eigen::SimplicialLLT<SparseMatrixType> directSolver;
#endif
        
        Eigen::ConjugateGradient<SparseMatrixType> iterativeSolver;
        
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
        Eigen::Matrix<double,6,6> get_C(const double& mu, const double& nu) const
        {
//            const double  mu=1.0;  // dimensionless
//            const double  nu=Material<dim,Isotropic>::nu;
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

        /**********************************************************************/
        Eigen::VectorXd solve(const Eigen::VectorXd& b,const TrialFunctionType& guess)
        {
            return solve(b,guess.dofVector());
        }
        
        /**********************************************************************/
        Eigen::VectorXd solve(const Eigen::VectorXd& b,const Eigen::VectorXd& y)
        {/*!@param[in] b the rhs of A*x=b
          * @param[in] y the guess for x
          */
            
            std::cout<<"Setting up null-space solver..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            
            const size_t cSize(displacement().dirichletConditions().size());
            const size_t tSize = gSize-cSize;
            
            std::vector<Eigen::Triplet<double> > tTriplets;
            tTriplets.reserve(tSize);
            
            Eigen::VectorXd g(Eigen::VectorXd::Zero(gSize));
            Eigen::VectorXd guess(Eigen::VectorXd::Zero(tSize));
            
            size_t startRow=0;
            size_t col=0;
            
            for (const auto& cIter : displacement().dirichletConditions())
            {
                const size_t& endRow = cIter.first;
                g(endRow)= cIter.second;
                for (size_t row=startRow;row!=endRow;++row)
                {
                    tTriplets.emplace_back(row,col,1.0);
                    guess(col)=y(row);
                    col++;
                }
                startRow=endRow+1;
            }
            for (size_t row=startRow;row!=gSize;++row)
            {
                tTriplets.emplace_back(row,col,1.0);
                guess(col)=y(row);
                col++;
            }
            
            SparseMatrixType T1(gSize,tSize);
            T1.setFromTriplets(tTriplets.begin(),tTriplets.end());
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            
            // Check if null-space has changed
            const auto t4= std::chrono::system_clock::now();
            std::cout<<"Checking if null-space has changed... "<<std::flush;
            bool sameT=false;
            if(T.cols()==T1.cols() && T.rows()==T1.rows())
            {
                const double normT2=(T1-T).squaredNorm();
                if(normT2==0.0)
                {
                    sameT=true;
                }
                
            }
            std::cout<<!sameT<<std::flush;
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t4)).count()<<" sec]"<<std::endl;
            
            // Update A1 if null-space has changed
            if(!sameT)
            { // need to re-factorized the LLT decomposition
                const auto t1= std::chrono::system_clock::now();
                if(use_directSolver)
                {
#ifdef _MODEL_PARDISO_SOLVER_
                    std::cout<<"PardisoLLT: factorizing..."<<std::flush;
#else
                    std::cout<<"SimplicialLLT: factorizing..."<<std::flush;
#endif
                }
                else
                {
                    std::cout<<"ConjugateGradient: factorizing..."<<std::flush;
                }
                
                T=T1; // store new T
                A1=T.transpose()*A*T; // store new A1
                if(use_directSolver)
                {
                    directSolver.compute(A1);
                    assert(directSolver.info()==Eigen::Success && "SOLVER  FAILED");
                }
                else
                {
                    iterativeSolver.compute(A1);
                    assert(iterativeSolver.info()==Eigen::Success && "SOLVER  FAILED");
                }

                std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<std::endl;
            }
            
            // Solve
            const auto t2= std::chrono::system_clock::now();
            Eigen::VectorXd b1(T.transpose()*(b-A*g));
            //            Eigen::VectorXd x(Eigen::VectorXd::Zero(gSize));
            Eigen::VectorXd x(Eigen::VectorXd::Zero(b1.rows()));
            
            
            if(use_directSolver)
            {
#ifdef _MODEL_PARDISO_SOLVER_
                std::cout<<"PardisoLLT: solving..."<<std::flush;
#else
                std::cout<<"SimplicialLLT: solving..."<<std::flush;
#endif
                //                x=T*directSolver.solve(b1)+g;
                x=directSolver.solve(b1);
                assert(directSolver.info()==Eigen::Success && "SOLVER  FAILED");
                const double b1Norm(b1.norm());
                if(b1Norm>0.0)
                {
                    const double axbNorm((A1*x-b1).norm());
                    if(axbNorm/b1Norm>tolerance)
                    {
                        std::cout<<"norm(A*x-b)/norm(b)="<<axbNorm/b1Norm<<std::endl;
                        std::cout<<"tolerance="<<tolerance<<std::endl;
                        assert(0 && "SOLVER FAILED");
                    }
                }
            }
            else
            {
                std::cout<<"ConjugateGradient: solving..."<<std::flush;
                iterativeSolver.setTolerance(tolerance);
                //                x=T*iterativeSolver.solveWithGuess(b1,guess)+g;
                x=iterativeSolver.solveWithGuess(b1,guess);
                std::cout<<" (relative error ="<<iterativeSolver.error()<<", tolerance="<<iterativeSolver.tolerance()<<")";
                assert(iterativeSolver.info()==Eigen::Success && "SOLVER  FAILED");
            }
            
            
            
            
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]"<<std::endl;
            return T*x+g;
        }
        
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
        
        
    public:
        
        double tolerance;
        const bool use_directSolver;
        const int stepsBetweenBVPupdates;
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        BVPsolver(const SimplicialMesh<dim>& mesh_in,const DislocationNetworkType& DN) :
        /* init  */ mesh(mesh_in)
        /* init  */,gSize(0)
        /* init  */,C(get_C(DN.poly.mu,DN.poly.nu))
        /* init  */,fe(mesh)
        /* init  */,cfe(mesh)
        /* init  */,u(fe.template trial<'u',dim>())
        /* init  */,b(grad(u))
        /* init  */,e(def(u))
        /* init  */,s(C*e)
        /* init  */,c(cfe.template trial<'c',1>())
        /* init  */,tolerance(TextFileParser(DN.simulationParameters.traitsIO.ddFile).readScalar<double>("solverTolerance",true))
        /* init  */,use_directSolver(TextFileParser(DN.simulationParameters.traitsIO.ddFile).readScalar<int>("use_directSolver_FEM",true))
        /* init  */,stepsBetweenBVPupdates(TextFileParser(DN.simulationParameters.traitsIO.ddFile).readScalar<int>("stepsBetweenBVPupdates",true))
        {
            
            
            std::cout<<"Initializing BVPsolver"<<std::endl;
            
            //            fe = new FiniteElementType(mesh);
            //            u  = new TrialFunctionType(fe.template trial<'u',dim>());
            //            b  = new TrialGradType(grad(u));
            //            e  = new TrialDefType(def(u));
            C=get_C(DN.poly.mu,DN.poly.nu); // Material<Isotropic>  may have changed since construction
//            s  = new TrialStressType(C*e);
            dV = fe.template domain<EntireDomain,4,GaussLegendre>();
            bWF = new BilinearWeakFormType((test(e),s),dV);
            gSize=TrialBase<TrialFunctionType>::gSize();
            
            // Compute and store stiffness matrix
            const auto t0= std::chrono::system_clock::now();
            std::cout<<"Computing stiffness matrix: gSize="<<gSize<<std::endl;
            std::vector<Eigen::Triplet<double> > globalTriplets(bWF->globalTriplets());
            A.resize(gSize,gSize);
            A.setFromTriplets(globalTriplets.begin(),globalTriplets.end());
            A.prune(A.norm()/A.nonZeros(),FLT_EPSILON);
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            
            // Initialize LoadController
            lc = new LoadControllerType(displacement());
            lc ->init(DN);
            
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
        const LoadControllerType& loadController() const
        {
            return *lc;
        }
        
        /**********************************************************************/
        template<typename DislocationNetworkType>
        void modifyDirichletConditions(const DislocationNetworkType& DN)
        {/*!@param[in] dc the Dirichlet condition
          * @param[in] the nodal dof to be constrained
          */
            
#ifdef _MODEL_TEST_DD_DISPLACEMENT_
            SequentialOutputFile<'Y',1>::set_count(DN.runningID());
            SequentialOutputFile<'Y',1> dd_disp_file;
#endif
            
            std::cout<<"subtracting DislocationDisplacement... "<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            
            // Add the (negative) dislocation displacement
            std::vector<FEMnodeEvaluation<ElementType,dim,1>> fieldPoints;
            fieldPoints.reserve(displacement().dirichletNodeMap().size());
            
            for (const auto& pair : displacement().dirichletNodeMap()) // range-based for loop (C++11)
            {
                const auto& gID(pair.first);
                                fieldPoints.emplace_back(gID,fe.node(gID).P0);
            }
  
            DN.displacement(fieldPoints);
            
            // Subtract the DislocationNetwork displacement from the Dirichlet conditions
            for(const auto& fieldPoint : fieldPoints)
            {
                const size_t& gID=fieldPoint.pointID;
                for(int dof=0;dof<dofPerNode;++dof)
                {
                    if(displacement().dirichletNodeMap()[gID][dof])
                    {
#ifdef _MODEL_TEST_DD_DISPLACEMENT_
                        dd_disp_file<<fieldPoint.P.transpose()<<" "<<fieldPoint.template field<DisplacementField>().transpose()<<"\t"<<fieldPoint.S.transpose()<<"\n";
#endif
                        u.dirichletConditions().at(gID*dofPerNode+dof) -= fieldPoint(dof);

                    }
                }
            }
            
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        
        
//#ifdef userBVPfile
#ifdef _MODEL_GREATWHITE_
#include <BVPsolverGreatWhite.h>
#else
        /**********************************************************************/
        template <typename DislocationNetworkType,int qOrder>
        void assembleAndSolve(const DislocationNetworkType& DN, const bool& )
        {
            // Clear exisintg Dirichlet condition
            u.clearDirichletConditions();
            
            // Update loadController
            lc->update(DN);
            
            // Add Dirichlet conditions by loadController
            lc->addDirichletConditions(DN);
            
            // Modify Dirichlet conditions by subtracting dislocation displacement
            modifyDirichletConditions(DN);
            
            // Compute dislocation traction
            std::cout<<"Computing DD boundary traction..."<<std::flush;
            auto ndA=fe.template boundary<ExternalBoundary,qOrder,GaussLegendre>();
            auto eb_list = ndA.template integrationList<FEMfaceEvaluation<ElementType,dim,dim>>(); // TO DO: make this a member data to be able to output
            
            const auto t0= std::chrono::system_clock::now();
            if(TypeTraits<DislocationNetworkType>::corder==0)
            {
                //            REIMPLEMENT DD TRACTION
                //                assert(0 && "REIMPLEMENT DD TRACTION");
                //
                DN.stress(eb_list);
            }
            else
            {
                assert(0 && "REIMPLEMENT THIS");
            }
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
            
            auto dislocationTraction=(test(u),eb_list);
            
            // Assemble loadController and dislocaiton tractions and solve
            displacement()=solve(lc->globalVector(DN)-dislocationTraction.globalVector(),displacement());
//                        displacement()=solveRemovingRigidBodyMotion(lc->globalVector(DN)-dislocationTraction.globalVector());
            
        }
#endif

//#else
//        /**********************************************************************/
//        template <typename DislocationNetworkType,int qOrder>
//        void assembleAndSolve(const DislocationNetworkType& )
//        {
//            assert(0 && "YOU MUST #define THE userBVPfile to use BVPsolver.");
//        }
//#endif
        
        

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
        
        /**********************************************************************/
        VectorDim bvpTraction(const Eigen::Matrix<double,dim-1,1>& a1, const ElementType& ele, const int& boundaryFace) const
        {
            const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(a1));
            const Eigen::Matrix<double,dim+1,1> bary(face2domainBary(b1,boundaryFace));
            return stress(ele,bary)*JGNselector<dim>::jGN(ele.jGN(bary,boundaryFace));
        }
        
        /**********************************************************************/
        VectorDim bvpMoment(const Eigen::Matrix<double,dim-1,1>& a1, const ElementType& ele, const int& boundaryFace,const VectorDim& x0) const
        {
            const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(a1));
            const Eigen::Matrix<double,dim+1,1> bary(face2domainBary(b1,boundaryFace));
            return (ele.position(bary)-x0).cross(stress(ele,bary)*JGNselector<dim>::jGN(ele.jGN(bary,boundaryFace)));
        }
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        VectorDim ddTraction(const Eigen::Matrix<double,dim-1,1>& a1, const ElementType& ele, const int& boundaryFace, const DislocationNetworkType& DN) const
        {
            const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(a1));
            const Eigen::Matrix<double,dim+1,1> bary(face2domainBary(b1,boundaryFace));
            return DN.stress(ele.position(bary))*JGNselector<dim>::jGN(ele.jGN(bary,boundaryFace));
        }
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        VectorDim ddMoment(const Eigen::Matrix<double,dim-1,1>& a1, const ElementType& ele, const int& boundaryFace, const VectorDim& x0,const DislocationNetworkType& DN) const
        {
            const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(a1));
            const Eigen::Matrix<double,dim+1,1> bary(face2domainBary(b1,boundaryFace));
            const VectorDim pos(ele.position(bary));
            return (pos-x0).cross(DN.stress(pos)*JGNselector<dim>::jGN(ele.jGN(bary,boundaryFace)));
        }
        
        /**********************************************************************/
        void output() const
        {
            //#ifdef _MODEL_TEST_DD_STRESS_
            //            SequentialOutputFile<'Z',1>::set_count(DN.runningID());
            //            SequentialOutputFile<'Z',1> dd_stress_file;
            //            for (const auto& point : eb_list)
            //            {
            //            dd_stress_file<<point.P.transpose()<<"\t"<< point.template field<StressField>().row(0)<<"\t"
            //                                                     << point.template field<StressField>().row(1)<<"\t"
            //                                                     << point.template field<StressField>().row(2)<<"\n";
            //            }
            //            
            //#endif
        }
        
    };
}
#endif

