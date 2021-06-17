/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_WeakProblem_H_
#define model_WeakProblem_H_

//#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <chrono>

//#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>


#include <unsupported/Eigen/IterativeSolvers>

#include <ExpressionRef.h>


#ifdef _MODEL_PARDISO_SOLVER_
#include <Eigen/PardisoSupport>
#endif

#ifdef _OPENMP
#include <omp.h> // for PARDISO
#endif



namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename LHS, typename RHS>
	class WeakProblem
    {
        
        //typedef WeakForm<T1,T2,TF> BilinearWeakFormType;
        
        typedef typename LHS::TrialFunctionType TrialFunctionType;
//        typedef typename TrialFunctionType::DirichletConditionContainerType DirichletConditionContainerType;
        typedef typename TrialBase<TrialFunctionType>::DirichletConditionContainerType DirichletConditionContainerType;
        
        
//        const LHS& lhs;
//        const RHS&   rhs;
        
        ExpressionRef<LHS> lhs;
        ExpressionRef<RHS> rhs;
        
        enum {NO_CONSTRAINTS,LAGRANGE,PENALTY,MASTERSLAVE};
        int assembleCase;
        
//        double maxAbsValue;
        
        
    public:
        const size_t gSize;
        
        
        typedef Eigen::SparseMatrix<double> SparseMatrixType;
        SparseMatrixType A;
        Eigen::VectorXd  b;
        Eigen::VectorXd  x;
        
        /**********************************************************************/
        WeakProblem(const LHS& lhs_in, const RHS& rhs_in) :
        /* init list */ lhs(lhs_in),
        /* init list */ rhs(rhs_in),
        /* init list */ assembleCase(NO_CONSTRAINTS),
//        /* init list */ maxAbsValue(0.0),
//        /* init list */ gSize(lhs.gSize)
        /* init list */ gSize(TrialBase<TrialFunctionType>::gSize())

        {
            std::cout<<greenColor<<"Creating WeakProblem "<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        WeakProblem(LHS&& lhs_in, const RHS& rhs_in) :
        /* init list */ lhs(std::move(lhs_in)),
        /* init list */ rhs(rhs_in),
        /* init list */ assembleCase(NO_CONSTRAINTS),
        //        /* init list */ maxAbsValue(0.0),
        //        /* init list */ gSize(lhs.gSize)
        /* init list */ gSize(TrialBase<TrialFunctionType>::gSize())
        
        {
            std::cout<<greenColor<<"Creating WeakProblem "<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        WeakProblem(const LHS& lhs_in, RHS&& rhs_in) :
        /* init list */ lhs(lhs_in),
        /* init list */ rhs(std::move(rhs_in)),
        /* init list */ assembleCase(NO_CONSTRAINTS),
        //        /* init list */ maxAbsValue(0.0),
        //        /* init list */ gSize(lhs.gSize)
        /* init list */ gSize(TrialBase<TrialFunctionType>::gSize())
        
        {
            std::cout<<greenColor<<"Creating WeakProblem "<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        WeakProblem(LHS&& lhs_in, RHS&& rhs_in) :
        /* init list */ lhs(std::move(lhs_in)),
        /* init list */ rhs(std::move(rhs_in)),
        /* init list */ assembleCase(NO_CONSTRAINTS),
        //        /* init list */ maxAbsValue(0.0),
        //        /* init list */ gSize(lhs.gSize)
        /* init list */ gSize(TrialBase<TrialFunctionType>::gSize())
        
        {
            std::cout<<greenColor<<"Creating WeakProblem "<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void assemble()
        {
            assembleCase=NO_CONSTRAINTS;
            
            // Get global triplets
//            maxAbsValue=0.0;
            std::vector<Eigen::Triplet<double> > globalTriplets(lhs().globalTriplets());
            
            // Create sparse matrix from global triplets
            std::cout<<"Creating matrix from triplets..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            A.resize(gSize,gSize);
            A.setFromTriplets(globalTriplets.begin(),globalTriplets.end());
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            
            // Assemble b
            b=rhs().globalVector();
            
            // Resize x
            x.setZero(gSize);
        }
        
        /**********************************************************************/
        void assembleWithMasterSlaveConstraints()
        {
            assemble();
            assembleCase=MASTERSLAVE;
        }
        
        
        /**********************************************************************/
        void assembleWithPenaltyConstraints(const double& pf)
        {/*!@param[in] pf the penalty multiplicative factor
          */
            assemble(); // this sets assembleCase=NO_CONSTRAINTS;
            assembleCase=PENALTY;
            
            std::cout<<"Assembling penalty constraints..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            const double K(pf*A.norm()/A.nonZeros());

            for (const auto& cIter : TrialBase<TrialFunctionType>::dirichletConditions())
            {
                const size_t& i(cIter.first);
                const double& val(cIter.second);
                A.coeffRef(i,i) += K;
                b(i) += K*val;
                x(i)=val;
            }
            
            A.makeCompressed();
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            
        }
        
        /**********************************************************************/
        void assembleWithLagrangeConstraints()
        {
            
            assembleCase=LAGRANGE;
            
            const size_t cSize(lhs.trialExpr.trial().dirichletConditions().size());
            
            if(cSize>0)
            {
                // Get globalTripletsreserve additional space
//                maxAbsValue=0.0;
                std::vector<Eigen::Triplet<double> > globalTriplets(lhs.assembleOnDomain());
                globalTriplets.reserve(globalTriplets.size()+2*cSize);
                // Resize b and x
                b.setZero(gSize+cSize);
                b.segment(0,gSize)=rhs.globalVector();
                x.setZero(gSize+cSize);
                
                std::cout<<"Assembling Lagrange constraints..."<<std::flush;
                const auto t0= std::chrono::system_clock::now();
                //                for (int c=0;c<cSize;++c)
                //                {
                //                    globalTriplets.emplace_back(gSize+c,lhs.trialExpr.trial().dcContainer[c].first,1.0);
                //                    globalTriplets.emplace_back(lhs.trialExpr.trial().dcContainer[c].first,gSize+c,1.0);
                //                    b(gSize+c)=lhs.trialExpr.trial().dcContainer[c].second;
                //                    x(lhs.trialExpr.trial().dcContainer[c].first)=lhs.trialExpr.trial().dcContainer[c].second;
                //                }
                size_t c=0;
//                for (typename DirichletConditionContainerType::const_iterator cIter =lhs.trialExpr.trial().dirichletConditions().begin();
//                     /*                                                    */ cIter!=lhs.trialExpr.trial().dirichletConditions().end();
//                     /*                                                    */ cIter++)
//                {
                for (const auto& cIter : TrialBase<TrialFunctionType>::dirichletConditions())
                {
                    const size_t& i(cIter.first);
                    const double& val(cIter.second);
                    globalTriplets.emplace_back(gSize+c,i,1.0);
                    globalTriplets.emplace_back(i,gSize+c,1.0);
                    b(gSize+c)=val;
                    x(i)=val;
                    c++;
                }
                
                std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
                
                std::cout<<"Creating matrix from triplets..."<<std::flush;
                const auto t1= std::chrono::system_clock::now();
                // Resize A and pupulate it from tempTriplets
                A.resize(gSize+cSize,gSize+cSize);
                A.setFromTriplets(globalTriplets.begin(),globalTriplets.end());
                std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<std::endl;
                
            }
            else // no Dirichlet Boundary conditions
            {
                std::cout<<"No constraints found"<<std::endl;
                assemble();
            }
            
        }

        
        
        /**********************************************************************/
        Eigen::VectorXd solveWithGuess(const double& tol, const TrialFunctionType& trial)
        {
            return solveWithGuess(tol,trial.dofVector());
        }
        
        /**********************************************************************/
        Eigen::VectorXd solve(const double& tol)
        {
            return solveWithGuess(tol,Eigen::VectorXd::Zero(gSize));
        }
        
        /**********************************************************************/
        Eigen::VectorXd solveWithGuess(const double& tol, const Eigen::VectorXd& y)
        {
            std::cout<<"pruning matrix: "<<A.nonZeros()<<"->";
            const auto t0= std::chrono::system_clock::now();
            A.prune(A.norm()/A.nonZeros(),FLT_EPSILON);
            std::cout<<A.nonZeros()<<" non-zeros ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            
            
            
            switch (assembleCase)
            {
                    /************/
                case NO_CONSTRAINTS:
                {// A is SPsD, so use ConjugateGradients sover
                    
                    const auto t0= std::chrono::system_clock::now();
                    std::cout<<"Solving (ConjugateGradient)..."<<std::flush;
                    Eigen::ConjugateGradient<SparseMatrixType> solver(A);
                    solver.setTolerance(tol);
                    x=solver.solveWithGuess(b,x);
//                    Eigen::SimplicialLDLT<SparseMatrixType> solver(A);
//                    x=solver.solve(b);
//                    Eigen::SparseLU<SparseMatrixType> solver(A);
//                    x=solver.solve(b);

                    
                    std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
                    assert(solver.info()==Eigen::Success && "SOLVER  FAILED");
                    break;
                }
                
                case LAGRANGE:
                {// A is SPsD, so use MINRES sover
                    
                    assert(0 && "FINISH HERE. MUST CHANGE GUESS");
                    
                    const auto t0= std::chrono::system_clock::now();
                    std::cout<<"Solving (MINRES)..."<<std::flush;
                    Eigen::MINRES<SparseMatrixType> solver(A);
                    solver.setTolerance(tol);
                    x=solver.solveWithGuess(b,x);
                    std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
                    assert(solver.info()==Eigen::Success && "SOLVER  FAILED");
                    break;
                }
                    
                    /************/
                case PENALTY:
                {// A is SPD, so use ConjugateGradient sover
                    
                    assert(0 && "FINISH HERE. MUST CHANGE GUESS");
                    
                    const auto t0= std::chrono::system_clock::now();
#ifdef _MODEL_PARDISO_SOLVER_
                    std::cout<<"Solving (PardisoLLT)..."<<std::flush;
                    Eigen::PardisoLLT<SparseMatrixType> solver(A);
                    x=solver.solve(b);
#else
                    std::cout<<"Solving (ConjugateGradient)..."<<std::flush;
                    Eigen::ConjugateGradient<SparseMatrixType> solver(A);
                    solver.setTolerance(tol);
                    x=solver.solveWithGuess(b,x);
#endif
                    std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
                    assert(solver.info()==Eigen::Success && "SOLVER  FAILED");
                    break;
                }
                    
                    /************/
                case MASTERSLAVE:
                {// A1 is SPD, so use ConjugateGradient sover
                    
                    std::cout<<"Setting up master-slave constraints..."<<std::flush;
                    const auto t0= std::chrono::system_clock::now();
                    
                    const size_t cSize(TrialBase<TrialFunctionType>::dirichletConditions().size());
                    const size_t tSize = gSize-cSize;
                    
                    std::vector<Eigen::Triplet<double> > tTriplets;
                    tTriplets.reserve(tSize);
                    
                    Eigen::VectorXd g(Eigen::VectorXd::Zero(gSize));
                    Eigen::VectorXd guess(Eigen::VectorXd::Zero(tSize));
                    
                    size_t startRow=0;
                    size_t col=0;
                    
//                    for (typename DirichletConditionContainerType::const_iterator cIter =lhs.trialExpr.trial().dirichletConditions().begin();
//                         /*                                                    */ cIter!=lhs.trialExpr.trial().dirichletConditions().end();
//                         /*                                                    */ cIter++)
                    for(const auto& cIter : TrialBase<TrialFunctionType>::dirichletConditions())
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
                        //                        const double& val(cIter->second); //! ??????????????
                    }
                    for (size_t row=startRow;row!=gSize;++row)
                    {
                        tTriplets.emplace_back(row,col,1.0);
                        guess(col)=y(row);
                        col++;
                    }
                    

                    
                    SparseMatrixType T(gSize,tSize);
                    T.setFromTriplets(tTriplets.begin(),tTriplets.end());
                    
                    SparseMatrixType A1(T.transpose()*A*T);
                    Eigen::VectorXd b1(T.transpose()*(b-A*g));
                    std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
                    
                    std::cout<<"pruning master-slave matrix: "<<A1.nonZeros()<<"->";
                    const auto t3= std::chrono::system_clock::now();
                    A1.prune(A1.norm()/A1.nonZeros(),FLT_EPSILON);
                    std::cout<<A1.nonZeros()<<" non-zeros ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t3)).count()<<" sec]"<<std::endl;
                    
//                    std::ofstream file("K.txt");
//                    file<<std::setprecision(15)<<std::scientific<<A1.toDense()<<std::endl;
//                    Eigen::SparseLU<SparseMatrixType>  lu(A1);
//                    std::cout<<"lu.det="<<lu.determinant()<<std::endl;
                    
                    const auto t1= std::chrono::system_clock::now();
#ifdef _MODEL_PARDISO_SOLVER_
                    std::cout<<"Solving (PardisoLLT)..."<<std::flush;
                    Eigen::PardisoLLT<SparseMatrixType> solver(A1);
                    x=T*solver.solve(b1).eval()+g;
#else
                    std::cout<<"Solving (ConjugateGradient)..."<<std::flush;
                    Eigen::ConjugateGradient<SparseMatrixType> solver(A1);
//                    Eigen::ConjugateGradient<SparseMatrixType> solver(T.transpose()*A*T); // this gives a segmentation fault
                    solver.setTolerance(tol);
                    //                    x=T*solver.solve(T.transpose()*(b-A*g))+g;
                    x=T*solver.solveWithGuess(b1,guess).eval()+g;
//                    x=T*solver.solve(b1)+g;
                    std::cout<<" (relative error ="<<solver.error()<<", tolerance="<<solver.tolerance();
#endif
                    std::cout<<") ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<std::endl;
                    assert(solver.info()==Eigen::Success && "SOLVER  FAILED");
//                    std::cout<<std::setprecision(15)<<std::scientific<<x<<std::endl;
//                    std::cout<<std::setprecision(15)<<std::scientific<<A.toDense()<<std::endl;
//                    std::cout<<std::setprecision(15)<<std::scientific<<b<<std::endl;
                    break;
                }
                    
                default:
                    break;
            }
            
            return x.segment(0,gSize);
        }
        
    };
    
}	// close namespace
#endif
