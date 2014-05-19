/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_WeakProblem_H_
#define model_WeakProblem_H_

#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

#include <Eigen/SparseCore>


#include <model/Utilities/SequentialOutputFile.h>

#include <model/Math/MINRES.h>



namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    template <typename BilinearWeakFormType, typename LinearWeakFormType>
	class WeakProblem
    {
        
        //typedef WeakForm<T1,T2,TF> BilinearWeakFormType;
        
        const BilinearWeakFormType& lhs;
        const LinearWeakFormType&   rhs;
        
        
        enum {NO_CONSTRAINTS,LAGRANGE,PENALTY};
        int assembleCase;
        
    public:
        const size_t gSize;
        
        
        typedef Eigen::SparseMatrix<double> SparseMatrixType;
        SparseMatrixType A;
        Eigen::VectorXd  b;
        Eigen::VectorXd  x;
        
        /**********************************************************************/
        WeakProblem(const BilinearWeakFormType& lhs_in, const LinearWeakFormType& rhs_in) :
        /* init list */ lhs(lhs_in),
        /* init list */ rhs(rhs_in),
        /* init list */ assembleCase(NO_CONSTRAINTS),
        /* init list */ gSize(lhs.gSize)
        {
            
            std::cout<<greenColor<<"Creating WeakProblem "<<defaultColor<<std::endl;
            
            
            
            if(lhs.gSize!=rhs.gSize)
            {
                std::cout<<"LHS and RHS OF WEAKPROBLEM HAVE DIFERENT GLOBAL SIZE! Exiting."<<std::endl;
                std::exit(EXIT_FAILURE);
            }
            
        }
        
        
        /**********************************************************************/
        void assemble()
        {
            assembleCase=NO_CONSTRAINTS;
            
            
            const double t0(clock());
            std::cout<<"    Assembling global matrix..."<<std::flush;
            
            // Create sparse stiffness matrix from global triplets
            A.resize(gSize,gSize);
            A.setFromTriplets(lhs.globalTriplets.begin(),lhs.globalTriplets.end());
            
            //            std::cout<<" (pruning: "<<A.nonZeros()<<"->";
            //            A.prune(lhs.maxAbsValue);
            //            std::cout<<A.nonZeros()<<" non-zeros)";
            
            // Resize b
            b=rhs.globalVector;
            x.resize(gSize);
            
            std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
        }
        
        
        /**********************************************************************/
        void assembleWithPenaltyConstraints(const double& pf)
        {/*!@param[in] pf the penalty multiplicative factor
          */
            assemble(); // this sets assembleCase=NO_CONSTRAINTS;
            assembleCase=PENALTY;
            
            const double K(pf*lhs.maxAbsValue);
            const size_t cSize(lhs.trialExp.trial().dcContainer.size());
            for (int c=0;c<cSize;++c)
            {
                A.coeffRef(lhs.trialExp.trial().dcContainer[c].first,lhs.trialExp.trial().dcContainer[c].first) += K;
                b(lhs.trialExp.trial().dcContainer[c].first) += K*lhs.trialExp.trial().dcContainer[c].second;
                x(lhs.trialExp.trial().dcContainer[c].first)=lhs.trialExp.trial().dcContainer[c].second;
            }
            A.makeCompressed();
        }
        
        /**********************************************************************/
        void assembleWithLagrangeConstraints()
        {
            assembleCase=LAGRANGE;
            
            double t0(clock());
            std::cout<<"Assembling global constrained matrix..."<<std::endl;
            
            const size_t cSize(lhs.trialExp.trial().dcContainer.size());
            
            //            if(lhs.trialExp.tf.dcContainer.size()>0)
            if(cSize>0)
            {
                std::cout<<"    creating triplets..."<<std::flush;
                // Copy lhs.globalTriplets into tempTriplets
                //                std::vector<Eigen::Triplet<double> > tempTriplets(lhs.globalTriplets);
                std::vector<Eigen::Triplet<double> > tempTriplets;
                
                // Compute number of constraints, and reserve space in tempTriplets
                tempTriplets.reserve(tempTriplets.size()+2*cSize);
                
                tempTriplets=lhs.globalTriplets;
                
                // Resize b and x
                b.setZero(gSize+cSize);
                b.segment(0,gSize)=rhs.globalVector;
                x.setZero(gSize+cSize);
                
                for (int c=0;c<cSize;++c)
                {
                    //                    tempTriplets.emplace_back(gSize+c,lhs.trialExp.tf.dcContainer[c].first,1.0);
                    //                    tempTriplets.emplace_back(lhs.trialExp.tf.dcContainer[c].first,gSize+c,1.0);
                    
                    tempTriplets.emplace_back(gSize+c,lhs.trialExp.trial().dcContainer[c].first,1.0);
                    tempTriplets.emplace_back(lhs.trialExp.trial().dcContainer[c].first,gSize+c,1.0);
                    
                    
                    b(gSize+c)=lhs.trialExp.trial().dcContainer[c].second;
                    x(lhs.trialExp.trial().dcContainer[c].first)=lhs.trialExp.trial().dcContainer[c].second;
                }
                std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
                
                std::cout<<"    creating matrix..."<<std::flush;
                t0=clock();
                // Resize A and pupulate it from tempTriplets
                A.resize(gSize+cSize,gSize+cSize);
                A.setFromTriplets(tempTriplets.begin(),tempTriplets.end());
                std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
                
                //                std::cout<<"    pruning... ("<<A.nonZeros()<<"->";
                //                t0=clock();
                //                A.prune(lhs.maxAbsValue);
                //                std::cout<<A.nonZeros()<<" non-zeros)";
                //                std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
            }
            else // no Dirichlet Boundary conditions
            {
                std::cout<<" (no constraints found)"<<std::endl;
                assemble();
            }
            
        }
        
        /**********************************************************************/
        void solve(const double& tol, const bool& enforceExactConstraints=false)
        {
            std::cout<<" (pruning: "<<A.nonZeros()<<"->";
            double t0(clock());
            A.prune(lhs.maxAbsValue);
            std::cout<<A.nonZeros()<<" non-zeros)";
            std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
            
            
            std::cout<<"    Solving..."<<std::flush;
            t0=clock();
            
            switch (assembleCase)
            {
                case LAGRANGE:
                {
                    MINRES<double> mRS(A,b,x,tol);
                    x=mRS.xMR;
                    break;
                }
                    
                case PENALTY:
                {
                    const size_t cSize(lhs.trialExp.trial().dcContainer.size());

                    Eigen::ConjugateGradient<SparseMatrixType> solver(A);
                    solver.setTolerance(tol);
                    x=solver.solveWithGuess(b,x);
                    assert(solver.info()==Eigen::Success && "SOLVER  FAILED");
                    if(enforceExactConstraints)
                    {
                        for (int c=0;c<cSize;++c)
                        {
                            x(lhs.trialExp.trial().dcContainer[c].first)=lhs.trialExp.trial().dcContainer[c].second;
                        }
                    }
                    break;
                }
    
                default:
                    break;
            }
            
//                       lhs.trialExp.trial().dofContainer=x.segment(0,lhs.gSize);
            
            std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
            
            
        }
        
        
//        /**********************************************************************/
//        void output() const
//        {
//            const double t0(clock());
//            std::cout<<"    Outputing..."<<std::flush;
//            
//            
//            
//            
//            
//            
//            
//            //Eigen::VectorXd x(mRS.xMR);
//            
//            SequentialOutputFile<'F',1> fFile;
//            fFile<<b;
//            
//            SequentialOutputFile<'X',1> xFile;
//            //            xFile<<x;
//            
//            for (size_t n=0;n<lhs.trialExp.nodeSize();++n)
//            {
//                xFile<<x.template segment<BilinearWeakFormType::dofPerNode>(BilinearWeakFormType::dofPerNode*n).transpose()<<"\n";
//            }
//            
//            
//            //            SequentialOutputFile<'K',1> kFile;
//            //            kFile<<Eigen::MatrixXd(lhs.K);
//            
//            SequentialOutputFile<'P',1> pFile;
//            for (size_t n=0;n<lhs.trialExp.nodeSize();++n)
//            {
//                pFile<<lhs.trialExp.node(n).p0.transpose()<<"\n";
//            }
//            
//            
//            std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
//            
//        }
        
    };
    
    
    
    
    
}	// close namespace
#endif