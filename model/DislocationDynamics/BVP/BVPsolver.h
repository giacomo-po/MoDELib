/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BVPsolver_H_
#define model_BVPsolver_H_

#include <array>
#include <deque>
#include <chrono>
#include <memory> // unique_ptr
#include <model/FEM/FiniteElement.h>
#include <model/DislocationDynamics/Materials/Material.h>
//#include <model/DislocationDynamics/BVP/DislocationNegativeFields.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/FEM/WeakForms/JGNselector.h>
#include <model/ParticleInteraction/SingleFieldPoint.h>
#include <model/DislocationDynamics/BVP/BoundaryDisplacementPoint.h>
#include <model/DislocationDynamics/BVP/BoundaryStressPoint.h>
#include <model/FEM/Domains/LinearWeakList.h>
#include <model/Utilities/RuntimeError.h>

namespace model
{
    
    template <int dim, int sfOrder>
    class BVPsolver
    {
        
    public:
        
        typedef LagrangeElement<dim,sfOrder> ElementType;
        typedef FiniteElement<ElementType> FiniteElementType;
        typedef TrialFunction<dim,FiniteElementType> TrialFunctionType;
        typedef typename TrialFunctionType::DirichletConditionContainerType DirichletConditionContainerType;
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
        
        
        const SimplicialMesh<dim>& mesh;
        size_t gSize;
        
    private:
        Eigen::Matrix<double,6,6> C; // matrix of elastic moduli
        FiniteElementType* fe;
        TrialFunctionType*  u;  // displacement field *u=[u1 u2 u3]'
        TrialGradType*  b;      // displacement gradient *b=[u11 u12 u13 u21 u22 u23 u31 u32 u33]'
        TrialDefType*  e;       // strain *e=[e11 e22 e33 e12 e23 e13]'
        TrialStressType* s;     // stress *s=[s11 s22 s33 s12 s23 s13]'
        
        IntegrationDomainType dV;
        BilinearWeakFormType* bWF;
        
        SparseMatrixType A;
        
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
        
        /**********************************************************************/
        Eigen::VectorXd solve(const Eigen::VectorXd& b,const Eigen::VectorXd& y) const
        {/*!@param[in] b the rhs of A*x=b
          * @param[in] y the guess for x
          */
            
            model::cout<<"Setting up null-space solver..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            
            const size_t cSize(displacement().dirichletConditions().size());
            const size_t tSize = gSize-cSize;
            
            std::vector<Eigen::Triplet<double> > tTriplets;
            tTriplets.reserve(tSize);
            
            Eigen::VectorXd g(Eigen::VectorXd::Zero(gSize));
            Eigen::VectorXd guess(Eigen::VectorXd::Zero(tSize));
            
            size_t startRow=0;
            size_t col=0;
            
            for (typename DirichletConditionContainerType::const_iterator cIter =displacement().dirichletConditions().begin();
                 /*                                                    */ cIter!=displacement().dirichletConditions().end();
                 /*                                                    */ cIter++)
            {
                const size_t& endRow = cIter->first;
                g(endRow)= cIter->second;
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
            
            SparseMatrixType T(gSize,tSize);
            T.setFromTriplets(tTriplets.begin(),tTriplets.end());
            
            SparseMatrixType A1(T.transpose()*A*T);
            Eigen::VectorXd b1(T.transpose()*(b-A*g));
            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            
            const auto t1= std::chrono::system_clock::now();
#ifdef _MODEL_PARDISO_SOLVER_
            model::cout<<"Solving (PardisoLLT)..."<<std::flush;
            Eigen::PardisoLLT<SparseMatrixType> solver(A1);
            const Eigen::VectorXd x=T*solver.solve(b1)+g;
#else
            model::cout<<"Solving (ConjugateGradient)..."<<std::flush;
            Eigen::ConjugateGradient<SparseMatrixType> solver(A1);
            //                    Eigen::ConjugateGradient<SparseMatrixType> solver(T.transpose()*A*T); // this gives a segmentation fault
            solver.setTolerance(tolerance);
            //                    x=T*solver.solve(T.transpose()*(b-A*g))+g;
            const Eigen::VectorXd x=T*solver.solveWithGuess(b1,guess)+g;
            //                    x=T*solver.solve(b1)+g;
            model::cout<<" (relative error ="<<solver.error()<<", tolerance="<<solver.tolerance();
#endif
            model::cout<<") ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<std::endl;
            assert(solver.info()==Eigen::Success && "SOLVER  FAILED");
            
            return x;
        }
        
    public:
        
        double tolerance;
        
        /**********************************************************************/
        BVPsolver(const SimplicialMesh<dim>& mesh_in) :
        /* init  */ mesh(mesh_in),
        /* init  */ gSize(0),
        /* init  */ C(get_C()),
        /* init  */ tolerance(0.0001)
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
        
        /**********************************************************************/
        const TrialStressType& stress() const
        {
            return *s;
        }
        
        /**********************************************************************/
        void init()
        {
            fe = new FiniteElementType(mesh);
            u  = new TrialFunctionType(fe->template trial<dim>());
            b  = new TrialGradType(grad(*u));
            e  = new TrialDefType(def(*u));
            C=get_C(); // Material<Isotropic>  may have changed since construction
            s  = new TrialStressType(C**e);
            dV = fe->template domain<EntireDomain,4,GaussLegendre>();
            bWF = new BilinearWeakFormType((e->test(),*s),dV);
            gSize=bWF->gSize;
            
            // Compute and store stiffness matrix
            const auto t0= std::chrono::system_clock::now();
            std::cout<<"Computing stiffness matrix: gSize="<<gSize<<std::endl;
            std::vector<Eigen::Triplet<double> > globalTriplets(bWF->assembleOnDomain());
            A.resize(gSize,gSize);
            A.setFromTriplets(globalTriplets.begin(),globalTriplets.end());
            A.prune(A.norm()/A.nonZeros(),FLT_EPSILON);
            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
        }

        /**********************************************************************/
        template<typename Condition,typename DislocationNetworkType>
        void addDirichletCondition(const size_t& nodeListID,
                                   const Condition& cond,
                                   const std::array<bool,dofPerNode>& constrainDof,
                                   const DislocationNetworkType& DN)
        {/*!@param[in] dc the Dirichlet condition
          * @param[in] the nodal dof to be constrained
          */
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"adding DirichletCondition... "<<std::flush;
            // Add the Dirichlet condition
            u->addDirichletCondition(nodeListID,cond,constrainDof);
            
            // Compute  the DislocationNetwork displacement at nodeList
            model::cout<<"subtracting DislocationDisplacement... "<<std::flush;
            typedef BoundaryDisplacementPoint<DislocationNetworkType> FieldPointType;
            typedef typename FieldPointType::DisplacementField DisplacementField;
            std::deque<FieldPointType> fieldPoints; // the container of field points
            
            for (auto node : fe->nodeList(nodeListID)) // range-based for loop (C++11)
            {
                // Compute S vector
                Eigen::Matrix<double,dim,1> s(Eigen::Matrix<double,dim,1>::Zero());
                for(auto ele : *node)
                {
                    const Eigen::Matrix<double,dim+1,1> bary(ele->simplex.pos2bary(node->P0));
                    //std::map<double,int> baryIDmap;
                    for(int k=0;k<dim+1;++k)
                    {
                        if (std::fabs(bary(k))<FLT_EPSILON && ele->simplex.child(k).isBoundarySimplex())
                        {
                            s += ele->simplex.nda.col(k).normalized();
                        }
                    }
                }
                const double sNorm(s.norm());
                RuntimeError(sNorm>FLT_EPSILON,"s-vector has zero norm.");
//                assert(sNorm>0.0 && "s-vector has zero norm.");
                fieldPoints.emplace_back(*node,s/sNorm);
            }
            DN.template computeField<FieldPointType,DisplacementField>(fieldPoints);
            
            // Subtract the DislocationNetwork displacement from the Dirichlet conditions
            for(int n=0;n<fieldPoints.size();++n)
            {
                for(int dof=0;dof<dofPerNode;++dof)
                {
                    if(constrainDof[dof])
                    {
                        u->dirichletConditions().at(fe->nodeList(nodeListID)[n]->gID*dofPerNode+dof) -= fieldPoints[n].template field<DisplacementField>()(dof);
                    }
                }
            }
            
            if(DN.shared.use_virtualSegments) // Add solid angle contribution
            {
                //std::cout<<"NEED TO COMPUTE DISPLACEMENT OF RADIAL SEGMENTS"<<std::endl;
                
                for(int n=0;n<fieldPoints.size();++n)
                {
                    VectorDim dispJump(VectorDim::Zero());
                    
                    for (typename DislocationNetworkType::NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
                    {
                        linkIter->second.addToSolidAngleJump(fieldPoints[n].P,fieldPoints[n].S,dispJump);
                    }
                    
                    for(int dof=0;dof<dofPerNode;++dof)
                    {
                        if(constrainDof[dof])
                        {
                            u->dirichletConditions().at(fe->nodeList(nodeListID)[n]->gID*dofPerNode+dof) -= dispJump(dof);
                        }
                    }
                }
            } // end virtual loops
            
            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
//        /**********************************************************************/
//        template<typename Condition,typename DislocationNetworkType>
//        void addDirichletCondition(const Condition& cond, const NodeList<FiniteElementType>& nodeList, const int& dof,
//                                   const DislocationNetworkType& DN)
//        {/*!@param[in] dc the Dirichlet condition
//          * @param[in] the nodal dof to be constrained
//          */
//            const auto t0= std::chrono::system_clock::now();
//            model::cout<<"adding DirichletCondition... "<<std::flush;
//            // Add the Dirichlet condition
//            u->addDirichletCondition(cond,nodeList,dof);
//            
//            // Compute  the DislocationNetwork displacement at nodeList
//            model::cout<<"subtracting DislocationDisplacement... "<<std::flush;
//            typedef BoundaryDisplacementPoint<DislocationNetworkType> FieldPointType;
//            typedef typename FieldPointType::DisplacementField DisplacementField;
//            std::deque<FieldPointType> fieldPoints; // the container of field points
//            
//            for (auto node : nodeList) // range-based for loop (C++11)
//            {
//                // Compute S vector
//                Eigen::Matrix<double,dim,1> s(Eigen::Matrix<double,dim,1>::Zero());
//                for(auto ele : *node)
//                {
//                    const Eigen::Matrix<double,dim+1,1> bary(ele->simplex.pos2bary(node->P0));
//                    //std::map<double,int> baryIDmap;
//                    for(int k=0;k<dim+1;++k)
//                    {
//                        if (std::fabs(bary(k))<FLT_EPSILON && ele->simplex.child(k).isBoundarySimplex())
//                        {
//                            s += ele->simplex.nda.col(k).normalized();
//                        }
//                    }
//                }
//                const double sNorm(s.norm());
//                assert(sNorm>0.0 && "s-vector has zero norm.");
//                fieldPoints.emplace_back(*node,s/sNorm);
//            }
//            DN.template computeField<FieldPointType,DisplacementField>(fieldPoints);
//            
//            // Subtract the DislocationNetwork displacement from the Dirichlet conditions
//            for(int n=0;n<fieldPoints.size();++n)
//            {
//                u->dirichletConditions().at(nodeList[n]->gID*dofPerNode+dof) -= fieldPoints[n].template field<DisplacementField>()(dof);
//            }
//            
//            if(DN.shared.use_virtualSegments) // Add solid angle contribution
//            {
//                //std::cout<<"NEED TO COMPUTE DISPLACEMENT OF RADIAL SEGMENTS"<<std::endl;
//                
//                for(int n=0;n<fieldPoints.size();++n)
//                {
//                    VectorDim dispJump(VectorDim::Zero());
//                    
//                    for (typename DislocationNetworkType::NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter)
//                    {
//                        linkIter->second.addToSolidAngleJump(fieldPoints[n].P,fieldPoints[n].S,dispJump);
//                    }
//                    
//                    u->dirichletConditions().at(nodeList[n]->gID*dofPerNode+dof) -= dispJump(dof);
//                }
//            } // end virtual loops
//            
//            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//        }
        
        
        /**********************************************************************/
        template <typename DislocationNetworkType,int qOrder>
        void assembleAndSolve(const DislocationNetworkType& DN)
        {
#ifdef userBVPfile
            
            typedef typename DislocationNetworkType::StressField StressField;
            typedef BoundaryStressPoint<DislocationNetworkType> FieldPointType;
            
            //LinearWeakFormType lwf(u->test(),ds);
            //            const auto t0= std::chrono::system_clock::now();
            
            auto ndA=fe->template boundary<ExternalBoundary,qOrder,GaussLegendre>();
            
            auto eb_list = ndA.template integrationList<FieldPointType>();
            
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"Computing DD boundary traction..."<<std::flush;
            //            if (DN.shared.use_virtualSegments)
            //            {
            //                DN.template computeField<FieldPointType,StressField>(eb_list,DN.shared.bdn);
            //            }
            //            else
            //            {
            //                DN.template computeField<FieldPointType,StressField>(eb_list);
            //            }
            DN.template computeField<FieldPointType,StressField>(eb_list);
            
            model::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            
            //
            //            if(DN.shared.use_virtualSegments)
            //            {
            //                std::cout<<"HERE NEED TO ADD STRESS OF RADIAL BOUNDARY SEGMENTS!!"<<std::endl;
            //            }
            
            
            auto dislocationTraction=(u->test(),eb_list);

            fe->clearNodeLists();
            u->clearDirichletConditions();

            
#include userBVPfile // userBVPfile defines additional loads, boundary conditions, and calls solver
#else
            assert(0 && "YOU MUST #define THE userBVPfile to use BVPsolver.");
#endif
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
        
    };
    
    
} // namespace model
#endif


