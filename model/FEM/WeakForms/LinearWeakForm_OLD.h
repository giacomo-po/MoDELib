/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LinearWeakForm_H_
#define model_LinearWeakForm_H_

#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <vector>
#include <Eigen/SparseCore>

#include <model/FEM/TrialOperators/EvalExpression.h>
#include <model/FEM/TrialOperators/TestExpression.h>
#include <model/Quadrature/Quadrature.h>
#include <model/Utilities/AreSameType.h>
#include <model/Utilities/TerminalColors.h>
#include <model/FEM/WeakForms/ElementaryDomain.h> // REMOVE THIS
//#include <model/FEM/WeakForms/DomainBoundary.h>
#include <model/FEM/Domains/IntegrationDomain.h>
#include <model/FEM/WeakForms/LinearWeakExpression.h>
#include <model/FEM/WeakForms/LinearWeakSum.h>


namespace model
{
    
    template <int cols>
    struct JGNselector
    {
		template <int dim>
        static Eigen::Matrix<double,dim,1> jGN(const Eigen::Matrix<double,dim,1>& jgn)
        {
            static_assert(cols==dim,"DerivedEval MUST HAVE EXACTLY either one or dim COLUMNS");
            return jgn;
		}
    
    };

    template <>
    struct JGNselector<1>
    {
		template <int dim>
        static double jGN(const Eigen::Matrix<double,dim,1>& jgn)
        {
            return jgn.norm();
		}
        
    };
    


    
    /**************************************************************************/
	/**************************************************************************/
    template <typename DerivedTest,typename DerivedEval>
	class LinearWeakForm : public LinearWeakExpression<LinearWeakForm<DerivedTest,DerivedEval>>
    {
        
        static_assert(DerivedTest::rows==DerivedEval::rows,"YOU ARE CREATING A LinearWeakForm BETWEEN EXPRESSIONS WITH DIFFERENT NUMBER OF ROWS");
//        static_assert(DerivedEval::cols==1,"DerivedEval MUST HAVE EXACTLY ONE COLUMN");
        
    public:
        
        typedef typename DerivedTest::TrialFunctionType  TrialFunctionType;
        
        
    private:
        typedef LinearWeakForm<DerivedTest,DerivedEval> LinearWeakFormType;
        typedef typename TypeTraits< TrialFunctionType>::ElementType ElementType;
        
        constexpr static int dim=TypeTraits< TrialFunctionType>::dim;
        constexpr static int dofPerNode=TypeTraits< TrialFunctionType>::dofPerNode;
        constexpr static int dofPerElement=TypeTraits< TrialFunctionType>::dofPerElement;
        
        
        typedef Eigen::Matrix<double,dofPerElement,1> ElementVectorType;
        
        
        /**********************************************************************/
        Eigen::Matrix<double,dim+1,1> face2domainBary(const Eigen::Matrix<double,dim,1>& b1, const int& boundaryFace) const
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

        
        
    public:
        
        const DerivedTest  testExp;
        const DerivedEval  evalExp;
        const size_t gSize;
        
        
        Eigen::Matrix<double,Eigen::Dynamic,1> _globalVector; // this should be private
        
        /**********************************************************************/
        LinearWeakForm(const TestExpression<DerivedTest>& testE, const EvalExpression<DerivedEval>& evalE) :
        /* init list */ testExp(testE.trial()),
        /* init list */ evalExp(evalE.derived()),
        /* init list */ gSize(testExp.nodeSize()*dofPerNode)
        {
            std::cout<<greenColor<<"Creating LinearWeakForm "<<defaultColor<<std::endl;
            _globalVector.resize(gSize);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,Eigen::Dynamic,1>& globalVector() const
        {
            
            // INTEGRATION SHOULD BE DONE HERE!!!!!
            
            return _globalVector;
        }
        
        /**********************************************************************/
        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        const LinearWeakFormType& operator*(const ElementaryDomain<dim,qOrder,QuadratureRule>& dV)
        {
            assembleOnDomain<qOrder,QuadratureRule>();
            return *this;
        }
        

//        /**********************************************************************/
//        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
//        const LinearWeakFormType& operator*(const IntegrationDomain<dim,0,qOrder,QuadratureRule>& bnd)
//        {
//            assembleOnDomain<qOrder,QuadratureRule>(bnd);
//            return *this;
//        }
        
        /**********************************************************************/
        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        void assembleOnDomain()
        {
            const double t0(clock());
            std::cout<<"Assembling LinearWeakForm on domain..."<<std::flush;
            
            _globalVector.setZero();
            
            typedef Quadrature< TrialFunctionType::dim,qOrder,QuadratureRule> QuadratureType;
            typedef typename QuadratureType::VectorDim AbscissaType;
            
            for (int n=0;n<testExp.elementSize();++n)
            {
                // Integrate element vector
                ElementVectorType ve(ElementVectorType::Zero());
                QuadratureType::integrate(this,ve,&LinearWeakFormType::domainAssemblyKernel<AbscissaType>,testExp.element(n));
                
                //Assemble element vector into global vector
                for (int i=0;i<dofPerElement;++i)
                {
                    const size_t  nodeID_I(i/dofPerNode);
                    const size_t nodeDof_I(i%dofPerNode);
                    const size_t gI= testExp.element(n).node(nodeID_I).gID*dofPerNode+nodeDof_I;
                    _globalVector(gI) += ve(i);
                }
            }
            
            std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
        }
        
        /**********************************************************************/
		template <typename AbscissaType>
        ElementVectorType domainAssemblyKernel(const AbscissaType& a, const ElementType& ele) const
        {
            const Eigen::Matrix<double,dim+1,1> bary(BarycentricTraits<dim>::x2l(a));
            return testExp.sfm(ele,bary).transpose()*evalExp(ele,bary)*ele.absJ(bary);
		}
        
//        /**********************************************************************/
//        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
//        const LinearWeakFormType& operator*(const ElementaryDomain<dim-1,qOrder,QuadratureRule>& dA)
//        {
//            assembleOnBoundary<qOrder,QuadratureRule>();
//            return *this;
//        }

        /**********************************************************************/
        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        const LinearWeakFormType& operator*(const IntegrationDomain<dim,1,qOrder,QuadratureRule>& bnd)
        {
            assembleOnBoundary<qOrder,QuadratureRule>(bnd);
            return *this;
        }
        
        /**********************************************************************/
        template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
        void assembleOnBoundary(const IntegrationDomain<dim,1,qOrder,QuadratureRule>& bnd)
        {
            const double t0(clock());
            std::cout<<"Assembling LinearWeakForm on faces ("<<bnd.size()<<" faces) ..."<<std::flush;
            
            _globalVector.setZero();
            
            typedef Quadrature< TrialFunctionType::dim-1,qOrder,QuadratureRule> QuadratureType;
            typedef typename QuadratureType::VectorDim AbscissaType;
            
//            for (int k=0;k<bnd.size();++k)
//            {
//                const int e(bnd[k].first);  // element ID
//                const int f(bnd[k].second); //    face ID
//                ElementVectorType ve(ElementVectorType::Zero());
//                QuadratureType::integrate(this,ve,&LinearWeakFormType::boundaryAssemblyKernel<AbscissaType>,testExp.element(e),f);
//                for (int i=0;i<dofPerElement;++i)
//                {
//                    const size_t  nodeID_I(i/dofPerNode);
//                    const size_t nodeDof_I(i%dofPerNode);
//                    const size_t gI= testExp.element(e).node(nodeID_I).gID*dofPerNode+nodeDof_I;
//                    _globalVector(gI) += ve(i);
//                }
//
//            }

            std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
        }
        
        /**********************************************************************/
		template <typename AbscissaType>
        ElementVectorType boundaryAssemblyKernel(const AbscissaType& a1, const ElementType& ele, const int& boundaryFace) const
        {
            
            const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(a1));
            const Eigen::Matrix<double,dim+1,1> bary(face2domainBary(b1,boundaryFace));
//            aFile<<ele.position(bary).transpose()<<" "<<ele.jGN(bary,boundaryFace).transpose()<<"\n"<<std::endl;
            return testExp.sfm(ele,bary).transpose()*evalExp(ele,bary)*JGNselector<DerivedEval::cols>::jGN(ele.jGN(bary,boundaryFace));

		}
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    // Operators
    
    template <typename T1,typename T2>
    LinearWeakForm<T1,T2> operator, (const TestExpression<T1>& testE, const EvalExpression<T2>& evalE)
    {
        return LinearWeakForm<T1,T2>(testE,evalE);
    }
    
    template <typename T1>
    LinearWeakForm<T1,Constant<double,1,1> > operator, (const TestExpression<T1>& testE, const double& c)
    {
        return LinearWeakForm<T1,Constant<double,1,1> >(testE,make_constant(c));
    }
    
    template <typename T1, int rows, int cols>
    LinearWeakForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> > operator, (const TestExpression<T1>& testE, const Eigen::Matrix<double,rows,cols>& c)
    {
        return LinearWeakForm<T1,Constant<Eigen::Matrix<double,rows,cols>,rows,cols> >(testE,make_constant(c));
    }
    
}	// close namespace
#endif


//
//        /**********************************************************************/
//		template <typename AbscissaType>
//        ElementVectorType orientedBoundaryAssemblyKernel(const AbscissaType& a1, const ElementType& ele, const int& boundaryFace) const
//        {
//            const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(a1));
//            const Eigen::Matrix<double,dim+1,1> bary(face2domainBary(b1,boundaryFace));
//            return testExp.sfm(ele,bary).transpose()*evalExp(ele,bary)*ele.jGN(bary,boundaryFace);
//		}


///**********************************************************************/
//template<int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
//void assembleOnBoundary(const IntegrationDomain<dim,1,qOrder,QuadratureRule>& bnd)
//{
//    const double t0(clock());
//    std::cout<<"Assembling LinearWeakForm on domain boundary..."<<std::flush;
//    
//    _globalVector.setZero();
//    
//    typedef Quadrature< TrialFunctionType::dim-1,qOrder,QuadratureRule> QuadratureType;
//    typedef typename QuadratureType::VectorDim AbscissaType;
//    
//    for (int n=0;n<testExp.elementSize();++n)
//    {
//        if(testExp.element(n).isBoundaryElement())
//        {
//            const std::vector<int> boundaryFaces=testExp.element(n).boundaryFaces();
//            
//            for (int f=0;f<boundaryFaces.size();++f)
//            {
//                // Integrate element vector
//                ElementVectorType ve(ElementVectorType::Zero());
//                
//                QuadratureType::integrate(this,ve,&LinearWeakFormType::boundaryAssemblyKernel<AbscissaType>,testExp.element(n),boundaryFaces[f]);
//                
//                //Assemble element vector into global vector
//                for (int i=0;i<dofPerElement;++i)
//                {
//                    const size_t  nodeID_I(i/dofPerNode);
//                    const size_t nodeDof_I(i%dofPerNode);
//                    const size_t gI= testExp.element(n).node(nodeID_I).gID*dofPerNode+nodeDof_I;
//                    _globalVector(gI) += ve(i);
//                }
//            }
//            
//        }
//        
//    }
//    
//    std::cout<<" done.["<<(clock()-t0)/CLOCKS_PER_SEC<<" sec]"<<std::endl;
//}


