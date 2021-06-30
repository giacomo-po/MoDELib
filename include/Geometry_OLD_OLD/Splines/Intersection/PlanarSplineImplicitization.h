/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PLANARSPLINEIMPLICITIZATION_H_
#define model_PLANARSPLINEIMPLICITIZATION_H_

#include <set>
#include <map>
#include <utility> // for std::pair
#include <cfloat> // for FLT_EPSILON
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Jacobi>
//#include <GeneralizedEigenSolver.h>

namespace model
{
    /**********************************************************************************************************************/
    /* PlanarSplineImplicitization<polyDegree> class template *************************************************************/
    /**********************************************************************************************************************/
    /*! \brief Class template that determines the implicit equation
     *  det\f$\left(\mathbf{M}_xX+\mathbf{M}_yY+\mathbf{M}_c\right)=0\f$
     *  of a planar spline of arbitrary degree and computes intersections with other splines of also arbitrary degree.
     *
     *  Adapted from "Algorithms for Intersecting Parametric and Algebraic Curves", Manocha, 1994.
     */
    template<short unsigned int polyDegree>
    class PlanarSplineImplicitization
    {
        
        enum {polyCoeff=polyDegree+1};
        typedef Eigen::Matrix<double,polyDegree,polyDegree> MatrixPolyDeg;
        
        /*********************************************************************/
        MatrixPolyDeg fill_Mx() const {
            MatrixPolyDeg temp;
            temp << coeffs(1,1), coeffs(1,2), coeffs(1,3),
            /*****/ coeffs(1,2), coeffs(1,3), 0.0,
            /*****/ coeffs(1,3), 0.0,         0.0;
            return temp;
        }
        
        /*********************************************************************/
        MatrixPolyDeg fill_My() const {
            MatrixPolyDeg temp;
            temp << -coeffs(0,1), -coeffs(0,2), -coeffs(0,3),
            /*****/ -coeffs(0,2), -coeffs(0,3),  0.0,
            /*****/ -coeffs(0,3),  0.0,          0.0;
            return temp;
        }
        
        /*********************************************************************/
        MatrixPolyDeg fill_Mc() const {
            MatrixPolyDeg temp;
            temp << coeffs(0,1)*coeffs(1,0)-coeffs(0,0)*coeffs(1,1),  coeffs(0,2)*coeffs(1,0)-coeffs(0,0)*coeffs(1,2),                              coeffs(0,3)*coeffs(1,0)-coeffs(0,0)*coeffs(1,3),
            /*****/ coeffs(0,2)*coeffs(1,0)-coeffs(0,0)*coeffs(1,2),  coeffs(0,3)*coeffs(1,0)-coeffs(0,0)*coeffs(1,3)+coeffs(0,2)*coeffs(1,1)-coeffs(0,1)*coeffs(1,2),  coeffs(0,3)*coeffs(1,1)-coeffs(0,1)*coeffs(1,3),
            /*****/ coeffs(0,3)*coeffs(1,0)-coeffs(0,0)*coeffs(1,3),  coeffs(0,3)*coeffs(1,1)-coeffs(0,1)*coeffs(1,3),                              coeffs(0,3)*coeffs(1,2)-coeffs(0,2)*coeffs(1,3);
            return temp;
        }
        
        
    public:
        
        const Eigen::Matrix<double,2,polyCoeff> coeffs;
        const MatrixPolyDeg Mx;
        const MatrixPolyDeg My;
        const MatrixPolyDeg Mc;
        
        static double compTol;
        static double physTol;
        
        
        /*********************************************************************/
        PlanarSplineImplicitization(const Eigen::Matrix<double,2,polyCoeff>& coeffs_in) : coeffs(coeffs_in),
        /*                                                                           */   Mx(fill_Mx()),
        /*                                                                           */   My(fill_My()),
        /*                                                                           */   Mc(fill_Mc())
        {/*! Consider a planar spline in explicit form:
          *  \f[
          *  \left[\begin{array}{c}
          *  X\\
          *  Y
          *  \end{array}\right]=\sum_{i=0}^P
          *  \left[\begin{array}{c}
          *  a_i\\
          *  b_i\\
          *  \end{array}\right]t^i
          *  \f]
          *
          *  where P is the degree of the spline. In order to derive the corresponding implicit equation of the spline, consider the polynomials:
          *  \f[
          *  \begin{array}{l}
          *	F(t)=a_it^i-X\\
          *	G(s)=b_is^i-Y
          *	\end{array}
          *  \f]
          *  and observe that each point \f$(X,Y)\f$ on the spline makes \f$F(t_0)\f$ and \f$G(s_0)\f$ vanish for some common root
          *  \f$t_0=s_0\f$. Now consider the bivariate polynomial:
          *  \f[
          *	P(t,s)\frac{F(t)G(s)-F(s)G(t)}{t-s}=s^i\left(M_{xij}X+M_{yij}Y+M_{0ij}\right)t^j
          *  \f]
          *  If \f$t_0\f$ is a common root then \f$P(t_0,s)=0\f$  \f$\forall s\f$. This is possible only if \f$\left(M_{xij}X+M_{yij}Y+M_{0ij}\right)t_0^j=0\f$.
          *  Therefore for non-trivial solutions
          *  we need:
          *  \f[
          *  det\left(\mathbf{M}_xX+\mathbf{M}_yY+\mathbf{M}_c\right)=0
          *  \f]
          *  which is the implicit equation of the spline.
          *
          *  We show that\f$P(t,s)\frac{F(t)G(s)-F(s)G(t)}{t-s}=s^i\left(M_{xij}X+M_{yij}Y+M_{0ij}\right)t^j\f$ and give expressions for
          *  \f$M_{xij}\f$, \f$M_{yij}\f$ and \f$M_{cij}\f$. Substituting the expressions for \f$F_{t}\f$ and \f$G_{s}\f$:
          *  \f[
          *	P(t,s)=X\sum_{j=0}^Pb_j\frac{t^j-s^j}{t-s}+Y\sum_{i=0}^Pa_i\frac{s^i-t^i}{t-s}+\sum_{i=0}^P\sum_{j=0}^P\frac{a_ib_jt^is^j-a_ib_js^it^j}{t-s}
          *  \f]
          *  Noticing that some terms vanish identically:
          *  \f[
          *	P(t,s)=X\sum_{j=1}^Pb_j\frac{t^j-s^j}{t-s}+Y\sum_{i=1}^Pa_i\frac{s^i-t^i}{t-s}+a_0\sum_{j=1}^P\frac{b_js^j-b_jt^j}{t-s}
          *  +b_0\sum_{i=1}^P\frac{a_it^i-a_is^i}{t-s}+\sum_{i=1}^P\sum_{j=1}^P\frac{a_ib_jt^is^j-a_ib_js^it^j}{t-s}
          *  \f]
          *  Now, recalling that \f$t^j-s^j=(t-s)\sum_{k=0}^{j-1}t^{j-1-k}s^k\f$, we obtain, for the first term:
          *  \f[
          *	X\sum_{j=1}^Pb_j\frac{t^j-s^j}{t-s}=X\sum_{j=1}^Pb_j\sum_{k=0}^{j-1}t^{j-1-k}s^k=X\sum_{l=0}^{P-1}b_{l+1}\sum_{k=0}^{l}t^{l-k}s^k
          *  =X\sum_{k=0}^{P-1}\sum_{l=k}^{P-1}b_{l+1}t^{l-k}s^k=X\sum_{k=0}^{P-1}\sum_{m=0}^{P-1-k}b_{m+k+1}t^{m}s^k=X\sum_{k=0}^{P-1}\sum_{m=0}^{P-1}b_{m+k+1}t^{m}s^k
          *  \f]
          *  where in the last step we used \f$b_{k}=0\f$ for \f$k>N\f$. Analogously:
          *  \f[
          *	a_0\sum_{j=1}^P\frac{b_js^j-b_jt^j}{t-s}=-a_0\sum_{k=0}^{P-1}\sum_{m=0}^{P-1}b_{m+k+1}t^{m}s^k
          *  \f]
          *  \f[
          *	b_0\sum_{i=1}^P\frac{a_it^i-a_is^i}{t-s}=b_0\sum_{k=0}^{P-1}\sum_{m=0}^{P-1}a_{m+k+1}t^{m}s^k
          *  \f]
          *
          *
          *  \f[
          *	P(t,s)=s^i\left(M_{xij}X+M_{yij}Y+M_{0ij}\right)t^j
          *  \f]
          
          *  We now give explicit expression for the matrices
          *  \f[
          *	\mathbf{M}_x=
          *  \left[\begin{array}{llll}
          *	b_1&b_2&\ldots&b_n\\
          *	b_2&\ddots&b_n&0\\
          *	\vdots&b_n&0&0\\
          *	b_n&0&0&0
          *	\end{array}\right]
          *  \hspace{2cm}
          *  M_{xij}=
          *  \left\{\begin{array}{ll}
          *	b_{i+j+1}& i+j+1\le N\\
          *	0& i+j+1 > N
          *	\end{array}\right.
          *  \f]
          *  \f[
          *	\mathbf{M}_y=-
          *  \left[\begin{array}{llll}
          *	a_1&a_2&\ldots&a_n\\
          *	a_2&\ddots&a_n&0\\
          *	\vdots&a_n&0&0\\
          *	a_n&0&0&0
          *	\end{array}\right]
          *  \hspace{2cm}
          *  M_{yij}=
          *  \left\{\begin{array}{ll}
          *	-a_{i+j+1}& i+j+1\le N\\
          *	0& i+j+1 > N
          *	\end{array}\right.
          *  \f]
          *  \f[
          *	\mathbf{M}_c=-a_0
          *  \left[\begin{array}{llll}
          *	b_1&b_2&\ldots&b_n\\
          *	b_2&\ddots&b_n&0\\
          *	\vdots&b_n&0&0\\
          *	b_n&0&0&0
          *	\end{array}\right]
          *  +b_0
          *  \left[\begin{array}{llll}
          *	a_1&a_2&\ldots&a_n\\
          *	a_2&\ddots&a_n&0\\
          *	\vdots&a_n&0&0\\
          *	a_n&0&0&0
          *	\end{array}\right]
          
          *  \f]
          */
            
            
            //! PlanarSplineImplicitization: HERE WE SHOULD ASSER THAT THIS SPLINE IS NOT DEGENERATE !!!!!!
            assert(polyDegree==3 && "GENERALIZE HERE IN THE WAY WE FILL Mx My M0!!!!");
            
        }
        
        //		/*********************************************************************/
        //		std::pair<bool,double> isOnSpline(const Eigen::Matrix<double,2,1>& P, const double& tol=FLT_EPSILON) const {
        //
        //			std::pair<bool,double> temp=std::make_pair(false,-1.0);
        //			MatrixPolyDeg M=Mx*P(0)+My*P(1)+M0;
        //			if (std::fabs(M.determinant())<tol){
        //				temp.first=true;
        //				//Eigen::SelfAdjointEigenSolver<MatrixPolyDeg> eigensolver(M);
        //				//temp.second=...
        //				assert(0 && "FINISH HERE!!!!");
        //			}
        //			return temp;
        //		}
        
        /*********************************************************************/
        template<short unsigned int otherPolyDegree>
        std::set<std::pair<double,double> > intersectWith(const Eigen::Matrix<double,2,otherPolyDegree+1>& otherCoeffs) const
        {/*!\brief Computes the intersection points between the implicitized spline and another spline
          * defined by the coefficients otherCoeffs.
          * @param[in]  otherCoeffs the matrix of coefficients of the other spline
          *
          *  \textbf{Theoretical Background}. Assume that the other spline has degree Q. The 2xQ coefficient matrix of the other spline defines the curve:
          *  \f[
          *  \left[\begin{array}{c}
          *  X\\
          *  Y
          *  \end{array}\right]=\sum_{i=0}^Q
          *  \left[\begin{array}{c}
          *  c_i\\
          *  d_i\\
          *  \end{array}\right]s^i
          *  \f]
          *  Substituting in the implicit form of the original spline we obtain:
          *  \f[
          *  det\left(\mathbf{M}_x\sum_{i=0}^Qc_is^i+\mathbf{M}_y\sum_{i=0}^Qd_is^i+\mathbf{M}_c\right)=0
          *  \f]
          * or, grouping powers of s:
          *  \f[
          *  det\left[\left(c_0\mathbf{M}_x+d_0\mathbf{M}_y+\mathbf{M}_c\right)+\sum_{i=1}^Q\left(c_i\mathbf{M}_x+d_i\mathbf{M}_y\right)s^i\right]
          *  =det\left[\mathbf{M}_0+\sum_{i=1}^Q\mathbf{M}_is^i\right]=0
          *  \f]
          *  where we defined \f$\mathbf{M}_0=c_0\mathbf{M}_x+d_0\mathbf{M}_y+\mathbf{M}_c\f$ and \f$\mathbf{M}_i=c_i\mathbf{M}_x+d_i\mathbf{M}_y\f$.
          *  Intersection points are therefore the roots of the above equation. We now transform the root finding problem into and eigenvalue problem.
          *  For this observe that if \f$det\left[\mathbf{M}_0+\sum_{i=1}^Q\mathbf{M}_is^i\right]=0\f$ then
          *  \f[
          *  \left[\mathbf{M}_0+\sum_{i=1}^Q\mathbf{M}_is^i\right] \mathbf{v}_0=\mathbf{0}
          *  \f]
          *  has non-trivial solution. We now introduce the vectors \f$\mathbf{v}_{i-1}=s^{i-1}\mathbf{v}_0\f$ and obtain:
          *  \f[
          *  -\mathbf{M}_0\mathbf{v}_0 = s\sum_{i=1}^Q\mathbf{M}_i\mathbf{v}_{i-1}
          *  \f]
          *  which, together with the recursive relation for the \f$\mathbf{v}_i\f$'s, leads to the generalized eigenvalue problem:
          *  \f[
          *  \left[\begin{array}{cccc}
          *	-\mathbf{M}_0& \mathbf{0}& \ldots& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{I}& \ldots& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{0}& \ddots& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{0}& \ldots& \mathbf{I}
          *  \end{array}\right]
          *  \left[\begin{array}{c}
          *	\mathbf{v}_0\\ \mathbf{v}_1\\ \vdots\\ \mathbf{v}_{N-1}\\
          *  \end{array}\right]=s
          *  \left[\begin{array}{cccc}
          *	\mathbf{M}_1& \mathbf{M}_2& \ldots& \mathbf{M}_N\\
          *	\mathbf{I}& \mathbf{0}& \ldots& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{I}& \ddots& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{0}& \mathbf{I}& \mathbf{0}
          *  \end{array}\right]
          *  \left[\begin{array}{c}
          *	\mathbf{v}_0\\ \mathbf{v}_1\\ \vdots\\ \mathbf{v}_{N-1}\\
          *  \end{array}\right]
          *  \f]
          * So we obtained the generalized eigenvalue problem:
          *  \f[
          *  \mathbf{A}\mathbf{x}=\lambda\mathbf{B}\mathbf{x}
          *  \f]
          *  where:
          *  \f[
          *  \mathbf{A}=\left[\begin{array}{cccc}
          *	-\mathbf{M}_0& \mathbf{0}& \ldots& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{I}& \ldots& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{0}& \ddots& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{0}& \ldots& \mathbf{I}
          *  \end{array}\right]
          *  \f]
          * is a symmetric matrix (since it's a linear combination of symmetric matrices) and
          *  \f[
          *  \mathbf{B}=\left[\begin{array}{cccc}
          *	\mathbf{M}_1& \mathbf{M}_2& \ldots& \mathbf{M}_N\\
          *	\mathbf{I}& \mathbf{0}& \ldots& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{I}& \ddots& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{0}& \mathbf{I}& \mathbf{0}
          *  \end{array}\right]
          *  \f]
          *  Note that \f$det(\mathbf{B})=det(\mathbf{M}_N)\f$ and that
          *  \f[
          *  \mathbf{B}^{-1}=\left[\begin{array}{cccc}
          *	\mathbf{0}& \mathbf{I}& \ldots& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{0}& \mathbf{I}& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{0}& \mathbf{0}& \mathbf{I}\\
          *	\mathbf{M}_N^{-1}& -\mathbf{M}_N^{-1}\mathbf{M}_1& \ldots& -\mathbf{M}_N^{-1}\mathbf{M}_{N-1}
          *  \end{array}\right]
          *  \f]
          *  Therefore if $\mathbf{M}_N$ one can also solve:
          *  \f[
          *  \mathbf{C}\mathbf{x}=\lambda\mathbf{x}
          *  \f]
          *  where:
          *  \f[
          *  \mathbf{C}=\mathbf{B}^{-1}\mathbf{A}=\left[\begin{array}{cccc}
          *	\mathbf{0}& \mathbf{I}& \ldots& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{0}& \mathbf{I}& \mathbf{0}\\
          *	\mathbf{0}& \mathbf{0}& \mathbf{0}& \mathbf{I}\\
          *	-\mathbf{M}_N^{-1}\mathbf{M}_0& -\mathbf{M}_N^{-1}\mathbf{M}_1& \ldots& -\mathbf{M}_N^{-1}\mathbf{M}_{N-1}
          *  \end{array}\right]
          *  \f]
          */
            
            //			assert((!SplineDegeneracy<dim,polyDegree>::isLine(coeffs,tol)) && "Spline1 is degenerate.");      // make sure that spline1 is not a line
            //			assert((!SplineDegeneracy<dim,polyDegree>::isLine(otherCoeffs,tol)) && "Spline2 is degenerate."); // make sure that spline2 is not a line
            
            //! Algorithm is:
            //! 1- Assemble the generalized eigenvalue problem At=sBt
            enum {eigenSize=polyDegree*otherPolyDegree};
            typedef Eigen::Matrix< double, eigenSize,eigenSize> EigenSizeMatrixType;
            
            
            std::set<std::pair<double,double> > intersectionParameters;

            
            const MatrixPolyDeg Mn(otherCoeffs(0,otherPolyDegree)*Mx+otherCoeffs(1,otherPolyDegree)*My);
            const Eigen::JacobiSVD<MatrixPolyDeg> jsvd(Mn, Eigen::ComputeFullU | Eigen::ComputeFullV);
            const double condNumber(jsvd.singularValues().maxCoeff()/jsvd.singularValues().minCoeff()); // singular values are positive
            
            
            //            std::cout<<"conditionNumber"<<condNumber<<std::endl;
            
            if (condNumber<1.0/FLT_EPSILON)
            {
                // Assemble the matrix C=inv(B)*A
                EigenSizeMatrixType C(Eigen::Matrix<double, eigenSize,eigenSize>::Zero());
                if (otherPolyDegree>1)
                {
                    C.template block<polyDegree*(otherPolyDegree-1),polyDegree*(otherPolyDegree-1)>(0,polyDegree).setIdentity();
                }
                C.template block<polyDegree,polyDegree>(polyDegree*(otherPolyDegree-1),0*polyDegree)= -jsvd.solve(otherCoeffs(0,0)*Mx+otherCoeffs(1,0)*My+Mc);
                for (int p=1; p<otherPolyDegree;++p)
                {
                    C.template block<polyDegree,polyDegree>(polyDegree*(otherPolyDegree-1),p*polyDegree)= -jsvd.solve(otherCoeffs(0,p)*Mx+otherCoeffs(1,p)*My);
                }
                
                //                std::cout<<"C="<<C<<std::endl;
                
                // Compute the eigenvalues and eigenvectors of C
                const Eigen::EigenSolver<EigenSizeMatrixType> es(C);
                
                //               std::cout<<std::setprecision(15)<<std::scientific<<es.eigenvalues()<<std::endl;
                //               std::cout<<std::setprecision(15)<<std::scientific<<es.eigenvectors()<<std::endl;
                //                std::cout<<es.eigenvalues()<<std::endl;
               
#ifdef _MODEL_BENCH_SPLINEINTERSECTION_

                EigenSizeMatrixType A(EigenSizeMatrixType::Identity());
                A.template block<polyDegree,polyDegree>(0,0)=-(otherCoeffs(0,0)*Mx+otherCoeffs(1,0)*My+Mc);
                
                EigenSizeMatrixType B(EigenSizeMatrixType::Zero());
                if (otherPolyDegree>1)
                {
                    B.template block<polyDegree*(otherPolyDegree-1),polyDegree*(otherPolyDegree-1)>(polyDegree,0).setIdentity();
                }
                for (int p=0; p<otherPolyDegree;++p)
                {
                    B.template block<polyDegree,polyDegree>(0,p*polyDegree)=otherCoeffs(0,p+1)*Mx+otherCoeffs(1,p+1)*My;
                }
                
                Eigen::GeneralizedEigenSolver<EigenSizeMatrixType> ges(A,B);
                //            ges.compute(A,B);
                
                std::cout<<"geav"<<ges.eigenvalues()<<std::endl<<std::endl;
                
                std::ofstream ABfile;
                ABfile.open ("AB.txt");
                //ABfile << "Writing this to a file.\n";
                ABfile<< std::scientific<<std::setprecision(15)<<A<<std::endl<<B<<std::endl;
                ABfile.close();
                
                std::ofstream coeffFile;
                coeffFile.open ("coeffFile.txt");
                //ABfile << "Writing this to a file.\n";
                coeffFile<< std::scientific<<std::setprecision(15)<<coeffs<<std::endl<<otherCoeffs<<std::endl;
                coeffFile.close();
                
                //! 2- Solve generalized eigenvalue problem At=sBt
                //            model::GeneralizedEigenSolver<double,eigenSize> GES;
                //            GES.solve(A,B);
                //
                //            //! 3- Clear and fill intersectionParameters if each real eigenvalue s in [0,1] has eigenvator ratio t=v(1)/v(0) also in [0,1]
                //            for (int i=0; i<eigenSize; ++i) {
                //                if (std::fabs(GES.D(i).imag())<tol) {					// consider only eigenvectors with 0 imaginary part
                //                    const double t=GES.V(1,i).real()/GES.V(0,i).real();
                //                    const double s=GES.D(i).real();
                //                    if ((s>=0.0 && s<=1.0) && (t>=0.0 && t<=1.0)) { // t and s must be in [0,1]
                //                        Eigen::Matrix<double,2,1> Pt(coeffs.col(0));
                //                        for(int k=1;k<polyDegree+1;++k){
                //                            Pt+=coeffs.col(k)*std::pow(t,k);
                //                        }
                //                        Eigen::Matrix<double,2,1> Ps(otherCoeffs.col(0));
                //                        for(int k=1;k<otherPolyDegree+1;++k){
                //                            Ps+=otherCoeffs.col(k)*std::pow(s,k);
                //                        }
                //                        if((Pt-Ps).norm()<physTol){
                //                            intersectionParameters.insert(std::make_pair(t,s));
                //                        }
                //                    }
                //                }
                //            }
                
                for(int i=0;i<eigenSize;++i)
                {
                    const std::complex<double> s(es.eigenvalues()(i));
                    const std::complex<double> t( es.eigenvectors().col(i)(1)/es.eigenvectors().col(i)(0) ); // the root for the spline defined by "coeffs"
                    std::cout<<t<<", "<<s<<std::endl;
                }
#endif
                
                if (es.info() == Eigen::Success)
                {
                    //                    std::cout<<"eigenvalues="<<std::endl<<std::setprecision(15)<<es.eigenvalues()<<std::endl;
                    //                    std::cout<<"eigenvectors="<<std::endl<<std::setprecision(15)<<es.eigenvectors()<<std::endl;
                    
                    // Sort real eigenvalues
                    std::map<double,const int> sMap;
                    for (int i=0; i<eigenSize; ++i)
                    {
                        if (std::fabs(es.eigenvalues()(i).imag())<compTol)
                        { // consider only eigenvectors with 0 imaginary part
                            //                            const double t( (es.eigenvectors().col(i)(1)/es.eigenvectors().col(i)(0)).real() ); // the root for the spline defined by "coeffs"
                            //                            const double s( es.eigenvalues()(i).real() ); // the root for the spline defined by "otherCoeffs"
                            //                            if ((s>=0.0 && s<=1.0) && (t>=0.0 && t<=1.0))
                            //                            if ((s>=0.0-compTol && s<=1.0+compTol) && (t>=0.0-compTol && t<=1.0+compTol))
                            //                            { // t and s must be in [0,1]
                            sMap.emplace(es.eigenvalues()(i).real(),i);
                            //                            const double s(  ); // the root for the spline defined by "otherCoeffs"
                            
                            //                                Eigen::Matrix<double,2,1> Pt(coeffs.col(0));
                            //                                for(int k=1;k<polyDegree+1;++k)
                            //                                {
                            //                                    Pt+=coeffs.col(k)*std::pow(t,k);
                            //                                }
                            //                                Eigen::Matrix<double,2,1> Ps(otherCoeffs.col(0));
                            //                                for(int k=1;k<otherPolyDegree+1;++k)
                            //                                {
                            //                                    Ps+=otherCoeffs.col(k)*std::pow(s,k);
                            //                                }
                            //                                if((Pt-Ps).norm()<physTol)
                            //                                { // TO DO: COMPARE THE squaredNorm to avoid computing the root
                            //                                    intersectionParameters.insert(std::make_pair(t,s));
                            //                                }
                            //intersectionParameters.insert(std::make_pair(t,s));
                            //                            }
                        }
                    }
                    
                    // Remove multiple eigenvalues
                    for(std::map<double,const int>::const_iterator last=sMap.begin();last!=sMap.end();)
                    {
                        std::map<double,const int>::const_iterator first = last;
                        std::map<double,const int>::const_iterator next = last;
                        next++;
                        last = sMap.upper_bound(first->first+FLT_EPSILON); // jump to the first element grater than current + FLT_EPSILON
                        //                        std::cout<<"first="<<&first<<std::endl;
                        //                        std::cout<<"next="<<&next<<std::endl;
                        //                        std::cout<<"last="<<&last<<std::endl;
                        //std::cout<<"end="<<&sMap.end()<<std::endl;
                        if(last!=next) // more than one equal elements exist
                        {
                            sMap.erase(first,last); // remove all ements between first (included) and last (excluded)
                        }
                    }
                    
                    // Check that simple eigenvlaues are in range
                    for (auto pair : sMap)
                    {
                        const double& s=pair.first;
                        const int& i=pair.second;
                        const double t( (es.eigenvectors().col(i)(1)/es.eigenvectors().col(i)(0)).real() ); // the root for the spline defined by "coeffs"
                        
                        if ((s>=0.0-compTol && s<=1.0+compTol) && (t>=0.0-compTol && t<=1.0+compTol))
                        { // t and s must be in [0,1]
                            Eigen::Matrix<double,2,1> Pt(coeffs.col(0));
                            for(int k=1;k<polyDegree+1;++k)
                            {
                                Pt+=coeffs.col(k)*std::pow(t,k);
                            }
                            Eigen::Matrix<double,2,1> Ps(otherCoeffs.col(0));
                            for(int k=1;k<otherPolyDegree+1;++k)
                            {
                                Ps+=otherCoeffs.col(k)*std::pow(s,k);
                            }
                            if((Pt-Ps).norm()<physTol)
                            { // TO DO: COMPARE THE squaredNorm to avoid computing the root
                                intersectionParameters.insert(std::make_pair(t,s));
                            }
                        }
                    }
                }
            }
            else{
                //    std::cout<<"ILL_CONDITIONED INTERSECTION"<<std::endl;
            }
            
            //! 4- Return the intersections
            return intersectionParameters;
        }
        
    };
    
    
    
    // Declare static data members;
    template<short unsigned int polyDegree>
    double PlanarSplineImplicitization<polyDegree>::compTol=FLT_EPSILON;
    
    // Declare static data members;
    template<short unsigned int polyDegree>
    double PlanarSplineImplicitization<polyDegree>::physTol=DBL_MAX;
    
    
    
    /**************************************************************************/
    /**************************************************************************/
}// namespace model
#endif




