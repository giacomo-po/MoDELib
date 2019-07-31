/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_UNIFORMOPEN_H_
#define model_UNIFORMOPEN_H_

#include <cmath>
#include <assert.h>
#include <Eigen/Dense>

namespace model
{
	
	/**************************************************/
	/* UniformOpen: general case                    */
	/**************************************************/
	/*! \brief Class template defining the GaussLegendre rules for determination 
	 *	of quadrature abscissas and weights.
	 */
	template<short unsigned int dim, size_t qOrder>
	struct UniformOpen
    {
		UniformOpen()
        {
			assert(0 && "UniformOpen: dimensionality not implemented.");
		}
		
	};
	
	
	/**************************************************/
	/* UniformOpen: template specialization dim=1   */
	/**************************************************/
	/*!\brief UniformOpen<1,N> places quadrature points
     * at the barycenter of equal subdivisions of size
     * 1/N of the unit interval, that is
     * q_i = 0.5/N + i/N (1=0...N-1)
	 */
	template<size_t qOrder>
	struct UniformOpen<1,qOrder>
    {
		static Eigen::Matrix<double,2,qOrder> abcsissasAndWeights()
        {
			Eigen::Matrix<double,2,qOrder> U;
            for(size_t k=0;k<qOrder;++k)
            {
                U(0,k)= 0.5/qOrder+k*1.0/qOrder;	//1.0 forces double conversion
                U(1,k)= 1.0/qOrder;
            }
			return U;
		}
		
	};
    
    /**************************************************/
    /* UniformOpen: template specialization dim=2   */
    /**************************************************/
    /*!\brief UniformOpen<2,N> places quadrature points
     * at the barycenter N^2 self-similar subdivisions 
     * of the unit triangle
     */
    template<size_t qOrder>
    struct UniformOpen<2,qOrder>
    {
//        static_assert(0,"THIS CLASS IS NOT COMPATIBLE WITH QUADRATURE.H YET, BECAUSE THE RETURNED MATRIX HAS qOrder^2 COLS, INSTEAD OF qOrder");
        
        static Eigen::Matrix<double,3,qOrder> abcsissasAndWeights()
        {
            std::cout<<"UniformOpen<2,"<<qOrder<<"> computing abcsissas and weights"<<std::endl;
            
            const double nd=sqrt(qOrder);
            const size_t N=round(nd);
            assert(nd==N && "qOrder must be a perfect sqaure (e.g. 1,4,9,16...)");
            
            const double w=0.5/qOrder; // constant weight
            
            Eigen::Matrix<double,3,qOrder> U;
            size_t m=0;
            // Triangles pointing up
            for(size_t i=0;i<N;++i)
            {
                for(size_t j=0;j<N-i;++j)
                {
                    U(0,m)= 1.0/N*(i+1.0/3.0);	//1.0 forces double conversion
                    U(1,m)= 1.0/N*(j+1.0/3.0);	//1.0 forces double conversion
                    U(2,m)= w;
                    m++;
                }
            }
            
            // Triangles pointing down
            for(size_t i=0;i<N-1;++i)
            {
                for(size_t j=0;j<N-1-i;++j)
                {
                    U(0,m)= 1.0/N*(i+2.0/3.0);	//1.0 forces double conversion
                    U(1,m)= 1.0/N*(j+2.0/3.0);	//1.0 forces double conversion
                    U(2,m)= w;
                    m++;
                }
            }
            
            assert(m==qOrder);
            
//            std::cout<<U.transpose()<<std::endl;

            return U;
        }
        

        
    };

} // namespace model
#endif

