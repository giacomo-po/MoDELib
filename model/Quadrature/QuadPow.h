/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_QUADPOW_H_
#define model_QUADPOW_H_

#include <assert.h>
#include <math.h>
#include <float.h>
#include <Eigen/Dense>



namespace model
{

	template<short unsigned int pOrder, short unsigned int qOrder, template<short unsigned int,short unsigned int> class QuadratureRule>
	class QuadPow
    {
	
    public:
        
		static const Eigen::Matrix<double,qOrder,pOrder+1>   uPow;
		static const Eigen::Matrix<double,qOrder,pOrder  >  duPow;
		static const Eigen::Matrix<double,qOrder,pOrder-1> dduPow;
		
        /**********************************************************************/
		static Eigen::Matrix<double,qOrder,pOrder+1> uPowFill()
        {
            const Eigen::Matrix<double,1,qOrder> abscissas=QuadratureRule<1,qOrder>::abcsissasAndWeights().template block<1,qOrder>(0,0);

            if ( abscissas.squaredNorm() < FLT_EPSILON)
            {
                std::cout<<"norm of quadrature="<<abscissas.squaredNorm()<<std::endl;
				assert(0);
			}
			

	
			Eigen::Matrix<double,qOrder,pOrder+1> temp;
			
			for (int i=0; i<qOrder; ++i)
            {
				for (int j=0; j<pOrder+1; ++j)
                {
					temp(i,j)=std::pow(abscissas(i),j);
				}
			}
			
			std::cout<< "Computing powers of abscissas:"<<std::endl;
			std::cout<<temp<<std::endl;
			
			
			return temp;
		}
		
        /**********************************************************************/
		static Eigen::Matrix<double,qOrder,pOrder> duPowFill()
        {
            
            const Eigen::Matrix<double,1,qOrder> abscissas=QuadratureRule<1,qOrder>::abcsissasAndWeights().template block<1,qOrder>(0,0);

			Eigen::Matrix<double,qOrder,pOrder> temp;
			
			for (int i=0; i<qOrder; ++i)
            {
				for (int j=1; j<pOrder+1; ++j)
                {
					temp(i,j-1)=j*std::pow(abscissas(i),j-1);
				}
			}
			std::cout<< "Computing 1-st derivative of powers of abscissas:"<<std::endl;
			std::cout<<temp<<std::endl;
			return temp;
		}
		
        /**********************************************************************/
		static Eigen::Matrix<double,qOrder,pOrder-1> dduPowFill()
        {
            
            const Eigen::Matrix<double,1,qOrder> abscissas=QuadratureRule<1,qOrder>::abcsissasAndWeights().template block<1,qOrder>(0,0);

			Eigen::Matrix<double,qOrder,pOrder-1> temp;
			
			for (int i=0; i<qOrder; ++i)
            {
				for (int j=2; j<pOrder+1; ++j)
                {
					temp(i,j-2)=j*(j-1)*std::pow(abscissas(i),j-2);
				}
			}
			std::cout<< "Computing 2-nd derivative of powers of abscissas:"<<std::endl;
			std::cout<<temp<<std::endl;
			return temp;
		}
	
	};
	
//    // Declare static data members
	
    template<short unsigned int pOrder, short unsigned int qOrder,  template <short unsigned int, short unsigned int> class QuadratureRule>
	const Eigen::Matrix<double,qOrder,pOrder+1> QuadPow<pOrder,qOrder,QuadratureRule>::uPow=QuadPow<pOrder,qOrder,QuadratureRule>::uPowFill();
	
	template<short unsigned int pOrder, short unsigned int qOrder,  template <short unsigned int, short unsigned int> class QuadratureRule>
	const Eigen::Matrix<double,qOrder,pOrder  > QuadPow<pOrder,qOrder,QuadratureRule>::duPow=QuadPow<pOrder,qOrder,QuadratureRule>::duPowFill();
	
	template<short unsigned int pOrder, short unsigned int qOrder,  template <short unsigned int, short unsigned int> class QuadratureRule>
	const Eigen::Matrix<double,qOrder,pOrder-1> QuadPow<pOrder,qOrder,QuadratureRule>::dduPow=QuadPow<pOrder,qOrder,QuadratureRule>::dduPowFill();
	
} // namespace model
#endif
