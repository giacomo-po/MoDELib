/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDtimeIntegrator_H_
#define model_DDtimeIntegrator_H_

#include <chrono>

#include <Eigen/Dense>


namespace model
{
	
	/**********************************************************************/
	/**********************************************************************/
	template <int N>
	struct DDtimeIntegrator
    {
		
		
	};
    
    /**********************************************************************/
    /**********************************************************************/
    template <>
    struct DDtimeIntegrator<0>
    {
        static constexpr auto tag="vMax integrator";
        static double dxMax;
        static double shearWaveSpeedFraction;
        
        
        /******************************************************************/
        static void initFromFile(const std::string& fileName)
        {
            
            dxMax=TextFileParser(fileName).readScalar<double>("dxMax",true);
            assert(dxMax>0.0);
            
        }
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        static void computeNodaVelocities(DislocationNetworkType& DN)
        {
            //! Compute and store DislocaitonNode velocities
            DN.assembleAndSolve();
            
            /*! Computes the time step size \f$dt\f$ for the current simulation step,
             *  based on maximum nodal velocity \f$v_{max}\f$.
             *
             *  The time step is calculated according to:
             *	\f[
             *  dt=
             *  \begin{cases}
             *		\frac{dx}{v_{max}} & v_{max} > fc_s\\
             *      \frac{dx}{fc_s} & v_{max} \le fc_s\\
             *  \end{cases}
             *	\f]
             *  where \f$c_s\f$ is the shear velocity and \f$f=0.1\f$ is a constant.
             */
            
            //			double vmax(0.0);
            int vMaxID=-1;
            double vmax=0.0;
            int nVmean=0;
            double vmean=0.0;
            double dt_mean=0.0;
            
//            std::cout<<"computing vMax for nodes: ";
            for (const auto& nodeIter : DN.nodes())
            {
                if(//   !nodeIter.second->isBoundaryNode()
                   //&& !nodeIter.second->isConnectedToBoundaryNodes()
                   //&&
                   nodeIter.second->meshPlanes().size()<3
//                   && !nodeIter.second->isOscillating()
                   )
                {
                    const double vNorm(nodeIter.second->get_V().norm());
                    vmean +=vNorm;
                    nVmean++;
                    if (vNorm>vmax)
                    {
                        vmax=vNorm;
                        vMaxID=nodeIter.first;
                    }
                }
            }
            vmean/=nVmean;
            
            if (vmax > DN.poly.cs*shearWaveSpeedFraction)
            {
                DN.set_dt(dxMax/vmax,vmax);
            }
            else
            {
                DN.set_dt(dxMax/(DN.poly.cs*shearWaveSpeedFraction),vmax);
            }
            
            if (vmean > DN.poly.cs*shearWaveSpeedFraction)
            {
                dt_mean=dxMax/vmean;
            }
            else
            {
                dt_mean=dxMax/(DN.poly.cs*shearWaveSpeedFraction);
            }
            
            model::cout<<std::setprecision(3)<<std::scientific<<" dt="<<DN.dt;
        }
        
    };

    //template <>
    double DDtimeIntegrator<0>::dxMax=10.0;

    //template <>
    double DDtimeIntegrator<0>::shearWaveSpeedFraction=1.0e-3;
	
//    /**********************************************************************/
//    /**********************************************************************/
//    template <>
//    struct DDtimeIntegrator<1>
//    {
//        //        DislocationNetworkType& DN;
//        
//        template <typename DislocationNetworkType>
//        static double step(DislocationNetworkType& DN)
//        {
//            
//            
//            
//            for(const auto& node : this->nodes())
//            {
//                node.second->storeP();
//            }
//            
//            double tol=10.0; // error of 10b with a time step of 1000
//            double error=2.0*tol;
//            while(error>tol)
//            {
//                assembleAndSolve();
//                for(const auto& node : this->nodes())
//                {
//                    node.second->stepA1(dt);
//                }
//                updateQuadraturePoints();
//                assembleAndSolve();
//                std::set<double> normSet;
//                for(const auto& node : this->nodes())
//                {
//                    normSet.insert(node.second->stepA2(dt));
//                }
//                
//                error=*normSet.rbegin();
//                
//                model::cout<<"      dt iteration dt="<<dt<<std::endl;
//                //                std::cout<<"normSet.size()="<<normSet.size()<<std::endl;
//                //                std::cout<<"        norm="<<*normSet.rbegin()<<std::endl;
//                std::cout<<"        error="<<error<<std::endl;
//                std::cout<<"        tol="<<tol<<std::endl;
//                
//                //                dt=0.9*tol/error*dt;
//                //
//                //                if(error < tol)
//                //                {
//                //                    break;
//                //                }
//                //                else
//                //                {
//                //
//                //                }
//                //if
//                
//                // compute dt for next time step
//                if(error>0.0)
//                {
//                    
//                    dt=0.9*tol/error*dt;
//                }
//                else
//                {
//                    dt*=2.0;
//                }
//                
//                
//            }
//            
//            
//            
//            return dt;
//            
//            
//        }
//        
//        
//        
//        
//    };

	
	/*********************************************************************/
	/*********************************************************************/
} // end namespace
#endif

