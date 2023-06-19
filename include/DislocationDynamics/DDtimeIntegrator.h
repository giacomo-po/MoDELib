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
#include <TextFileParser.h>

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
        const double dxMax;
        const double shearWaveSpeedFraction;
        const int timeIntegrationMethod;
        const double dtMax;
        
        
        /******************************************************************/
        DDtimeIntegrator(const std::string& fileName);
        
        /**********************************************************************/
         /**********************************************************************/
        template <typename DislocationNetworkType>
        double getGlideTimeIncrement(const DislocationNetworkType& DN)
        {
            double dt_temp(0.0);
            switch (timeIntegrationMethod)
            {
                case 0: //Constant time step
                {
                    dt_temp = dtMax;
                    break;
                }
                    
                case 1: //Adaptive time stepper
                {
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

                    double vmax = 0.0;
                    for (const auto &nodeIter : DN.networkNodes())
                    {
                        if (    !nodeIter.second.lock()->isBoundaryNode()
                            && nodeIter.second.lock()->glidePlanes().size() < 3
                        )
                        {
                            const double vNorm(nodeIter.second.lock()->get_V().norm());
                            if (vNorm > vmax)
                            {
                                vmax = vNorm;
                            }
                        }
                    }
                    const double Av(vmax+DN.poly.cs*std::exp(-vmax/(DN.poly.cs * shearWaveSpeedFraction)));
                    dt_temp = std::ceil(dxMax / Av); // ceil to truncate accumulation of floating point errors
                    break;
                }

                default:
                {
                    assert(0 && "Time stepper not implemented");
                }
            }
            return dt_temp;

        }
        
        /**********************************************************************/
        template <typename DislocationNetworkType>
        double getClimbTimeIncrement(const DislocationNetworkType& DN)
        {
            //! Compute and store DislocaitonNode velocities
            //            DN.assembleAndSolve(runID,straightSegmentsDeq);
            
            /*! Computes the time step size \f$dt\f$ for the current simulation step,
             *  based on maximum nodal velocity \f$v_{max}\f$.
             *
             *  The time step is calculated according to:
             *    \f[
             *  dt=
             *  \begin{cases}
             *        \frac{dx}{v_{max}} & v_{max} > fc_s\\
             *      \frac{dx}{fc_s} & v_{max} \le fc_s\\
             *  \end{cases}
             *    \f]
             *  where \f$c_s\f$ is the shear velocity and \f$f=0.1\f$ is a constant.
             */
            
            //            double vmax(0.0);
            //            int vMaxID=-1;
            double vmax=0.0;
            //            int nVmean=0;
            //            double vmean=0.0;
            //            double dt_mean=0.0;
            
            //            std::cout<<"computing vMax for nodes: ";
            for (const auto& nodeIter : DN.networkNodes())
            {
                if(//   !nodeIter.second->isBoundaryNode()
                   //&& !nodeIter.second->isConnectedToBoundaryNodes()
                   //&&
                   nodeIter.second.lock()->glidePlanes().size()<3
                   //                   && !nodeIter.second->isOscillating()
                   )
                {
                    //                    std::cout<<nodeIter.second->sID<<" ";
                    const double vNorm(nodeIter.second.lock()->get_V().norm());
                    //                    vmean +=vNorm;
                    //                    nVmean++;
                    if (vNorm>vmax)
                    {
                        vmax=vNorm;
                        //                        vMaxID=nodeIter.first;
                    }
                }
            }
            //            std::cout<<std::endl;
            //            vmean/=nVmean;
            
            //            double dt(dxMax/vmax);
            //            if(dt<dtMin)
            //            {
            //                return dtMin;
            //            }
            //            else if(dt>dtMax)
            //            {
            //                return dtMax;
            //            }
            //            else
            //            {
            //                return dt;
            //            }
            
            return dxMax/vmax;
            //vmax > DN.poly.cs*shearWaveSpeedFraction? dxMax/vmax : dxMax/(DN.poly.cs*shearWaveSpeedFraction);
            
            //            if (vmax > DN.poly.cs*shearWaveSpeedFraction)
            //            {
            //                DN.set_dt(dxMax/vmax,vmax);
            //            }
            //            else
            //            {
            //                DN.set_dt(dxMax/(DN.poly.cs*shearWaveSpeedFraction),vmax);
            //            }
            //
            //            if (vmean > DN.poly.cs*shearWaveSpeedFraction)
            //            {
            //                dt_mean=dxMax/vmean;
            //            }
            //            else
            //            {
            //                dt_mean=dxMax/(DN.poly.cs*shearWaveSpeedFraction);
            //            }
            
            //            model::cout<<std::setprecision(3)<<std::scientific<<" dt="<<DN.dt;
        }
        
    };

    //template <>
// #ifndef _MODEL_GREATWHITE_
//     double DDtimeIntegrator<0>::dxMax=10.0;
//     double DDtimeIntegrator<0>::shearWaveSpeedFraction=1.0e-3;
//     int DDtimeIntegrator<0>::timeIntegrationMethod=0;
//     double DDtimeIntegrator<0>::dt=1000.0;
// #endif

    //template <>
	
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
//                std::cout<<"      dt iteration dt="<<dt<<std::endl;
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

