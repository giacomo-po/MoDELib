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
	
//	/**********************************************************************/
//	/**********************************************************************/
//	template <int N>
//	struct DDtimeIntegrator
//    {
//
//
//	};
    
    /**********************************************************************/
    /**********************************************************************/
//    template <>
    struct DDtimeIntegrator
    {
        static constexpr auto tag="vMax integrator";
        const double dxMax;
        const double shearWaveSpeedFraction;
        const int timeIntegrationMethod;
        const double dtMax;
        const double dpdMax;

        
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
                    
                case 1: //Adaptive time stepper, velocity base
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
                    long int vmaxID=-1;
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
                                vmaxID=nodeIter.second.lock()->sID;
                            }
                        }
                    }
//                    const double Av(vmax+DN.poly.cs*std::exp(-vmax/(DN.poly.cs * shearWaveSpeedFraction)));
                    const double vRef(1.0);
                    const double Av(vmax+vRef*std::exp(-vmax/(vRef*shearWaveSpeedFraction)));
                    dt_temp = std::ceil(dxMax / Av); // ceil to truncate accumulation of floating point errors
                    std::cout<<"v_max="<<vmax<<", Av="<<Av<<", vmaxID="<<vmaxID<<std::endl;
//                    std::cout<<"v_max_ID="<<vmaxID<<std::endl;
                    break;
                }
                    
                case 2: //Adaptive time stepper, plastic distortion based
                {
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

                    const double pdr(DN.plasticDistortionRate().norm());
//                    const double pdrMax(1.0e-7);

                    const std::tuple<double,double,double,double> length(DN.networkLength());
                    const double glissileLength(std::get<0>(length)/DN.mesh.volume());
                    const double pdrRef(glissileLength); // rho*b*cs, but b and cs are 1 in code units
                    

                    
                    const double Av(pdr+pdrRef*std::exp(-pdr/(pdrRef * shearWaveSpeedFraction)));
                    dt_temp = std::ceil(dpdMax / Av); // ceil to truncate accumulation of floating point errors
                    
//                    std::cout<<"pdr="<<pdr<<std::endl;
//                    std::cout<<"pdrRef="<<pdrRef<<std::endl;
//                    std::cout<<"Av="<<Av<<std::endl;
//                    std::cout<<"dpdMax="<<dpdMax<<std::endl;
//                    std::cout<<"dt_temp="<<dt_temp<<std::endl;

                    
//                    std::cout<<"v_max_ID="<<vmaxID<<std::endl;
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
                   && nodeIter.second.lock()->velocityReduction()==1.0
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
                    }
                }
            }
            return dxMax/vmax;
        }
        
    };

} // end namespace
#endif

