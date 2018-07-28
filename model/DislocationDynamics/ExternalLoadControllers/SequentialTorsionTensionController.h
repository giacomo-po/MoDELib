/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *                       Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SequentialTorsionTensionController_H_
#define model_SequentialTorsionTensionController_H_

#include <iostream>
#include <sstream>      // std::stringstream
#include <cmath>
#include <cfloat>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/IO/EigenDataReader.h>
#include <model/IO/IDreader.h>
#include <model/IO/UniqueOutputFile.h>



namespace model
{

    template <int dim>
    class SequentialTorsionTensionController;
    
    template <int dim>
    using ExternalLoadController=SequentialTorsionTensionController<dim>;

    
    template <int dim>
    class SequentialTorsionTensionController
    {
        
//        static constexpr int voigtSize=dim*(dim+1)/2;
        static constexpr int outputCols=2;
        
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;

        const std::string inputFileName;
        double last_update_time;
        double lambda;  //unit is mu, lambda=2*v/(1-2*v)
        double sample_volume;
        int relaxSteps;
        int torsionSteps;
        double e33;
        double s33;
        double strainRate;
        double tau;
        double R;
        
        
    public:
        /**************************************************************************/
        SequentialTorsionTensionController():
        /* init list */ inputFileName("./externalLoadControl/SequentialTorsionTensionController.txt")
        /* init list */,last_update_time(0.0)
        /* init list */,lambda(1.0)
        /* init list */,sample_volume(0.0)
        /* init list */,relaxSteps(0)
        /* init list */,torsionSteps(0)
        /* init list */,e33(0.0)
        /* init list */,s33(0.0)
        /* init list */,strainRate(0.0)
        /* init list */,tau(0)
        /* init list */,R(0.0)
        {
            
        }
        
        /**************************************************************************/
        template <typename DislocationNetworkType>
        void init(DislocationNetworkType& DN)
        {
            const long int runID=DN.runningID();
            model::cout<<greenColor<<"Initializing SequentialTorsionTensionController at runID="<<runID<<defaultColor<<std::endl;


            
            
            model::EigenDataReader EDR;
            EDR.readScalarInFile(inputFileName,"use_externalStress",DN.use_externalStress);
            
            if (DN.use_externalStress)
            {
                assert(DN.use_boundary && "A boundary mesh is necessary. Set use_boundary=1; in DDinput.txt");
                
                lambda=2.0*Material<Isotropic>::nu/(1.0-2.0*Material<Isotropic>::nu);
                sample_volume=DN.mesh.volume();
                R=0.5*(DN.mesh.xMax()(0)-DN.mesh.xMin()(0));

                
                DN._userOutputColumn+=outputCols;  //put here in order for right bvp restart
            
                EDR.readScalarInFile(inputFileName,"relaxSteps",relaxSteps);
                EDR.readScalarInFile(inputFileName,"torsionSteps",torsionSteps);
                EDR.readScalarInFile(inputFileName,"strainRate",strainRate);
                EDR.readScalarInFile(inputFileName,"tau",tau);
                
                assert(torsionSteps>=relaxSteps);
                
                
                IDreader<'F',1,200,double> vReader;
                if (vReader.isGood(0,true))
                {
                    vReader.read(0,true);
                    const auto iter=vReader.find(runID);
                    //Eigen::Matrix<double,1,200> temp(Eigen::Matrix<double,1,200>::Zero());
                    if (iter!=vReader.end())
                    {
                        Eigen::Map<Eigen::Matrix<double,1,200>> temp(iter->second.data());
                        last_update_time=temp(0);
                        
                        size_t curCol=DN.userOutputColumn()-outputCols-1;
                        e33=temp(curCol+0);
                        s33=temp(curCol+1);
                        
                        model::cout<<"last_update_time="<<last_update_time<<std::endl;
                        model::cout<<"e33="<<e33<<std::endl;
                        model::cout<<"s33="<<s33<<std::endl;
                        
                    }
                    else
                    {
                        assert(0 && "SequentialTorsionTensionController::init runID not found inf F file");
                    }
                }
                else
                {
                    model::cout<<"SequentialTorsionTensionControllerController: F/F_0.txt cannot be opened. Inititializing from initial loading."<<std::endl;
                }
            }
            
        }
        
        /**************************************************************************/
        const MatrixDim externalStress(const VectorDim& x) const
        {
            
            
            const double theta=atan2 (x(1),x(0));
            const double r=sqrt(x(1)*x(1)+x(0)*x(0));
            
            const double tau_31=-tau*sin(theta);
            const double tau_32=+tau*cos(theta);
            
            MatrixDim temp(MatrixDim::Zero());
            temp(2,0)=tau_31*r/R;
            temp(2,1)=tau_32*r/R;
            
            temp(0,2)=temp(2,0);
            temp(1,2)=temp(2,1);
            temp(2,2)=s33;

            return temp;
        }
        
        /*************************************************************************/
        template <typename DislocationNetworkType>
        void update(const DislocationNetworkType& DN)
        {
            
            
            if(DN.runningID()>=torsionSteps)
            {
                const double deltaT = DN.get_totalTime() - last_update_time;
                last_update_time += deltaT;
                
                e33+=strainRate*deltaT;
                
                const MatrixDim plasticDistortionRate(DN.plasticDistortionRate()/sample_volume);
                const MatrixDim elasticDistortionRate((MatrixDim()<<0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,strainRate).finished()-plasticDistortionRate);

                // mu =1.0 by default
                s33+= (lambda*elasticDistortionRate.trace()+2.0*elasticDistortionRate(2,2))*deltaT;
            }
            
        }
        
        /**************************************************************************/
        void output(const long int& runID,
                           UniqueOutputFile<'F'>& f_file,
                    std::ofstream& F_labels,
                    int& labelCol) const
        {
            
            f_file<<" "<<e33<<" "<<s33<<" ";
            
            if(runID==0)
            {
                F_labels<<labelCol+0<<"    e33\n";
                F_labels<<labelCol+1<<"    s33\n";
                labelCol+=2;
            }
        }
        
    };
}
#endif
