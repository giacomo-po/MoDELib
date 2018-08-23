/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *                       Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ClockIndentationController_H_
#define model_ClockIndentationController_H_

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
    class ClockIndentationController;
    
    template <int dim>
    using ExternalLoadController=ClockIndentationController<dim>;

    
    template <int dim>
    class ClockIndentationController
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
        int indentationStepsPerSector;
        double sectorAngleDeg;
        double maxAngleDeg;
        double F;

        double e33;
        double s33;
        double strainRate;
//        double R;
//        int sectorID;
        double lineAngle;
        MatrixDim L2G; // local-to-global rotation matrix
        VectorDim lineDir;
        double zMax;
        VectorDim C;
        
        
    public:
        /**************************************************************************/
        ClockIndentationController():
        /* init list */ inputFileName("./externalLoadControl/ClockIndentationController.txt")
        /* init list */,last_update_time(0.0)
        /* init list */,lambda(1.0)
        /* init list */,sample_volume(0.0)
        /* init list */,relaxSteps(0)
        /* init list */,indentationStepsPerSector(0)
        /* init list */,sectorAngleDeg(0.0)
        /* init list */,maxAngleDeg(360.0)
        /* init list */,F(0)
        /* init list */,e33(0.0)
        /* init list */,s33(0.0)
        /* init list */,strainRate(0.0)
//        /* init list */,R(0.0)
//                /* init list */,sectorID(0)
        /* init list */,lineAngle(0.0)
        /* init list */,L2G(MatrixDim::Identity())
        /* init list */,lineDir(L2G*VectorDim::UnitY()) // equations are based on a line load along y direction
        /* init list */,zMax(0.0)
        /* init list */,C(VectorDim::Zero())
        {
            
        }
        
        /**************************************************************************/
        template <typename DislocationNetworkType>
        void init(DislocationNetworkType& DN)
        {
            const long int runID=DN.runningID();
            model::cout<<greenColor<<"Initializing ClockIndentationController at runID="<<runID<<defaultColor<<std::endl;


            
            
            model::EigenDataReader EDR;
            EDR.readScalarInFile(inputFileName,"use_externalStress",DN.use_externalStress);
            
            if (DN.use_externalStress)
            {
                assert(DN.use_boundary && "A boundary mesh is necessary. Set use_boundary=1; in DDinput.txt");
                
                lambda=2.0*Material<Isotropic>::nu/(1.0-2.0*Material<Isotropic>::nu);
                sample_volume=DN.mesh.volume();
//                R=0.5*(DN.mesh.xMax()(0)-DN.mesh.xMin()(0));
                zMax=DN.mesh.xMax()(2);
                C<<0.5*(DN.mesh.xMin()(0)+DN.mesh.xMax()(0)),0.5*(DN.mesh.xMin()(1)+DN.mesh.xMax()(1)),zMax; // point onf the line load
//                std::cout<<"zMax="<<zMax<<std::endl;


                
                DN._userOutputColumn+=outputCols;  //put here in order for right bvp restart
            
                EDR.readScalarInFile(inputFileName,"relaxSteps",relaxSteps);
                EDR.readScalarInFile(inputFileName,"indentationStepsPerSector",indentationStepsPerSector);
                                EDR.readScalarInFile(inputFileName,"sectorAngleDeg",sectorAngleDeg);
                EDR.readScalarInFile(inputFileName,"maxAngleDeg",maxAngleDeg);
                EDR.readScalarInFile(inputFileName,"forceDensity",F);

                
//                sectorID=DN.runningID()/indentationStepsPerSector;
//                lineAngle=sectorID*sectorAngleDeg;
//                Q=Eigen::AngleAxis(lineAngle*M_PI/180.0,VectorDim::UnitZ()).toRotationMatrix();
                
                EDR.readScalarInFile(inputFileName,"strainRate",strainRate);
                
//                assert(indentationStepsPerSector>=relaxSteps);
                
                update(DN);
                
                
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
                        assert(0 && "ClockIndentationController::init runID not found inf F file");
                    }
                }
                else
                {
                    model::cout<<"ClockIndentationControllerController: F/F_0.txt cannot be opened. Inititializing from initial loading."<<std::endl;
                }
            }
            
        }
        
        /**************************************************************************/
        const MatrixDim externalStress(const VectorDim& X) const
        {
            
            
            
            
            MatrixDim temp(MatrixDim::Zero());

            
            if(lineAngle>=maxAngleDeg)
            {
                temp(2,2)=s33;
            }
            else
            {// see Arnold Verruijt "An Introduction to Soil Dynamics", section 8.3.1
                // Apply Flamant formulas for a line load in the y direction
                
                const VectorDim L(C+(X-C).dot(lineDir)*lineDir); // closest point to X on the line of load
                const VectorDim LX(X-L); // vector from L to X
                const VectorDim lx(L2G.transpose()*LX); // same as above in a coordinate system with y axis along the line of load
                assert(std::fabs(lx(1))<FLT_EPSILON && " IN LOCAL Y COORIDNATE SHOULD BE 0");
                
                
//                const VectorDim x0(Q.transpose()*(VectorDim()<<X(0),X(1),0.0).finished());
                const double& x(lx(0));
                const double& z(lx(2));
                double r2=std::pow(x*x+z*z+1.0,2);
                const MatrixDim s0(-2.0*F/M_PI/r2*(MatrixDim()<<x*x*z,0.0,x*z*z,
                                                   /*         */0.0  ,0.0,0.0,
                                                   /*         */z*z*z,0.0,z*z*z
                                                   ).finished());
                
                temp=L2G*s0*L2G.transpose();
                
            }
            
            return temp;
        }
        
        /*************************************************************************/
        template <typename DislocationNetworkType>
        void update(const DislocationNetworkType& DN)
        {
            const double deltaT = DN.get_totalTime() - last_update_time;
            last_update_time += deltaT;

            
            if(DN.runningID()>relaxSteps)
            {
                int sectorID=(DN.runningID()-relaxSteps)/indentationStepsPerSector;
                lineAngle=sectorID*sectorAngleDeg;
                
                
                if(lineAngle>=maxAngleDeg)
                {
                    
                    lineAngle=maxAngleDeg; // stop clock intenter
                    
                    e33+=strainRate*deltaT;
                    
                    const MatrixDim plasticDistortionRate(DN.plasticDistortionRate()/sample_volume);
                    const MatrixDim elasticDistortionRate((MatrixDim()<<0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,strainRate).finished()-plasticDistortionRate);
                    
                    // mu =1.0 by default
                    s33+= (lambda*elasticDistortionRate.trace()+2.0*elasticDistortionRate(2,2))*deltaT;
                }
                
                std::cout<<"ClockIndentationController: lineAngle="<<lineAngle<<std::endl;
                L2G=Eigen::AngleAxisd(lineAngle*M_PI/180.0,VectorDim::UnitZ()).toRotationMatrix();
                lineDir=L2G*VectorDim::UnitY();
                std::cout<<"ClockIndentationController: lineDir="<<lineDir.transpose()<<std::endl;



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
