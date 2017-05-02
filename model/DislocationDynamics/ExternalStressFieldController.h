/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *                       Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ExternalStressFieldController_H_
#define model_ExternalStressFieldController_H_

#include <iostream>
#include <sstream>      // std::stringstream
#include <cmath>
#include <cfloat>
#include <model/DislocationDynamics/Materials/CrystalOrientation.h>
#include <model/Utilities/EigenDataReader.h>


using namespace model;

template <int dim> 
class ExternalStressFieldController
{
    typedef Eigen::Matrix<double,dim,dim> MatrixType;
    MatrixType ExternalStress; 

    //External stress control parameter
    MatrixType ExternalStress0; 
    MatrixType ExternalStressRate; 

    MatrixType ExternalStrain; 
    //External strain control parameter
    MatrixType ExternalStrain0; 
    MatrixType ExternalStrainRate;
    MatrixType plasticStrain;
    //finite machine stiffness effect
    double MachineStiffnessRatio;    //0 is stress control; infinity is pure strain control.
    double last_update_time;
    double lambda;  //unit is mu, lambda=2*v/(1-2*v)
    // Number of DD steps before applying strain.
    double nu_use;
    double sample_volume;
    int relaxSteps;


    public:   
    /**************************************************************************/
    ExternalStressFieldController():
    /* init list */ ExternalStress(MatrixType::Zero()),
    /* init list */ ExternalStress0(MatrixType::Zero()),
    /* init list */ ExternalStressRate(MatrixType::Zero()),
    /* init list */ ExternalStrain(MatrixType::Zero()),
    /* init list */ ExternalStrain0(MatrixType::Zero()),
    /* init list */ ExternalStrainRate(MatrixType::Zero()),
    /* init list */ plasticStrain(MatrixType::Zero()),
    /* init list */ MachineStiffnessRatio(0.0),
    /* init list */ last_update_time(0.0),
    /* init list */ lambda(1.0),
    /* init list */ nu_use(0.12),
    /* init list */ sample_volume(1.0),
    /* init list */ relaxSteps(0)
    {

    }
    
    /**************************************************************************/
    template <typename DislocationNetworkType>
    void init(const DislocationNetworkType& DN)
    {
        const long int runID=DN.runningID(); 
        const unsigned int userOutputColumn=DN.userOutputColumn(); 
        model::cout<<greenColor<<"Initializing External Stress Field Controller at runID="<<runID<<defaultColor<<std::endl; 


        double nu=Material<Isotropic>::nu;
        nu_use=nu/(1.0+nu)/2.0;
        lambda=2.0*nu/(1.0-2.0*nu);
  
        model::EigenDataReader EDR;
        EDR.readMatrixInFile("./loadInput.txt","ExternalStress0",ExternalStress0);
	assert((ExternalStress0-ExternalStress0.transpose()).norm()<DBL_EPSILON && "ExternalStress0 is not symmetric.");
        EDR.readMatrixInFile("./loadInput.txt","ExternalStressRate",ExternalStressRate);
        assert((ExternalStressRate-ExternalStressRate.transpose()).norm()<DBL_EPSILON && "ExternalStressRate is not symmetric.");
        EDR.readMatrixInFile("./loadInput.txt","ExternalStrain0",ExternalStressFieldController::ExternalStrain0);
        assert((ExternalStrain0-ExternalStrain0.transpose()).norm()<DBL_EPSILON && "ExternalStrain0 is not symmetric.");
        EDR.readMatrixInFile("./loadInput.txt","ExternalStrainRate",ExternalStrainRate);
        assert((ExternalStrainRate-ExternalStrainRate.transpose()).norm()<DBL_EPSILON && "ExternalStrainRate is not symmetric.");
        EDR.readScalarInFile("./loadInput.txt","MachineStiffnessRatio",ExternalStressFieldController::MachineStiffnessRatio);
        EDR.readScalarInFile("./loadInput.txt","relaxSteps",ExternalStressFieldController::relaxSteps);
  
        if (DN.shared.use_boundary)
        {
            sample_volume=DN.shared.mesh.volume(); 
        }
        else
        { 
            MachineStiffnessRatio=0.0;
            model::cout<<"Notice: There is no boundary, can only use pure stress control!"<<std::endl;
        }
      
       // if (std::abs(ExternalStressRate.maxCoeff())<1e-20 && std::abs(ExternalStressRate.maxCoeff())<1e-20)
       /// {
       //    USEchangingExternalStress=false;   //donot calculate the update of the external Stress field.
       // }
        

        VertexReader<'F',201,double> vReader;
        if (vReader.isGood(0,true))
        {
            vReader.read(0,true);
            const auto iter=vReader.find(runID);
            if (iter!=vReader.end())
            {
                Eigen::Matrix<double,1,200> temp=iter->second;
                last_update_time=iter->second(0);
                
                size_t curCol=userOutputColumn-19;
                model::cout<<"last_update_time="<<last_update_time<<std::endl;
                for(int r=0;r<3;++r)
                {
                    for(int c=0;c<3;++c)
                    {
                        ExternalStrain(r,c)=temp(curCol);
                        curCol+=1;
                    }
                }  

                for(int r=0;r<3;++r)
                {
                    for(int c=0;c<3;++c)
                    {
                        ExternalStress(r,c)=temp(curCol);
                        curCol+=1;
                    }
                }   
                 std::cout<<"reading External Strain=\n "<<ExternalStrain<<std::endl; 
                 std::cout<<"reading External StressField\n "<<ExternalStress<<std::endl;
                 plasticStrain=ExternalStrain-elasticstrain(ExternalStress,nu_use);
                 std::cout<<"Initial plastic Strain=\n "<<plasticStrain<<std::endl;
            }
            else
            {
                // assert(0 && "LoadController::init runID not found inf F file");
            }
        }
        else
        {
            MatrixType pdr(DN.plasticDistortion);
            plasticStrain=(pdr+pdr.transpose())*0.5/sample_volume;
	    MatrixType dstrain(ExternalStrain0+ExternalStrainRate*last_update_time-plasticStrain); 
	    MatrixType S_strain(straininducedStress(dstrain,lambda));
	    ExternalStress=(ExternalStress0+ExternalStressRate*last_update_time)/(1.0+MachineStiffnessRatio)+S_strain*MachineStiffnessRatio/(1.0+MachineStiffnessRatio); 
            model::cout<<"ExternalStressFieldControllerController: F/F_0.txt cannot be opened."<<std::endl;
        }
 
    }
    
    /**************************************************************************/
    MatrixType straininducedStress(const MatrixType& dstrain, const double& lambda) const
    {
	/*!for isotropic linear case, 
	  * \f[
	  * \sigma_{ij}=\lambda \epsilon_{kk} \delta_{ij}+ \mu (\epsilon_{ij}+\epsilon_{ji})
	  * \f]
        */
        return lambda*dstrain.trace()*MatrixType::Identity()+2.0*dstrain;     
    }
    /**************************************************************************/    
    MatrixType elasticstrain(const MatrixType& _externalStress,const double& nu_use) const
    {
	/*!for isotropic linear case, 
	  * \f[
	  * \epsilon_{ij}=\frac{1}{E}((1+\nu)\sigma_{ij}-\nu \sigma_{kk}\delta_{ij})
	  * \f]
        */
	return -nu_use*_externalStress.trace()*MatrixType::Identity()+_externalStress*0.5;
    } 
    /**************************************************************************/    
    const MatrixType& externalStress() const
    {
	return ExternalStress;
    }

    
    /*************************************************************************/    
    template <typename DislocationNetworkType>
    void update(const DislocationNetworkType& DN)
    {
	/*!For stress updating, 
	  * \f[
	  * \dot{\sigma}=\frac{\alpha}{1+\alpha} C::(\dot{\epsilon}-\dot{\epsilon}^p)+ \frac{1}{1+\alpha} \dot{\sigma} 
	  * \f]
          * \alpha=MachineStiffnessRatio, which is the stiffness ratio between machine stiffness and sample stiffness.
          * More details see <Physical Review Letters 117, 155502, 2016>.
        */
 
		if(DN.runningID()>=relaxSteps)
		{
		    const double deltaT = DN.get_totalTime() - last_update_time;
		    last_update_time += deltaT;
		    MatrixType PSR=DN.plasticStrainRate()/sample_volume;
		    MatrixType S_Strain(straininducedStress(ExternalStrainRate*deltaT-PSR*deltaT,lambda));
		    ExternalStress+=(ExternalStressRate*deltaT)/(1+MachineStiffnessRatio)+S_Strain*MachineStiffnessRatio/(1+MachineStiffnessRatio);
		    if (DN.shared.use_boundary)
		    {
                        plasticStrain+=PSR*deltaT;
                        ExternalStrain=elasticstrain(ExternalStress,nu_use)+plasticStrain;
                    }
                    else
                    {
                        /*The plastic strain cannot be calculated without the information of volume*/
                        ExternalStrain=elasticstrain(ExternalStress,nu_use); 
                     }
                
		}


    }
    
    /**************************************************************************/
    std::string output() const
    {

        std::stringstream os;
        os<<ExternalStrain.row(0)<<" "<<ExternalStrain.row(1)<<" "<<ExternalStrain.row(2)<<" "<<ExternalStress.row(0)<<" "<<ExternalStress.row(1)<<" "<<ExternalStress.row(2)<<" ";
        return os.str();
    }
    
};
#endif
