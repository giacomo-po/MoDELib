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
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/IO/EigenDataReader.h>
#include <model/IO/IDreader.h>


using namespace model;

template <int dim> 
class ExternalStressFieldController
{

    static constexpr int voigtSize=dim*(dim+1)/2;
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
    Eigen::Matrix<double,1,voigtSize> MachineStiffnessRatio;    //0 is stress control; infinity is pure strain control.
    Eigen::Matrix<size_t,voigtSize,2>     voigtorder;
    Eigen::Matrix<double,voigtSize,voigtSize> stressmultimachinestiffness;
    Eigen::Matrix<double,voigtSize,voigtSize> strainmultimachinestiffness;
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
    /* init list */ MachineStiffnessRatio(Eigen::Matrix<double,1,voigtSize>::Zero()),
    /* init list */ voigtorder(Eigen::Matrix<size_t,voigtSize,2>::Zero()),
    /* init list */ stressmultimachinestiffness(Eigen::Matrix<double,voigtSize,voigtSize>::Zero()),
    /* init list */ strainmultimachinestiffness(Eigen::Matrix<double,voigtSize,voigtSize>::Zero()),
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
        std::cout<<" nu="<<nu<<std::endl;
        nu_use=nu/(1.0+nu)/2.0;
        lambda=2.0*nu/(1.0-2.0*nu);

        model::EigenDataReader EDR;
        EDR.readMatrixInFile("./loadInput.txt","ExternalStress0",ExternalStress0);
	assert((ExternalStress0-ExternalStress0.transpose()).norm()<DBL_EPSILON && "ExternalStress0 is not symmetric.");
        EDR.readMatrixInFile("./loadInput.txt","ExternalStressRate",ExternalStressRate);
        assert((ExternalStressRate-ExternalStressRate.transpose()).norm()<DBL_EPSILON && "ExternalStressRate is not symmetric.");
        EDR.readMatrixInFile("./loadInput.txt","ExternalStrain0",ExternalStrain0);
        assert((ExternalStrain0-ExternalStrain0.transpose()).norm()<DBL_EPSILON && "ExternalStrain0 is not symmetric.");
        EDR.readMatrixInFile("./loadInput.txt","ExternalStrainRate",ExternalStrainRate);
        assert((ExternalStrainRate-ExternalStrainRate.transpose()).norm()<DBL_EPSILON && "ExternalStrainRate is not symmetric.");
        EDR.readMatrixInFile("./loadInput.txt","MachineStiffnessRatio",MachineStiffnessRatio);
        EDR.readScalarInFile("./loadInput.txt","relaxSteps",relaxSteps);
  
        if (DN.use_boundary)
        {
            sample_volume=DN.mesh.volume(); 
        }
        else
        { 
            MachineStiffnessRatio.setZero();
            model::cout<<"Notice: There is no boundary, can only use pure stress control!"<<std::endl;
        }
      
       // if (std::abs(ExternalStressRate.maxCoeff())<1e-20 && std::abs(ExternalStressRate.maxCoeff())<1e-20)
       // {
       //    USEchangingExternalStress=false;   //donot calculate the update of the external Stress field.
       // }


        // initilize the default voigt order 
        size_t vnum=0;
        for (size_t i=0;i<dim;i++)
        {
			for (size_t j=i;j<dim;j++)
			{			
			     voigtorder.row(vnum)<<j-i,j;
			     vnum++;
			}
        }	
        model::cout<<" voigtorder="<<voigtorder<<std::endl;
        
        // initilize the default voigt machinestiffness
        Eigen::Matrix<double,voigtSize,voigtSize>  machinestiffness=MachineStiffnessRatio.asDiagonal();
         
        stressmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness).inverse();
        strainmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness).inverse()*machinestiffness;
       /*    for (size_t i=0;i<voigtSize;i++)
        { 
           double alpha=MachineStiffnessRatio.row(i)[0];
           if (alpha>0.999e20)
           {
               stressmultimachinestiffness.row(i)[i]=0.0;
               strainmultimachinestiffness.row(i)[i]=1.0;
           }
           else
           {
               stressmultimachinestiffness.row(i)[i]=1.0/(1.0+alpha);
               strainmultimachinestiffness.row(i)[i]=alpha/(1.0+alpha);
           }
        } */
         model::cout<<" stressmultimachinestiffness="<<stressmultimachinestiffness<<std::endl;      
         model::cout<<" strainmultimachinestiffness="<<strainmultimachinestiffness<<std::endl;   
                 

        IDreader<'F',1,200,double> vReader;
        if (vReader.isGood(0,true))
        {
            vReader.read(0,true);
            const auto iter=vReader.find(runID);
            Eigen::Matrix<double,1,200> temp(Eigen::Matrix<double,1,200>::Zero());
            if (iter!=vReader.end())
            {
                temp=Eigen::Map<Eigen::Matrix<double,1,200>>(iter->second.data());
                last_update_time=temp(0);
                
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
            MatrixType pdr(DN.plasticDistortion());
          // MatrixType pdr(MatrixType::Zero());
            plasticStrain=(pdr+pdr.transpose())*0.5/sample_volume;
	        MatrixType dstrain(ExternalStrain0+ExternalStrainRate*last_update_time-plasticStrain); 
	        MatrixType S_strain(straininducedStress(dstrain,lambda));
            MatrixType S_stress(ExternalStress0+ExternalStressRate*last_update_time);
	        ExternalStress=stressconsidermachinestiffness(S_strain,S_stress);
            model::cout<<"ExternalStressFieldControllerController: F/F_0.txt cannot be opened."<<std::endl;
        }
 
    }

    /**************************************************************************/ 
	MatrixType v2m(const Eigen::Matrix<double,voigtSize,1>& voigtvector)
    {
		 //from voigt to matrix format
		 MatrixType temp(MatrixType::Zero());
         for (size_t i=0;i<voigtSize;i++)
         {
                temp.row(voigtorder.row(i)[0])[voigtorder.row(i)[1]]=voigtvector[i];
                temp.row(voigtorder.row(i)[1])[voigtorder.row(i)[0]]=voigtvector[i];
         }  
         return (temp);     
    }
    /**************************************************************************/ 
	Eigen::Matrix<double,voigtSize,1> m2v(const MatrixType& input_matrix)
    {
		//from matrix to voigt format
        Eigen::Matrix<double,voigtSize,1> temp(Eigen::Matrix<double,voigtSize,1>::Zero());
        for (size_t i=0;i<voigtSize;i++)
        { 
             temp.row(i)[0]=input_matrix.row(voigtorder.row(i)[0])[voigtorder.row(i)[1]];
        } 
        return (temp);      
    }
    /**************************************************************************/    
    MatrixType stressconsidermachinestiffness(const MatrixType& S_strain,const MatrixType& S_stress)
    {

	/*!For stress updating, 
	  * \f[
	  * \dot{\sigma}=\frac{\alpha}{1+\alpha} C::(\dot{\epsilon}-\dot{\epsilon}^p)+ \frac{1}{1+\alpha} \dot{\sigma} 
	  * \f]
          * \alpha=MachineStiffnessRatio, which is the stiffness ratio between machine stiffness and sample stiffness.
          * More details see <Physical Review Letters 117, 155502, 2016>.
        */  
            return (v2m(strainmultimachinestiffness*m2v(S_strain)+stressmultimachinestiffness*m2v(S_stress)));

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

 
		if(DN.runningID()>=relaxSteps)
		{
		    const double deltaT = DN.get_totalTime() - last_update_time;
		    last_update_time += deltaT;
		    MatrixType PSR=DN.plasticStrainRate()/sample_volume;
		    MatrixType S_Strain(straininducedStress(ExternalStrainRate*deltaT-PSR*deltaT,lambda));
		    ExternalStress+=stressconsidermachinestiffness(S_Strain,ExternalStressRate*deltaT);
		    if (DN.use_boundary)
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
