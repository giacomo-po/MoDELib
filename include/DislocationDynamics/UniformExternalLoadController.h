/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_UniformExternalLoadController_H_
#define model_UniformExternalLoadController_H_

#include <iostream>
#include <sstream>      // std::stringstream
#include <cmath>
#include <cfloat>

//#include <Material.h>
#include <ExternalLoadControllerBase.h>
//#include <TextFileParser.h>
#include <IDreader.h>
//#include <UniqueOutputFile.h>



namespace model
{
    
    template <typename DefectiveCrystalType>
    class UniformExternalLoadController : public ExternalLoadControllerBase<DefectiveCrystalType::dim>
    {
        
        static constexpr int dim=DefectiveCrystalType::dim;
        static constexpr int voigtSize=dim*(dim+1)/2;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
        const DefectiveCrystalType& DN;
        
//        MatrixDim ExternalStress;
        
        //External stress control parameter
//        const MatrixDim ExternalStress0;
// //       MatrixDim ExternalStressRate;
//
//  //      MatrixDim ExternalStrain;
//        //External strain control parameter
//        const MatrixDim ExternalStrain0;
//        MatrixDim ExternalStrainRate;
    //    MatrixDim plasticStrain;
        //finite machine stiffness effect
      //  Eigen::Matrix<double,1,voigtSize> this->MachineStiffnessRatio;    //0 is stress control; infinity is pure strain control.
        Eigen::Matrix<size_t,voigtSize,2>     voigtorder;
//        Eigen::Matrix<double,voigtSize,voigtSize> this->stressmultimachinestiffness;
//        Eigen::Matrix<double,voigtSize,voigtSize> this->strainmultimachinestiffness;
//        double last_update_time;
//        double this->lambda;  //unit is mu, this->lambda=2*v/(1-2*v)
//        double this->nu_use;
        const bool enable;
        const int relaxSteps;
        
        
    public:
        
        UniformExternalLoadController(const DefectiveCrystalType& _DN,const long int& runID) :
        /* init list */ ExternalLoadControllerBase<DefectiveCrystalType::dim>(_DN.simulationParameters.traitsIO.simulationFolder+"/inputFiles/uniformExternalLoadController.txt")
        /* init list */,DN(_DN)
//        /* init list */,ExternalStress0(TextFileParser(this->inputFileName).readMatrix<double>("ExternalStress0",dim,dim,true))
//        /* init list */,ExternalStrain0(TextFileParser(this->inputFileName).readMatrix<double>("ExternalStrain0",dim,dim,true))
        /* init list */,voigtorder(Eigen::Matrix<size_t,voigtSize,2>::Zero())
//        /* init list */,last_update_time(0.0)
        /* init list */,enable(TextFileParser(this->inputFileName).readScalar<int>("enable",true))
        /* init list */,relaxSteps(TextFileParser(this->inputFileName).readScalar<int>("relaxSteps",true))
        {
            std::cout<<greenColor<<"Initializing UniformExternalLoadController at runID="<<runID<<defaultColor<<std::endl;
            
            
            
            double nu=DN.poly.nu;
            std::cout<<" nu="<<nu<<std::endl;
            this->nu_use=nu/(1.0+nu)/2.0;
            this->lambda=2.0*nu/(1.0-2.0*nu);
            
//            TextFileParser parser(this->inputFileName);
//            assert((ExternalStress0-ExternalStress0.transpose()).norm()<DBL_EPSILON && "ExternalStress0 is not symmetric.");
//            assert((this->ExternalStressRate-this->ExternalStressRate.transpose()).norm()<DBL_EPSILON && "ExternalStressRate is not symmetric.");
//            assert((ExternalStrain0-ExternalStrain0.transpose()).norm()<DBL_EPSILON && "ExternalStrain0 is not symmetric.");
//            assert((this->ExternalStrainRate-this->ExternalStrainRate.transpose()).norm()<DBL_EPSILON && "ExternalStrainRate is not symmetric.");
            
            
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
            std::cout<<" voigtorder="<<voigtorder<<std::endl;
            
            // initilize the default voigt machinestiffness
            this->MachineStiffnessRatio.block(0,0,1,dim)=this->MachineStiffnessRatio.block(0,0,1,dim)*(2+this->lambda); // using strategy 2 in test.m file  alpha*diag(C)  (first three)*(2*mu+this->lambda)  (first three)*mu
            Eigen::Matrix<double,voigtSize,voigtSize>  Cinv=Eigen::Matrix<double,voigtSize,voigtSize>::Identity();
            Cinv.block(0,0,dim,dim)<<this->nu_use/nu, -this->nu_use,       -this->nu_use,
            -this->nu_use,   this->nu_use/nu,     -this->nu_use,
            -this->nu_use,   -this->nu_use,       this->nu_use/nu;
            Eigen::Matrix<double,voigtSize,voigtSize>  machinestiffness=this->MachineStiffnessRatio.asDiagonal();;
            //this->stressmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness).inverse();
            //this->strainmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness).inverse()*machinestiffness;
            this->stressmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness*Cinv).inverse();
            this->strainmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness*Cinv).inverse()*machinestiffness;
                       
            update(runID);
            
//            std::cout<<" this->stressmultimachinestiffness="<<this->stressmultimachinestiffness<<std::endl;
//            std::cout<<" this->strainmultimachinestiffness="<<this->strainmultimachinestiffness<<std::endl;
//
//
//            IDreader<1,200,double> vReader(DN.simulationParameters.traitsIO.fFolder+"/F");
//            if (vReader.isGood(0,true))
//            {// F_0.txt found
//                vReader.readLabelsFile(DN.simulationParameters.traitsIO.fFolder+"/F_labels.txt");
//                vReader.read(0,true);
//                const auto iter=vReader.find(runID);
//                if (iter!=vReader.end())
//                {// runID found in F_0.txt
//                    last_update_time=vReader(runID,"time [b/cs]");
//                    std::cout<<"last_update_time="<<last_update_time<<std::endl;
//                    for(int r=0;r<dim;++r)
//                    {
//                        for(int c=0;c<dim;++c)
//                        {
//                            this->ExternalStrain(r,c)=vReader(runID,"e_"+std::to_string(r+1)+std::to_string(c+1));
//                        }
//                    }
//
//                    for(int r=0;r<dim;++r)
//                    {
//                        for(int c=0;c<dim;++c)
//                        {
//                            this->ExternalStress(r,c)=vReader(runID,"s_"+std::to_string(r+1)+std::to_string(c+1)+" [mu]");
//                        }
//                    }
//                    std::cout<<"reading External Strain=\n "<<this->ExternalStrain<<std::endl;
//                    std::cout<<"reading External StressField\n "<<this->ExternalStress<<std::endl;
//                    this->plasticStrain=this->ExternalStrain-elasticstrain(this->ExternalStress,this->nu_use);
//                    std::cout<<"Initial plastic Strain=\n "<<this->plasticStrain<<std::endl;
//                }
//                else
//                {// runID not found in F_0.txt
//                    if(runID==0)
//                    {// F may just have been created, but runID must be zero
//                        MatrixDim pdr(DN.plasticDistortion());
//                        this->plasticStrain=(pdr+pdr.transpose())*0.5;
//                        MatrixDim dstrain(ExternalStrain0+this->ExternalStrainRate*last_update_time-this->plasticStrain);
//                        MatrixDim S_stress(ExternalStress0+this->ExternalStressRate*last_update_time);
//                        this->ExternalStress=stressconsidermachinestiffness(dstrain,S_stress);
//                        this->ExternalStrain=elasticstrain(this->ExternalStress,this->nu_use)+this->plasticStrain;
//                    }
//                    else
//                    {
//                        throw std::runtime_error("UniformExternalLoadController:: runID "+std::to_string(runID)+" not found in F file.");
//                    }
//                }
//            }
//            else
//            {// F_0.txt not found
//                std::cout<<"UniformExternalLoadControllerController: F/F_0.txt cannot be opened."<<std::endl;
//                MatrixDim pdr(DN.plasticDistortion());
//                // MatrixDim pdr(MatrixDim::Zero());
//                this->plasticStrain=(pdr+pdr.transpose())*0.5;
//                MatrixDim dstrain(ExternalStrain0+this->ExternalStrainRate*last_update_time-this->plasticStrain);
//                MatrixDim S_stress(ExternalStress0+this->ExternalStressRate*last_update_time);
//                this->ExternalStress=stressconsidermachinestiffness(dstrain,S_stress);
//                this->ExternalStrain=elasticstrain(this->ExternalStress,this->nu_use)+this->plasticStrain;
//            }
//            std::cout<<"Initial ExternalStress=\n"<<this->ExternalStress<<std::endl;
//            std::cout<<"Initial ExternalStrain=\n"<<this->ExternalStrain<<std::endl;
        }
        
        MatrixDim v2m(const Eigen::Matrix<double,voigtSize,1>& voigtvector, const bool& is_strain) const
        {
            //from voigt to matrix format
            MatrixDim temp(MatrixDim::Zero());
            for (size_t i=0;i<voigtSize;i++)
            {
                if(is_strain && voigtorder.row(i)[0] != voigtorder.row(i)[1])
                {
                    temp.row(voigtorder.row(i)[0])[voigtorder.row(i)[1]]=0.5*voigtvector[i];
                    temp.row(voigtorder.row(i)[1])[voigtorder.row(i)[0]]=0.5*voigtvector[i];
                }
                else
                {
                    temp.row(voigtorder.row(i)[0])[voigtorder.row(i)[1]]=voigtvector[i];
                    temp.row(voigtorder.row(i)[1])[voigtorder.row(i)[0]]=voigtvector[i];
                }
            }
            return temp;
        }
        
        Eigen::Matrix<double,voigtSize,1> m2v(const MatrixDim& input_matrix, const bool& is_strain) const
        {
            //from matrix to voigt format
            Eigen::Matrix<double,voigtSize,1> temp(Eigen::Matrix<double,voigtSize,1>::Zero());
            for (size_t i=0;i<voigtSize;i++)
            {
                if(is_strain && voigtorder.row(i)[0] != voigtorder.row(i)[1])
                {
                    temp.row(i)[0]=2.0*input_matrix.row(voigtorder.row(i)[0])[voigtorder.row(i)[1]];
                }
                else
                {
                    temp.row(i)[0]=input_matrix.row(voigtorder.row(i)[0])[voigtorder.row(i)[1]];
                }
            }
            return temp;
        }
        
        MatrixDim stressconsidermachinestiffness(const MatrixDim& S_strain,const MatrixDim& S_stress) const override
        {/*!For stress updating,
             * \f[
             * \dot{\sigma}=\frac{\alpha}{1+\alpha} C::(\dot{\epsilon}-\dot{\epsilon}^p)+ \frac{1}{1+\alpha} \dot{\sigma}
             * \f]
             * \alpha=this->MachineStiffnessRatio, which is the stiffness ratio between machine stiffness and sample stiffness.
             * More details see <Physical Review Letters 117, 155502, 2016>.
             */
            return v2m(this->strainmultimachinestiffness*m2v(S_strain,true)+this->stressmultimachinestiffness*m2v(S_stress,false),false);
            
        }
        
        MatrixDim straininducedStress(const MatrixDim& dstrain, const double& lambda) const
        {/*!for isotropic linear case,
             * \f[
             * \sigma_{ij}=\this->lambda \epsilon_{kk} \delta_{ij}+ \mu (\epsilon_{ij}+\epsilon_{ji})
             * \f]
             */
            return this->lambda*dstrain.trace()*MatrixDim::Identity()+2.0*dstrain;
        }
        
        MatrixDim elasticstrain(const MatrixDim& _externalStress,const double& nu_use) const override
        {/*!for isotropic linear case,
             * \f[
             * \epsilon_{ij}=\frac{1}{E}((1+\nu)\sigma_{ij}-\nu \sigma_{kk}\delta_{ij})
             * \f]
             */
            return -nu_use*_externalStress.trace()*MatrixDim::Identity()+_externalStress*0.5;
        }
        
        MatrixDim stress(const VectorDim&) const override
        {
            return enable? this->ExternalStress : MatrixDim::Zero();
        }
        
        MatrixDim strain(const VectorDim&) const override
        {
            return enable? this->ExternalStrain : MatrixDim::Zero();
        }
        
//        void update(const long int& runID) override
//        {
//            if(runID>=relaxSteps)
//            {
//                const double deltaT = DN.simulationParameters.totalTime - last_update_time;
//                last_update_time += deltaT;
//                MatrixDim PSR=DN.plasticStrainRate();
//                this->ExternalStress+=stressconsidermachinestiffness(this->ExternalStrainRate*deltaT-PSR*deltaT,this->ExternalStressRate*deltaT);  //2017-12-7
//                this->plasticStrain+=PSR*deltaT;
//                this->ExternalStrain=elasticstrain(this->ExternalStress,this->nu_use)+this->plasticStrain;
//            }
//        }
        
        void update(const long int& runID) override
        {
            if(runID>=relaxSteps)
            {
//                const double deltaT = DN.simulationParameters.totalTime - last_update_time;
//               last_update_time += deltaT;
                this->plasticStrain=DN.plasticStrain();
//                MatrixDim PSR=DN.plasticStrainRate();
                this->ExternalStress=stressconsidermachinestiffness(this->ExternalStrain0+this->ExternalStrainRate*DN.simulationParameters.totalTime-this->plasticStrain,this->ExternalStress0+this->ExternalStressRate*DN.simulationParameters.totalTime);  //2017-12-7
                this->ExternalStrain=elasticstrain(this->ExternalStress,this->nu_use)+this->plasticStrain;
            }
        }
        
        void output(const long int& runID,
                    std::ofstream& f_file,
                    std::ofstream& F_labels) const override
        {
            f_file<<this->ExternalStrain.row(0)<<" "<<this->ExternalStrain.row(1)<<" "<<this->ExternalStrain.row(2)<<" "<<this->ExternalStress.row(0)<<" "<<this->ExternalStress.row(1)<<" "<<this->ExternalStress.row(2)<<" ";
            
            if(runID==0)
            {
                F_labels<<"e_11\n";
                F_labels<<"e_12\n";
                F_labels<<"e_13\n";
                F_labels<<"e_21\n";
                F_labels<<"e_22\n";
                F_labels<<"e_23\n";
                F_labels<<"e_31\n";
                F_labels<<"e_32\n";
                F_labels<<"e_33\n";
                F_labels<<"s_11 [mu]\n";
                F_labels<<"s_12 [mu]\n";
                F_labels<<"s_13 [mu]\n";
                F_labels<<"s_21 [mu]\n";
                F_labels<<"s_22 [mu]\n";
                F_labels<<"s_23 [mu]\n";
                F_labels<<"s_31 [mu]\n";
                F_labels<<"s_32 [mu]\n";
                F_labels<<"s_33 [mu]"<<std::endl;
            }
        }
        
    };
}
#endif
