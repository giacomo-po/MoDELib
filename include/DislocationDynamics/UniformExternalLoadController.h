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
#include <TextFileParser.h>
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
        
        MatrixDim ExternalStress;
        
        //External stress control parameter
        MatrixDim ExternalStress0;
        MatrixDim ExternalStressRate;
        
        MatrixDim ExternalStrain;
        //External strain control parameter
        MatrixDim ExternalStrain0;
        MatrixDim ExternalStrainRate;
        MatrixDim plasticStrain;
        //finite machine stiffness effect
        Eigen::Matrix<double,1,voigtSize> MachineStiffnessRatio;    //0 is stress control; infinity is pure strain control.
        Eigen::Matrix<size_t,voigtSize,2>     voigtorder;
        Eigen::Matrix<double,voigtSize,voigtSize> stressmultimachinestiffness;
        Eigen::Matrix<double,voigtSize,voigtSize> strainmultimachinestiffness;
        double last_update_time;
        double lambda;  //unit is mu, lambda=2*v/(1-2*v)
        double nu_use;
        const bool enable;
        const int relaxSteps;
        
        
    public:
        /**************************************************************************/
        template <typename DislocationNetworkType>
        UniformExternalLoadController(const DislocationNetworkType& _DN,const long int& runID) :
        //        /* init list */ this->inputFileName("./externalLoadControl/UniformExternalLoadController.txt")
        /* init list */ ExternalLoadControllerBase<DefectiveCrystalType::dim>(_DN.simulationParameters.traitsIO.simulationFolder+"/inputFiles/uniformExternalLoadController.txt")
        /* init list */,DN(_DN)
        /* init list */,ExternalStress(MatrixDim::Zero())
        /* init list */,ExternalStress0(TextFileParser(this->inputFileName).readMatrix<double>("ExternalStress0",dim,dim,true))
        /* init list */,ExternalStressRate(TextFileParser(this->inputFileName).readMatrix<double>("ExternalStressRate",dim,dim,true))
        /* init list */,ExternalStrain(MatrixDim::Zero())
        /* init list */,ExternalStrain0(TextFileParser(this->inputFileName).readMatrix<double>("ExternalStrain0",dim,dim,true))
        /* init list */,ExternalStrainRate(TextFileParser(this->inputFileName).readMatrix<double>("ExternalStrainRate",dim,dim,true))
        /* init list */,plasticStrain(MatrixDim::Zero())
        /* init list */,MachineStiffnessRatio(TextFileParser(this->inputFileName).readMatrix<double>("MachineStiffnessRatio",1,voigtSize,true))
        /* init list */,voigtorder(Eigen::Matrix<size_t,voigtSize,2>::Zero())
        /* init list */,stressmultimachinestiffness(Eigen::Matrix<double,voigtSize,voigtSize>::Zero())
        /* init list */,strainmultimachinestiffness(Eigen::Matrix<double,voigtSize,voigtSize>::Zero())
        /* init list */,last_update_time(0.0)
        /* init list */,lambda(1.0)
        /* init list */,nu_use(0.12)
        /* init list */,enable(TextFileParser(this->inputFileName).readScalar<int>("enable",true))
        /* init list */,relaxSteps(TextFileParser(this->inputFileName).readScalar<int>("relaxSteps",true))
        {
            //            const long int runID=DN.runningID();
            //            const unsigned int userOutputColumn=DN.userOutputColumn();
            std::cout<<greenColor<<"Initializing UniformExternalLoadController at runID="<<runID<<defaultColor<<std::endl;
            
            
            
            double nu=DN.poly.nu;
            std::cout<<" nu="<<nu<<std::endl;
            nu_use=nu/(1.0+nu)/2.0;
            lambda=2.0*nu/(1.0-2.0*nu);
            
            //            model::EigenDataReader EDR;
            TextFileParser parser(this->inputFileName);
            assert((ExternalStress0-ExternalStress0.transpose()).norm()<DBL_EPSILON && "ExternalStress0 is not symmetric.");
            assert((ExternalStressRate-ExternalStressRate.transpose()).norm()<DBL_EPSILON && "ExternalStressRate is not symmetric.");
            assert((ExternalStrain0-ExternalStrain0.transpose()).norm()<DBL_EPSILON && "ExternalStrain0 is not symmetric.");
            assert((ExternalStrainRate-ExternalStrainRate.transpose()).norm()<DBL_EPSILON && "ExternalStrainRate is not symmetric.");
            
            
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
            MachineStiffnessRatio.block(0,0,1,dim)=MachineStiffnessRatio.block(0,0,1,dim)*(2+lambda); // using strategy 2 in test.m file  alpha*diag(C)  (first three)*(2*mu+lambda)  (first three)*mu
            Eigen::Matrix<double,voigtSize,voigtSize>  Cinv=Eigen::Matrix<double,voigtSize,voigtSize>::Identity();
            Cinv.block(0,0,dim,dim)<<nu_use/nu, -nu_use,       -nu_use,
            -nu_use,   nu_use/nu,     -nu_use,
            -nu_use,   -nu_use,       nu_use/nu;
            Eigen::Matrix<double,voigtSize,voigtSize>  machinestiffness=MachineStiffnessRatio.asDiagonal();;
            //stressmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness).inverse();
            //strainmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness).inverse()*machinestiffness;
            stressmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness*Cinv).inverse();
            strainmultimachinestiffness=(Eigen::Matrix<double,voigtSize,voigtSize>::Identity()+machinestiffness*Cinv).inverse()*machinestiffness;
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
            std::cout<<" stressmultimachinestiffness="<<stressmultimachinestiffness<<std::endl;
            std::cout<<" strainmultimachinestiffness="<<strainmultimachinestiffness<<std::endl;
            
            
            IDreader<1,200,double> vReader(DN.simulationParameters.traitsIO.fFolder+"/F");
            if (vReader.isGood(0,true))
            {
//                vReader.readLabelsFile("F/F_labels.txt");
                vReader.readLabelsFile(DN.simulationParameters.traitsIO.fFolder+"/F_labels.txt");
                vReader.read(0,true);
                const auto iter=vReader.find(runID);
                //                    Eigen::Matrix<double,1,200> temp(Eigen::Matrix<double,1,200>::Zero());
                if (iter!=vReader.end())
                {
                    //                        temp=Eigen::Map<Eigen::Matrix<double,1,200>>(iter->second.data());
                    //                        last_update_time=temp(0);
                    last_update_time=vReader(runID,"time [b/cs]");
                    
                    //                        size_t curCol=DN.userOutputColumn()-19;
                    //                        std::cout<<"userOutputColumn="<<DN.userOutputColumn()<<std::endl;
                    //                        std::cout<<"curCol="<<curCol<<std::endl;
                    std::cout<<"last_update_time="<<last_update_time<<std::endl;
                    for(int r=0;r<dim;++r)
                    {
                        for(int c=0;c<dim;++c)
                        {
                            //                                ExternalStrain(r,c)=temp(curCol);
                            ExternalStrain(r,c)=vReader(runID,"e_"+std::to_string(r+1)+std::to_string(c+1));
                            //                              curCol+=1;
                        }
                    }
                    
                    for(int r=0;r<dim;++r)
                    {
                        for(int c=0;c<dim;++c)
                        {
                            //                                ExternalStress(r,c)=temp(curCol);
                            ExternalStress(r,c)=vReader(runID,"s_"+std::to_string(r+1)+std::to_string(c+1)+" [mu]");
                            //                                curCol+=1;
                        }
                    }
                    std::cout<<"reading External Strain=\n "<<ExternalStrain<<std::endl;
                    std::cout<<"reading External StressField\n "<<ExternalStress<<std::endl;
                    plasticStrain=ExternalStrain-elasticstrain(ExternalStress,nu_use);
                    std::cout<<"Initial plastic Strain=\n "<<plasticStrain<<std::endl;
                }
                else
                {
                    if(runID==0)
                    {// F may just have been created, but runID must be zero
                        MatrixDim pdr(DN.plasticDistortion());
                        // MatrixDim pdr(MatrixDim::Zero());
                        plasticStrain=(pdr+pdr.transpose())*0.5;
                        MatrixDim dstrain(ExternalStrain0+ExternalStrainRate*last_update_time-plasticStrain);
                        //MatrixDim S_strain(straininducedStress(dstrain,lambda));
                        MatrixDim S_stress(ExternalStress0+ExternalStressRate*last_update_time);
                        ExternalStress=stressconsidermachinestiffness(dstrain,S_stress);
                    }
                    else
                    {
                        throw std::runtime_error("UniformExternalLoadController:: runID "+std::to_string(runID)+" not found in F file.");
                    }
                    // assert(0 && "LoadController::init runID not found inf F file");
                }
            }
            else
            {
                MatrixDim pdr(DN.plasticDistortion());
                // MatrixDim pdr(MatrixDim::Zero());
                plasticStrain=(pdr+pdr.transpose())*0.5;
                MatrixDim dstrain(ExternalStrain0+ExternalStrainRate*last_update_time-plasticStrain);
                //MatrixDim S_strain(straininducedStress(dstrain,lambda));
                MatrixDim S_stress(ExternalStress0+ExternalStressRate*last_update_time);
                ExternalStress=stressconsidermachinestiffness(dstrain,S_stress);
                std::cout<<"UniformExternalLoadControllerController: F/F_0.txt cannot be opened."<<std::endl;
            }
            std::cout<<"Initial ExternalStress=\n"<<ExternalStress<<std::endl;
        }
        
        
        
        /**************************************************************************/
        MatrixDim v2m(const Eigen::Matrix<double,voigtSize,1>& voigtvector)
        {
            //from voigt to matrix format
            MatrixDim temp(MatrixDim::Zero());
            for (size_t i=0;i<voigtSize;i++)
            {
                temp.row(voigtorder.row(i)[0])[voigtorder.row(i)[1]]=voigtvector[i];
                temp.row(voigtorder.row(i)[1])[voigtorder.row(i)[0]]=voigtvector[i];
            }
            return (temp);
        }
        /**************************************************************************/
        Eigen::Matrix<double,voigtSize,1> m2v(const MatrixDim& input_matrix)
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
        MatrixDim stressconsidermachinestiffness(const MatrixDim& S_strain,const MatrixDim& S_stress)
        {/*!For stress updating,
             * \f[
             * \dot{\sigma}=\frac{\alpha}{1+\alpha} C::(\dot{\epsilon}-\dot{\epsilon}^p)+ \frac{1}{1+\alpha} \dot{\sigma}
             * \f]
             * \alpha=MachineStiffnessRatio, which is the stiffness ratio between machine stiffness and sample stiffness.
             * More details see <Physical Review Letters 117, 155502, 2016>.
             */
            return (v2m(strainmultimachinestiffness*m2v(S_strain)+stressmultimachinestiffness*m2v(S_stress)));
            
        }
        /**************************************************************************/
        MatrixDim straininducedStress(const MatrixDim& dstrain, const double& lambda) const
        {/*!for isotropic linear case,
             * \f[
             * \sigma_{ij}=\lambda \epsilon_{kk} \delta_{ij}+ \mu (\epsilon_{ij}+\epsilon_{ji})
             * \f]
             */
            return lambda*dstrain.trace()*MatrixDim::Identity()+2.0*dstrain;
        }
        /**************************************************************************/
        MatrixDim elasticstrain(const MatrixDim& _externalStress,const double& nu_use) const
        {/*!for isotropic linear case,
             * \f[
             * \epsilon_{ij}=\frac{1}{E}((1+\nu)\sigma_{ij}-\nu \sigma_{kk}\delta_{ij})
             * \f]
             */
            return -nu_use*_externalStress.trace()*MatrixDim::Identity()+_externalStress*0.5;
        }
        /**************************************************************************/
        MatrixDim stress(const VectorDim&) const override
        {
            return enable? ExternalStress : MatrixDim::Zero();
        }
        
        /**************************************************************************/
        MatrixDim strain(const VectorDim&) const override
        {
            return enable? ExternalStrain : MatrixDim::Zero();
        }
        
        /*************************************************************************/
        void update(const long int& runID) override
        {
            if(runID>=relaxSteps)
            {
                const double deltaT = DN.simulationParameters.totalTime - last_update_time;
                last_update_time += deltaT;
                MatrixDim PSR=DN.plasticStrainRate();
                ExternalStress+=stressconsidermachinestiffness(ExternalStrainRate*deltaT-PSR*deltaT,ExternalStressRate*deltaT);  //2017-12-7
                plasticStrain+=PSR*deltaT;
                ExternalStrain=elasticstrain(ExternalStress,nu_use)+plasticStrain;
            }
        }
        
        /**************************************************************************/
        void output(const long int& runID,
                    std::ofstream& f_file,
                    std::ofstream& F_labels) const override
        {
            f_file<<ExternalStrain.row(0)<<" "<<ExternalStrain.row(1)<<" "<<ExternalStrain.row(2)<<" "<<ExternalStress.row(0)<<" "<<ExternalStress.row(1)<<" "<<ExternalStress.row(2)<<" ";
            
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
