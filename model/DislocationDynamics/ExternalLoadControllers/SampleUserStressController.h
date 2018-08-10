/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *                       Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SampleUserStressController_H_
#define model_SampleUserStressController_H_

#include <iostream>
#include <sstream>      // std::stringstream
#include <cmath>
#include <cfloat>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/IO/EigenDataReader.h>
#include <model/IO/IDreader.h>


namespace model
{

    template <int dim>
    class SampleUserStressController;
    
    template <int dim>
    using ExternalLoadController=SampleUserStressController<dim>;

    template <int dim>
    class SampleUserStressController
    {
        
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
        const std::string inputFileName;
        MatrixDim _externalStress;
        VectorDim P0;
        int relaxSteps;
        double nu;  //unit is mu, lambda=2*v/(1-2*v)
        double lambda;  //unit is mu, lambda=2*v/(1-2*v)
        VectorDim meshSize;
        double tau_xz;
        size_t dislocationImages_x;
        size_t dislocationImages_y;
                size_t dislocationImages_z;
        
    public:
        
        /**************************************************************************/
        SampleUserStressController():
        /* init list */ inputFileName("./externalLoadControl/SampleUserStressController.txt")
        /* init list */,_externalStress(MatrixDim::Zero())
        /* init list */,P0(VectorDim::Zero())
        /* init list */,relaxSteps(0)
        /* init list */,nu(0.33)
        /* init list */,lambda(2.0*nu/(1.0-2.0*nu))
        /* init list */,meshSize(VectorDim::Zero())
        /* init list */,tau_xz(0.0)
        /* init list */,dislocationImages_x(0)
        /* init list */,dislocationImages_y(0)
        /* init list */,dislocationImages_z(0)
        {
            
        }
        
        /**************************************************************************/
        template <typename DislocationNetworkType>
        void init(DislocationNetworkType& DN)
        {
            const long int runID=DN.runningID();
            const unsigned int userOutputColumn=DN.userOutputColumn();
            model::cout<<greenColor<<"Initializing SampleUserStressController at runID="<<runID<<defaultColor<<std::endl;
            
                        model::EigenDataReader EDR;
            EDR.readScalarInFile(inputFileName,"use_externalStress",DN.use_externalStress);
            if(DN.use_externalStress)
            {
                EDR.readScalarInFile(inputFileName,"relaxSteps",relaxSteps);

                EDR.readMatrixInFile(inputFileName,"externalStress",_externalStress);
                assert((_externalStress-_externalStress.transpose()).norm()<DBL_EPSILON && "externalStress is not symmetric.");
                
                EDR.readMatrixInFile(inputFileName,"P0",P0);

                EDR.readScalarInFile(inputFileName,"tau_xz",tau_xz);

                
                DN._userOutputColumn+=0;

                meshSize=DN.mesh.xMax()-DN.mesh.xMin();
                          std::cout<<"meshSize="<<meshSize.transpose()<<std::endl;
                          
                          nu=Material<Isotropic>::nu;
                          lambda=2.0*nu/(1.0-2.0*nu);

                          std::cout<<"nu="<<nu<<std::endl;
                          std::cout<<"lambda="<<lambda<<std::endl;
                
                dislocationImages_x=DN.dislocationImages_x;
                dislocationImages_y=DN.dislocationImages_y;
                dislocationImages_z=DN.dislocationImages_z;

            }
            
        }
        

        /**************************************************************************/
        MatrixDim externalStress(const VectorDim& x) const
        {
            
            double tau_13=0.0;
            
            // Consider stress of images
            for(int i=-dislocationImages_x;i<=dislocationImages_x;++i)
            {
                for(int j=-dislocationImages_y;j<=dislocationImages_y;++j)
                {
                    for(int k=-dislocationImages_z;k<=dislocationImages_z;++k)
                    {
                        const Eigen::Matrix<int,3,1> cellID((Eigen::Matrix<int,3,1>()<<i,j,k).finished());
                        const VectorDim P=P0+(meshSize.array()*cellID.cast<double>().array()).matrix(); // image of P0 in cell (i,j,k)
                        tau_13+=tau_xz/(x-P).norm(); // dummy stress using x and P
                    }
                }
            }
            
            
            MatrixDim temp(MatrixDim::Zero());
            temp(0,2)=tau_13;
            temp(2,0)=temp(0,2);
            
            return _externalStress+temp;
        }
        /*************************************************************************/
        template <typename DislocationNetworkType>
        void update(const DislocationNetworkType& )
        {
            
        }
        
        /**************************************************************************/
        void output(const long int& ,
                    UniqueOutputFile<'F'>& ,
                    std::ofstream& ,
                    int& ) const
        {
            

        }
        
    };
}
#endif
