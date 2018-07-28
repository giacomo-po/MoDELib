/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *                       Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ExtraStraightSegmentsController_H_
#define model_ExtraStraightSegmentsController_H_

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
    class ExtraStraightSegmentsController;
    
    template <int dim>
    using ExternalLoadController=ExtraStraightSegmentsController<dim>;

    
    template <int dim>
    class ExtraStraightSegmentsController
    {
        
//        static constexpr int voigtSize=dim*(dim+1)/2;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
        const std::string inputFileName;
        MatrixDim _externalStress;
//        VectorDim P0;
        int relaxSteps;

//        //External stress control parameter
//        MatrixDim _externalStress0;
//        MatrixDim _externalStressRate;
//        
//        MatrixDim ExternalStrain;
//        //External strain control parameter
//        MatrixDim ExternalStrain0;
//        MatrixDim ExternalStrainRate;
//        MatrixDim plasticStrain;
//        //finite machine stiffness effect
//        Eigen::Matrix<double,1,voigtSize> MachineStiffnessRatio;    //0 is stress control; infinity is pure strain control.
//        Eigen::Matrix<size_t,voigtSize,2>     voigtorder;
//        Eigen::Matrix<double,voigtSize,voigtSize> stressmultimachinestiffness;
//        Eigen::Matrix<double,voigtSize,voigtSize> strainmultimachinestiffness;
//        double last_update_time;
        double nu;  //unit is mu, lambda=2*v/(1-2*v)
        double lambda;  //unit is mu, lambda=2*v/(1-2*v)
//        // Number of DD steps before applying strain.
//        double nu_use;
//        double sample_volume;
        
        
    public:
        /**************************************************************************/
        ExtraStraightSegmentsController():
        /* init list */ inputFileName("./externalLoadControl/ExtraStraightSegmentsController.txt")
        /* init list */,_externalStress(MatrixDim::Zero())
//        /* init list */,P0(VectorDim::Zero())
//        /* init list */ _externalStress0(MatrixDim::Zero()),
//        /* init list */ _externalStressRate(MatrixDim::Zero()),
//        /* init list */ ExternalStrain(MatrixDim::Zero()),
//        /* init list */ ExternalStrain0(MatrixDim::Zero()),
//        /* init list */ ExternalStrainRate(MatrixDim::Zero()),
//        /* init list */ plasticStrain(MatrixDim::Zero()),
//        /* init list */ MachineStiffnessRatio(Eigen::Matrix<double,1,voigtSize>::Zero()),
//        /* init list */ voigtorder(Eigen::Matrix<size_t,voigtSize,2>::Zero()),
//        /* init list */ stressmultimachinestiffness(Eigen::Matrix<double,voigtSize,voigtSize>::Zero()),
//        /* init list */ strainmultimachinestiffness(Eigen::Matrix<double,voigtSize,voigtSize>::Zero()),
//        /* init list */ last_update_time(0.0),
//        /* init list */ lambda(1.0),
//        /* init list */ nu_use(0.12),
//        /* init list */ sample_volume(1.0),
        /* init list */,relaxSteps(0)
        /* init list */,nu(0.33)
        /* init list */,lambda(2.0*nu/(1.0-2.0*nu))

        {
            
        }
        
        /**************************************************************************/
        template <typename DislocationNetworkType>
        void init(const DislocationNetworkType& DN)
        {
            const long int runID=DN.runningID();
            const unsigned int userOutputColumn=DN.userOutputColumn();
            model::cout<<greenColor<<"Initializing External Stress Field Controller at runID="<<runID<<defaultColor<<std::endl;
            
            
            
            model::EigenDataReader EDR;
            EDR.readMatrixInFile(inputFileName,"externalStress",_externalStress);
            assert((_externalStress-_externalStress.transpose()).norm()<DBL_EPSILON && "externalStress is not symmetric.");
            
//            EDR.readMatrixInFile(inputFileName,"P0",P0);
            EDR.readScalarInFile(inputFileName"relaxSteps",relaxSteps);
            
            
            EDR.readScalarInFile(inputFileName,"use_extraStraightSegments",DN.use_extraStraightSegments);
            if (DN.use_extraStraightSegments)
            {
                DN.ssdeq.clear();
                typedef IDreader<'B',1,10,double> IDreaderType;
                IDreaderType vReader;
                if (vReader.isGood(0,true))
                {
                    vReader.read(0,true);
                    for (const auto& vIter : vReader)
                    {
                        Eigen::Map<const Eigen::Matrix<double,1,9>> row(vIter.second.data());
                        VectorDimD P0(row.template segment<dim>(0));// P0 position
                        VectorDimD P1(row.template segment<dim>(dim)); // P1 position
                        VectorDimD B(row.template segment<dim>(dim*2));  // Burgers vector
                        DN.ssdeq.emplace_back(StressStraight<dim>(P0,P1,B));
                    }
                }
                else
                {
                    model::cout<<"could not read runID from B/B_0.txt"<<std::endl;
                }
            }
            
            
            EDR.readScalarInFile(inputFileName,"use_externalStress",DN.use_externalStress);
            if (DN.use_externalStress)
            {
                DN._userOutputColumn+=0; 
            }
            
            
            nu=Material<Isotropic>::nu;
//            std::cout<<" nu="<<nu<<std::endl;
            lambda=2.0*nu/(1.0-2.0*nu);

        }
        

        /**************************************************************************/
        const MatrixDim& externalStress(const VectorDim& x) const
        {
            return _externalStress;
        }
        /*************************************************************************/
        template <typename DislocationNetworkType>
        void update(const DislocationNetworkType& )
        {
            
        }
        
        /**************************************************************************/
        std::string output() const
        {
            return std::string();
        }
        
    };
}
#endif
