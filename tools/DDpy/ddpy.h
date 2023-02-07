
#ifndef MODELIB2PY_H
#define MODELIB2PY_H

#define _MODEL_NON_SINGULAR_DD_ 1 // 0 classical theory, 1 Cai's regularization method, 2 Lazar's regularization method

#include <string>
#include <tuple>
#include <list>
#include <map>
#include <stdlib.h> // EXIT_SUCCESS, EXIT_FAILURE

#include <pybind11/pybind11.h>
//#include <pybind11/smart_holder.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include <Polycrystal.h>
//#include <SimplicialMesh.h>
//#include <DDconfigIO.h>
//#include <Grain.h>
#include <DefectiveCrystal.h>
#include <MicrostructureGenerator.h>
#include <MicrostructureGeneratorBase.h>
#include <DislocationDynamicsBase.h>
//#include <DislocationSegment.h>
//#include <DislocationDynamicsModule.h>
#include <UniformExternalLoadController.h>
#include <IDreader.h>
#include <string>
#include <Eigen/Core>

namespace ddpy
{
class DDInterface
{
   private:
      typedef model::DefectiveCrystal<3,0> DefectiveCrystalType;
      std::string folderPath;
      model::DislocationDynamicsBase<3> ddBase;
      DefectiveCrystalType DC;
   public:
      typedef Eigen::Matrix<double,3,1> VectorDim;
      typedef Eigen::Matrix<double,3,3> MatrixDim;
      DDInterface( const std::string& folderName):
          folderPath( folderName)
          ,ddBase( model::DislocationDynamicsBase<3>(folderPath))
          ,DC( ddBase)
      {};

      std::list<std::tuple< size_t, size_t, double>>
         getResolvedShearStresses();
      std::list<std::tuple< size_t, size_t, double>>
         getResolvedShearStrains();
      std::list<std::tuple< size_t, size_t, double>>
         getPlasticStrains();
      ////std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
      //std::map< std::pair<size_t,size_t>, VectorDim>
      //   getSlipSystemNormals() const;
      //std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
      //std::map< std::pair<size_t,size_t>, Eigen::Matrix<double,3,1> >
      //py::dict 
      //   getSlipSystemBurgersVectors() const;

      std::string getFolderPath(){ return folderPath;}
      void setEndingStep( const long int& endingStep)
      {
         //if ( DC == NULL)
         //{
         //   std::cout << "Error: cannot assign ending step until "
         //      << "DefectiveCrystal is instantiated" << std::endl;
         //   return;
         //}
         if ( endingStep <= DC.simulationParameters.Nsteps)
         {
            std::cout << "Warning: assigning endingStep=" << endingStep
               << " <= currentStep:" << DC.simulationParameters.Nsteps
               << std::endl;
         }
         DC.simulationParameters.Nsteps = endingStep;
         return;
      }

      size_t getCurrentStep()
      {
         //if ( DC == NULL)
         //{
         //   std::cout << "Error: cannot getCurrentStep until "
         //      << "DefectiveCrystal is instantiated" << std::endl;
         //   return 0;
         //}
         return DC.simulationParameters.runID;
      }

      void setCurrentStep( const long int& step)
      {
         //if ( DC == NULL)
         //{
         //   std::cout << "Error: cannot setCurrentStep until "
         //      << "DefectiveCrystal is instantiated" << std::endl;
         //   return;
         //}

         DC.simulationParameters.runID = step;
         DC.externalLoadController = DC.getExternalLoadController(
                  DC.simulationParameters, DC, step
                  );
         DC.simulationParameters.manageRestart();
         return;
      }

      void runGlideSteps( size_t Nsteps)
      {
         //if ( DC == NULL)
         //{
         //   DC = new DefectiveCrystalType( ddBase);
         //}
         setEndingStep( getCurrentStep() + Nsteps);
         DC.runGlideSteps();
         return;
      }

      //void regenerate()
      //{
      //   //if ( DC != NULL)
      //   //{
      //   //   std::cout << "new configuration evl/evl_0.txt generated, "
      //   //      << "but I won't guarantee that the new DefectiveCrystal"
      //   //      << " will load the 0'th timestep if others are present"
      //   //      << " in " << folderPath << "/evl/ and "
      //   //      << folderPath << "/F/" << std::endl;
      //   //   //std::cout << "Instantiating DefectiveCrystal with new "
      //   //   //   << "configuration" << std::endl;
      //   //   delete DC;
      //   //}
      //   //DefectiveCrystalType::~DefectiveCrystal( DC);
      //   DC.~DefectiveCrystal();
      //   model::StaticID<model::
      //      >::set_count(0);
      //   model::MicrostructureGenerator mg( ddBase);
      //   DefectiveCrystalType DC( ddBase);
      //   // TODO: set the time step of DC to 0
      //   return;
      //}
};

//   //std::map< size_t, VectorDim> printSlipSystemNormals()
//   //std::map< std::pair<size_t,size_t>, model::ReciprocalLatticeVector<3>>
//   // printSlipSystemNormals()
//   //std::map< std::pair<size_t,size_t>, VectorDim> printSlipSystemNormals()
//   //std::map< std::pair<size_t,size_t>, std::string> printSlipSystemNormals()



//   void generateMicrostructure();
//}; // class DefectiveCrystalInterface

} // namespace ddpy

#endif
