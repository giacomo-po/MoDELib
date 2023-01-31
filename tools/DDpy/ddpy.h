
#ifndef MODELIB2PY_H
#define MODELIB2PY_H

#define _MODEL_NON_SINGULAR_DD_ 1 // 0 classical theory, 1 Cai's regularization method, 2 Lazar's regularization method

#include <string>
#include <tuple>
#include <list>
#include <stdlib.h> // EXIT_SUCCESS, EXIT_FAILURE

#include <pybind11/pybind11.h>
//#include <pybind11/smart_holder.h>

#include <Polycrystal.h>
//#include <SimplicialMesh.h>
//#include <DDconfigIO.h>
//#include <Grain.h>
#include <DefectiveCrystal.h>
//#include <DislocationSegment.h>
//#include <DislocationDynamicsModule.h>
#include <UniformExternalLoadController.h>
#include <IDreader.h>
#include <VTKsegments.h>
#include <string>
#include <Eigen/Core>

namespace model2py
{
//template< std::string folderName>
class DefectiveCrystalInterface : public model::DefectiveCrystal<3,0>
{
   private:
      typedef model::DefectiveCrystal<3,0> DefectiveCrystalType;
      typedef Eigen::Matrix<double,3,1> VectorDim;
      typedef Eigen::Matrix<double,3,3> MatrixDim;
   public:
   DefectiveCrystalInterface(
         const std::string& folderName
         ) : DefectiveCrystalType( folderName)
   {}

   size_t getCurrentStep()
   {
      return this->simulationParameters.runID;
   }
   void setCurrentStep( const long int& currentStep)
   {
      this->simulationParameters.runID = currentStep;
      this->externalLoadController = getExternalLoadController(
               simulationParameters, *this, currentStep
               );
      this->simulationParameters.manageRestart(); // TODO: is this necessary or better than ^^^?
   }

   void setEndingStep( const long int& endingStep)
   {
      if ( endingStep <= this->simulationParameters.Nsteps)
      {
         std::cout << "Warning: assigning endingStep=" << endingStep
            << " <= currentStep:" << this->simulationParameters.Nsteps
            << std::endl;
      }
      this->simulationParameters.Nsteps = endingStep;
      return;
   }

   long int getEndingStep()
   {
      return this->simulationParameters.Nsteps;
   }

   void printSlipSystemNormals()
   {
      size_t grainCount = 0;
      size_t slipSystemCount = 0;
      for ( const auto& grain : this->poly.grains)
      {
         ++grainCount;
         std::cout << "grain " << grainCount << std::endl;
         for ( const auto& nn : grain.second.singleCrystal->planeNormals())
         {
            ++slipSystemCount;
            std::cout << "slip system " << slipSystemCount 
               << " planeNormal: " << std::endl
               << " (" << nn->row(0) 
               << ", " << nn->row(1) 
               << ", " << nn->row(2) 
               << ")" << std::endl;
         }
      }
      return;
   }

   void printSlipSystemBurgersVectors()
   {
      size_t grainCount = 0;
      size_t slipSystemCount = 0;
      for ( const auto& grain : this->poly.grains)
      {
         ++grainCount;
         std::cout << "grain " << grainCount << std::endl;
         for ( const auto& ss : grain.second.singleCrystal->slipSystems())
         { // loop over slip system
            ++slipSystemCount;
            std::cout << "slip system " << slipSystemCount 
               << " burgers vector:" << std::endl
               << " (" << ss->s.cartesian().row(0) 
               << ", " << ss->s.cartesian().row(1) 
               << ", " << ss->s.cartesian().row(2) 
               << ")" << std::endl;
         }
      }
   }

   //std::tuple<double> printResolvedShearStress();

   //void model2py::setCurrentStep( size_t currentStep)
}; // class DefectiveCrystalInterface 

//class DDsegments
//{
//   public:
//      DDsegments( const std::string& folderName)
//
//      {
//         VTKsegments vs( folderName);
//         vs.readVTK( vtkFile);
//      }
//}; // class DDsegments

} // namespace model2py

#endif
