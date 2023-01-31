
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
//#include <DislocationSegment.h>
//#include <DislocationDynamicsModule.h>
#include <UniformExternalLoadController.h>
#include <IDreader.h>
#include <string>
#include <Eigen/Core>

namespace ddpy
{
class MicrostructureGeneratorInterface : public model::MicrostructureGenerator
{
   private:
      //typedef model::DefectiveCrystal<3,0> DefectiveCrystalType;
      //typedef Eigen::Matrix<double,3,1> VectorDim;
      //typedef Eigen::Matrix<double,3,3> MatrixDim;
   public:
      MicrostructureGeneratorInterface(
            const std::string& folderName
            ) : MicrostructureGenerator( folderName)
      {}
};

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

   //std::map< size_t, VectorDim> printSlipSystemNormals()
   //std::map< std::pair<size_t,size_t>, model::ReciprocalLatticeVector<3>>
   // printSlipSystemNormals()
   //std::map< std::pair<size_t,size_t>, VectorDim> printSlipSystemNormals()
   std::list<std::tuple< size_t, size_t, VectorDim>> printSlipSystemNormals()
   //std::map< std::pair<size_t,size_t>, std::string> printSlipSystemNormals()
   {
      size_t grainCount = 0;
      size_t slipSystemCount = 0;
      //std::map< std::pair<size_t,size_t>, model::ReciprocalLatticeVector<3>> normals;
      //std::map< std::pair<size_t,size_t>, VectorDim> normals;
      std::list<std::tuple< size_t, size_t, VectorDim>> normals;
      //std::map< std::pair<size_t,size_t>, std::string> normals;
      // ((grain number, slip system number), plane normal)
      for ( const auto& grain : this->poly.grains)
      {
         std::cout << "grain " << grainCount << std::endl;
         for ( const auto& ss : grain.second.singleCrystal->slipSystems())
         { // loop over slip system
            std::cout << "slip system " << slipSystemCount
               << " plane normal:" << std::endl
               << " (" << std::endl << ss->unitNormal
               << ")" << std::endl;
            normals.push_back(
                  std::tuple( grainCount, slipSystemCount, ss->unitNormal)
                  );
            ++slipSystemCount;
         }
         ++grainCount;
      }
      for ( const auto& nn : normals)
      {
         std::cout << "grain " << std::get<0>( nn)
            << " slip system " << std::get<1>( nn)
            << std::endl << std::get<2>( nn) << std::endl;
      }
      return normals;
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
               << " (" << std::endl << ss->unitSlip
               << ")" << std::endl;
         }
      }
   }

   std::list<std::tuple< size_t, size_t, double>>
      getResolvedShearStresses();
   std::list<std::tuple< size_t, size_t, double>>
      getResolvedShearStrains();
   std::list<std::tuple< size_t, size_t, double>>
      getPlasticStrains();

   void generateMicrostructure();
}; // class DefectiveCrystalInterface

} // namespace ddpy

#endif
