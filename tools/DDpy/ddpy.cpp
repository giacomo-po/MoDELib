
#ifndef MODELIB2PY_CPP
#define MODELIB2PY_CPP

#include <ddpy.h>

//ddpy::DefectiveCrystalInterface::DefectiveCrystalInterface( const std::string& folderName) : model::DefectiveCrystal<3,0>( folderName) {}

PYBIND11_MODULE( ddpy, m) {
   namespace py = pybind11;
   m.doc() = "TODO: revise m.doc() in src/modelib2py.cpp";
   //py::class_<model::DefectiveCrystal<3,0>>( m, "DefectiveCrystal")
   //   .def( py::init([](){
   //            return std::unique_ptr<model::DefectiveCrystal<3,0>>( new model::DefectiveCrystal<3,0>(std::string("./")));
   //            }))
   py::class_<ddpy::DefectiveCrystalInterface>( m, "DefectiveCrystalInterface")
      .def( py::init([](const std::string& folderName){ // lambda function that returns an instantiation
               return std::unique_ptr< ddpy::DefectiveCrystalInterface>(
                     new ddpy::DefectiveCrystalInterface( folderName)
                     );
               }), py::arg("folderPath"))
      .def("getCurrentStep",
            &ddpy::DefectiveCrystalInterface::getCurrentStep)
      .def("setCurrentStep",
            &ddpy::DefectiveCrystalInterface::setCurrentStep,
            py::arg("currentStep")
          )
      .def("getEndingStep",
            &ddpy::DefectiveCrystalInterface::getEndingStep
          )
      .def("setEndingStep",
            &ddpy::DefectiveCrystalInterface::setEndingStep,
            py::arg("endingStep")
          )
      .def("runGlideSteps",
            &ddpy::DefectiveCrystalInterface::runGlideSteps
          )
      .def("printSlipSystemNormals",
            &ddpy::DefectiveCrystalInterface::printSlipSystemNormals
          )
      .def("printSlipSystemBurgersVectors",
            &ddpy::DefectiveCrystalInterface::printSlipSystemBurgersVectors
          )
      //.def("printResolvedShearStress",
      //      &ddpy::DefectiveCrystalInterface::printResolvedShearStress
      //    )
      ;
   py::class_<ddpy::MicrostructureGeneratorInterface>( m, "MicrostructureGeneratorInterface")
      .def( py::init([](const std::string& folderName){ // lambda function that returns an instantiation
               return std::unique_ptr< ddpy::MicrostructureGeneratorInterface>(
                     new ddpy::MicrostructureGeneratorInterface( folderName)
                     );
               }), py::arg("folderPath"))
      // TODO: expose any other class members to pybind here
      //.def("printDisplacements", &ddpy::DefectiveCrystalInterface::printDisplacements)
      ;
} //  PYBIND11_MODULE

//std::tuple<double>
//ddpy::DefectiveCrystalInterface::printResolvedShearStress()
//   {
//      size_t grainCount = 0;
//      size_t slipSystemCount = 0;
//      //Eigen::Matrix< double, 3, 3> stress; 
//      MatrixDim stress;
//      MatrixDim strain;
//      std::list<double> rss;
//      if ( this->externalLoadController != NULL) 
//      {
//         stress = this->externalLoadController->ExternalStress;
//      }
//      else
//      {
//         std::cout << "error: printResolvedShearStress(), "
//           << " this->externalLoadController == NULL" << std::endl;
//         return std::make_tuple(rss);
//      }
//      for ( const auto& grain : this->poly.grains)
//      {
//         ++grainCount;
//         rss.clear();
//         std::cout << "grain " << grainCount << std::endl;
//         for ( const auto& ss : grain.second.singleCrystal->slipSystems())
//         { // loop over slip system
//            ++slipSystemCount;
//            std::cout << "slip system " << slipSystemCount 
//               << " planeNormal: " << std::endl
//               << ss->unitNormal << std::endl
//               << " burgers vector: " << std::endl
//               << ss->s.cartesian() << std::endl;
//
//            std::cout << "stress: " << std::endl << stress << std::endl;
//            std::cout << "stress * planeNormal: " << std::endl
//               << stress * (ss->unitNormal) << std::endl;
//            std::cout << "rss: "
//               <<  (stress * (ss->unitNormal)).dot(ss->s.cartesian())
//               << std::endl;
//            rss.push_back(
//                  (stress * (ss->unitNormal)).dot(ss->s.cartesian()));
//         } // loop over slip system
//         return std::make_tuple(rss);
//      }
//      return;
//   }

//   //////model::Grain<3> myGrain = DC.poly.grain( 1);
//   //////return myGrain.slipSystems().size();
//   ////Eigen::Matrix<double,3,1> nn;
//   ////std::cout << "number of grains: " 
//   ////   << DC.poly.grains.size() << std::endl;
//   ////std::cout << "number of slip systems: " 
//   ////   << DC.poly.grain( 1).singleCrystal->slipSystems().size() << std::endl;
//   ////Eigen::Matrix<double,3,1> bb;
//   //////Eigen::Matrix<double,3,1> bbTemp;
//   ////Eigen::Matrix< double, 3, 3> stress; 
//   //////std::vector<double> rss;
//   ////double rss=0.0;
//   ////Eigen::Matrix< double, 3, 1> 
//   ////   position( Eigen::Matrix< double, 3, 1>::Zero());
//
//   ////for ( size_t ssIdx=0; ssIdx < DC.poly.grain( 1).singleCrystal->slipSystems().size();
//   ////      ++ssIdx)
//   ////{
//   ////   bb = DC.poly.grain( 1).singleCrystal->slipSystems()[ssIdx]->s.cartesian();
//   ////   nn = DC.poly.grain( 1).singleCrystal->slipSystems()[ssIdx]->unitNormal;
//   ////   std::cout << "slip system: " << ssIdx << std::endl; // debug
//   ////   std::cout << "  Burgers b : " << bb[0] << ", " << bb[1] << ", " << bb[2] << std::endl; // debug
//   ////   std::cout << "  slip plane normal n : " << nn[0] << ", " << nn[1] << ", " << nn[2] << std::endl; // debug
//
//   ////   if ( DC.DN->externalLoadController != NULL)
//   ////      stress = DC.DN->externalLoadController->stress( position);
//   ////   else
//   ////   {
//   ////      std::cout << "DC.DN->externalLoadConstroller is null" << std::endl; // debug
//   ////      continue;
//   ////   }
//
//   ////   // TODO: calculate (sigma * nn) dot b/|b|
//   ////   rss = ( stress * nn).dot( bb);  
//   ////   std::cout << "  resolved shear stress: " << rss << std::endl;
//   ////}
//
//
//
//   //// read step beginStep
//   ////DC.simulationParameters.runID = beginStep;
//   //////  time step delta t is accessible from DN.simulationParameters.dt
//   ////std::cout << "dt : " << DC.simulationParameters.dt << std::endl;
//   ////std::cout << "runID : " << DC.simulationParameters.runID << std::endl;
//   ////std::cout << "stress tensor: " << std::endl
//   ////   << " " << stress(0,0) 
//   ////   << " " << stress(0,1) 
//   ////   << " " << stress(0,2) << std::endl
//   ////   << " " << stress(1,0) 
//   ////   << " " << stress(1,1) 
//   ////   << " " << stress(1,2) << std::endl
//   ////   << " " << stress(2,0) 
//   ////   << " " << stress(2,1) 
//   ////   << " " << stress(2,2) << std::endl;
//
//   ////std::cout << "n[1] : " << nn[1] << std::endl;
//   ////std::cout << "setting Nsteps to " << runSteps << std::endl;
//   ////DC.simulationParameters.Nsteps = runSteps;
//   ////std::cout << "running " << DC.simulationParameters.Nsteps
//   ////   << " steps" << std::endl;
//   //DC.runGlideSteps();
//   //std::cout << "finished running " << runSteps << " steps" << std::endl;
//   //return 0;
//   ////return myGrain.slipSystems()[0].s.cartesian().normalized;
//}
//
////double get_something_from_modelib( const size_t& timeStep)
////{
////   std::string folderName = "evl";
////   std::string suffix;
////   Eigen::Matrix<double, 3, 3> rotationMatrix;
////   rotationMatrix( 0, 0) = 1.0;
////   rotationMatrix( 0, 1) = 0.0;
////   rotationMatrix( 0, 2) = 0.0;
////   rotationMatrix( 1, 0) = 0.0;
////   rotationMatrix( 1, 1) = 1.0;
////   rotationMatrix( 1, 2) = 0.0;
////   rotationMatrix( 2, 0) = 0.0;
////   rotationMatrix( 2, 0) = 0.0;
////   rotationMatrix( 2, 1) = 0.0;
////   rotationMatrix( 2, 2) = 1.0;
////
////   Eigen::Matrix<double, 3, 1> displacement;
////
////   // model::SingleCrystal
////   // TODO: something that instantiates a PolyCrystal and DDconfigIO
////   std::string polycrystalFilePath = "inputFiles/polycrystal.txt";
////   std::string meshFileName = "MeshLibrary/small_block_structured1_fine_scaled_2order.msh";
////   // TODO: read periodicFaceIDs from polycrystal.txt
////   std::set<int> periodicFaceIDs( {0, 1, 2, 3, 4, 5});
////   
////   //model::SimplicialMesh<3> mesh;
////   model::SimplicialMesh<3> mesh( polycrystalFilePath, rotationMatrix, 
////         displacement, periodicFaceIDs);
////   mesh.readMesh( meshFileName, rotationMatrix, displacement, periodicFaceIDs);
////   model::Polycrystal<3> polyCrystal( polycrystalFilePath, mesh);
////   model::DDconfigIO<3> modelibReader( folderName, suffix);
////   std::cout << "reading time step " << timeStep << std::endl; // debug
////   modelibReader.read( timeStep); // read in simulated step 1
////   //model::Grain<3> myGrain = polyCrystal.grains().at( timeStep);
////   //model::Grain<3> myGrain = polyCrystal.grain( 1);
////   //model::Grain<3> myGrain = polyCrystal.grain( timeStep);
////   //size_t grainSize = polyCrystal.grain( timeStep)
////   //size_t ssSize = polyCrystal.grain( timeStep).slipSystems().size();
////   size_t grainsSize = polyCrystal.grains.size();
////   //myGrain.s
////   
////   return grainsSize; // arbitrary return value because I'm just messing around
////   //return myGrain.slipSlipSystems().size();
////}
//

   //py::object world = py::cast("World"); // explicit conversion to py object
   //m.attr("what") = world;

   // notes:
   //  instantiate a class Polycrystal, providing a container of 
   //  SlipSystems
   //  SlipSystem.unitNormal provides \hat{n} normal to the slip plane
   //  SlipSystem.s is a lattice vector in units of lattice spacing
   //  SlipSystem.s.cartesian() returns a unit vector parallel to the
   //   slip system's Burger's vector
   //  resolved shear stress would then be ((total stress tensor)*(slip system unit normal \hat{n}))*(SlipSystem.s.cartesian())
//}
#endif
