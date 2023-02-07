
#ifndef MODELIB2PY_CPP
#define MODELIB2PY_CPP

#include <ddpy.h>

std::list<std::tuple< size_t, size_t, double>>
ddpy::DDInterface::getResolvedShearStresses()
{
   size_t grainCount = 0;
   size_t slipSystemCount = 0;
   //Eigen::Matrix< double, 3, 3> stress;
   MatrixDim stress;
   std::list<std::tuple<size_t,size_t,double>> rss;
   //if ( DC == NULL) 
   //{
   //   std::cout << "error: getResolvedShearStresses(), "
   //     << " DefectiveCrystal not yet initialized" << std::endl;
   //   return rss;
   //}
   if ( DC.externalLoadController != NULL)
   {
      stress = DC.externalLoadController->stress( VectorDim::Zero());
   }
   else
   {
      std::cout << "error: getResolvedShearStresses(), "
        << " DC.externalLoadController == NULL" << std::endl;
      return rss;
   }

   for ( const auto& grain : DC.poly.grains)
   {
      rss.clear();
      //std::cout << "grain " << grainCount << std::endl;
      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
      { // loop over slip system
         //std::cout << "rss " << slipSystemCount << ": "
         //   <<  (stress * (ss->unitNormal)).dot(ss->unitSlip)
         //   << std::endl;
         rss.push_back(
                  std::tuple(
                     grainCount,
                     slipSystemCount,
                     (stress * (ss->unitNormal)).dot(ss->unitSlip)
                     )
               );

         ++slipSystemCount;
      } // loop over slip system
      ++grainCount;
   }

   //for ( const auto& itr : rss)
   //{
   //   std::cout << "(grain,slipSystem,stress): " << "("
   //      << std::get<0>(itr) << ","
   //      << std::get<1>(itr) << ","
   //      << std::get<2>(itr) << ")" << std::endl;
   //}

   return rss;
}

std::list<std::tuple< size_t, size_t, double>>
ddpy::DDInterface::getResolvedShearStrains()
{
   size_t grainCount = 0;
   size_t slipSystemCount = 0;
   MatrixDim strain;
   std::list<std::tuple<size_t,size_t,double>> rss;
   //if ( DC == NULL) 
   //{
   //   std::cout << "error: getResolvedShearStrains(), "
   //     << " DefectiveCrystal not yet initialized" << std::endl;
   //   return rss;
   //}
   if ( DC.externalLoadController != NULL)
   {
      strain = DC.externalLoadController->strain( VectorDim::Zero());
   }
   else
   {
      std::cout << "error: getResolvedShearStrains(), "
        << " DC.externalLoadController == NULL" << std::endl;
      return rss;
   }
   for ( const auto& grain : DC.poly.grains)
   {
      rss.clear();
      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
      { // loop over slip system
         rss.push_back(
                  std::tuple(
                     grainCount,
                     slipSystemCount,
                     (strain * (ss->unitNormal)).dot(ss->unitSlip)
                     )
               );
         ++slipSystemCount;
      } // loop over slip system
      ++grainCount;
   }
   //for ( const auto& itr : rss)
   //{
   //   std::cout << "(grain,slipSystem,strain): " << "("
   //      << std::get<0>(itr) << ","
   //      << std::get<1>(itr) << ","
   //      << std::get<2>(itr) << ")" << std::endl;
   //}
   return rss;
}

std::list<std::tuple< size_t, size_t, double>>
ddpy::DDInterface::getPlasticStrains()
{
   std::list<std::tuple< size_t, size_t, double>> plasticStrains;
   //if ( DC == NULL) 
   //{
   //   std::cout << "error: getPlasticStrains(), "
   //     << " DefectiveCrystal not yet initialized" << std::endl;
   //   return plasticStrains;
   //}
   std::map<std::pair<int,int>,double> sspd(
         DC.DN->slipSystemPlasticDistortion());
   // [[grain ID, slip system ID], strain value in that slip system]

   //std::cout << "sspd.size: " << sspd.size() << std::endl;

   // copy into list of tuples to be returned
   for (const auto& itr : sspd)
   {
       plasticStrains.push_back(
             std::tuple( itr.first.first, itr.first.second, itr.second)
             );
   }

   //for ( const auto& itr : plasticStrains)
   //{
   //   std::cout << "(grain,slipSystem,strain): " << "("
   //      << std::get<0>(itr) << ","
   //      << std::get<1>(itr) << ","
   //      << std::get<2>(itr) << ")" << std::endl;
   //}

   return plasticStrains;
}

//std::map<std::pair< size_t, size_t>, ddpy::DDInterface::VectorDim>
//   ddpy::DDInterface::getSlipSystemNormals() const
//{
//   size_t grainCount = 0;
//   size_t slipSystemCount = 0;
//   //std::map< std::pair<size_t,size_t>, model::ReciprocalLatticeVector<3>> normals;
//   std::map< std::pair<size_t,size_t>, VectorDim> normals;
//   //std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
//      normals;
//   if ( DC == NULL) 
//   {
//      std::cout << "error: getSlipSystemNormals(), "
//        << " DefectiveCrystal not yet initialized" << std::endl;
//      return normals;
//   }
//   //std::map< std::pair<size_t,size_t>, std::string> normals;
//   // ((grain number, slip system number), plane normal)
//   for ( const auto& grain : DC->poly.grains)
//   {
//      std::cout << "grain " << grainCount << std::endl;
//      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
//      { // loop over slip system
//         //std::cout << "slip system " << slipSystemCount
//         //   << " plane normal:" << std::endl
//         //   << " (" << std::endl << ss->unitNormal
//         //   << ")" << std::endl;
//         normals.push_back(
//               std::tuple( grainCount, slipSystemCount, ss->unitNormal)
//               );
//         ++slipSystemCount;
//      }
//      ++grainCount;
//   }
//   for ( const auto& nn : normals)
//   {
//      std::cout << "grain " << std::get<0>( nn)
//         << " slip system " << std::get<1>( nn)
//         << std::endl << std::get<2>( nn) << std::endl;
//   }
//   return normals;
//}

//std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
//std::map< std::pair<size_t,size_t>, Eigen::Matrix<double,3,1> >
//   ddpy::DDInterface::getSlipSystemBurgersVectors() const
//{
//   size_t grainCount = 0;
//   size_t slipSystemCount = 0;
//   //std::list<std::tuple< size_t, size_t, Eigen::Matrix<double,3,1>>>
//   std::map< std::pair<size_t,size_t>, Eigen::Matrix<double,3,1>>
//      burgersVectors;
//   if ( DC == NULL) 
//   {
//      std::cout << "error: getSlipSystemBurgersVectors(), "
//        << " DefectiveCrystal not yet initialized" << std::endl;
//      return burgersVectors;
//   }
//   for ( const auto& grain : DC->poly.grains)
//   {
//      ++grainCount;
//      std::cout << "grain " << grainCount << std::endl;
//      for ( const auto& ss : grain.second.singleCrystal->slipSystems())
//      { // loop over slip system
//         ++slipSystemCount;
//         burgersVectors.emplace(
//               std::make_pair(
//                  std::make_pair( grainCount, slipSystemCount),
//                  ss->unitSlip)
//               );
//         //std::cout << "slip system " << slipSystemCount
//         //   << " burgers vector:" << std::endl
//         //   << " (" << std::endl << ss->unitSlip
//         //   << ")" << std::endl;
//      }
//   }
//   return burgersVectors;
//}

PYBIND11_MODULE( ddpy, m) {
   namespace py = pybind11;
   m.doc() = "TODO: revise m.doc() in src/modelib2py.cpp";

   //py::class_<ddpy::DDInterface>( m, "DDInterface")
   //   .def( py::init([](const std::string& folderName){ // lambda function that returns an instantiation
   //            return std::unique_ptr< ddpy::DDInterface>(
   //                  new ddpy::DDInterface( folderName)
   //                  );
   //            }), py::arg("folderPath"))
   //   .def("",
   //         &ddpy::DDInterface::getPlasticStrains
   //       )
   //   ;
   py::class_<ddpy::DDInterface>( m, "DDInterface")
      .def(
            py::init(
               // using a lambda function that returns an instantiation
               [](const std::string& folderName)
               {
                  return std::unique_ptr< ddpy::DDInterface>(
                        new ddpy::DDInterface( folderName)
                        );
               }), py::arg("folderPath")
            )
      //.def("generate",
      //      &ddpy::DDInterface::generate
      //    )
      .def("getCurrentStep",
            &ddpy::DDInterface::getCurrentStep
          )
      .def("setCurrentStep",
            &ddpy::DDInterface::setCurrentStep
          )
      .def("runGlideSteps",
            &ddpy::DDInterface::runGlideSteps,
            py::arg("Nsteps").none(false)
          )
      .def("getResolvedShearStresses",
            &ddpy::DDInterface::getResolvedShearStresses
          )
      .def("getResolvedShearStrains",
            &ddpy::DDInterface::getResolvedShearStrains
          )
      .def("getPlasticStrains",
            &ddpy::DDInterface::getPlasticStrains
          )
      //.def("getSlipSystemNormals"
      //      &ddpy::DDInterface::getSlipSystemNormals
      //    )
      //.def("getSlipSystemBurgersVectors"
      //      &ddpy::DDInterface::getSlipSystemBurgersVectors
      //    )
      ;
   //py::class_<ddpy::DDInterface>( m, "DDInterface")
   //   .def(
   //         py::init(
   //            // using a lambda function that returns an instantiation
   //            [](const std::string& folderName)
   //            {
   //               model::DislocationDynamicsBase ddBase( folderName);
   //               return std::unique_ptr< model::DefectiveCrystal>(
   //                     new model::DefectiveCrystal( ddBase)
   //                     );
   //            }
   //            ), py::arg("folderPath")
   //         )
   //   .def("getCurrentStep",
   //         &ddpy::DDInterface::getCurrentStep
   //         )
   //   .def("setCurrentStep",
   //         &ddpy::DDInterface::setCurrentStep,
   //         py::arg("currentStep")
   //       )
   //   .def("runGlideSteps",
   //         &ddpy::DDInterface::runGlideSteps
   //       );
   //   // TODO: expose any other class members to pybind here
   //   //.def("printDisplacements", &ddpy::DDInterface::printDisplacements)
   //   ;
} //  PYBIND11_MODULE

#endif
