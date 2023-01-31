
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
      .def("getResolvedShearStresses",
            &ddpy::DefectiveCrystalInterface::getResolvedShearStresses
          )
      .def("getResolvedShearStrains",
            &ddpy::DefectiveCrystalInterface::getResolvedShearStrains
          )
      .def("getPlasticStrains",
            &ddpy::DefectiveCrystalInterface::getPlasticStrains
          )
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

std::list<std::tuple< size_t, size_t, double>>
ddpy::DefectiveCrystalInterface::getResolvedShearStresses()
{
   size_t grainCount = 0;
   size_t slipSystemCount = 0;
   //Eigen::Matrix< double, 3, 3> stress;
   MatrixDim stress;
   std::list<std::tuple<size_t,size_t,double>> rss;
   if ( this->externalLoadController != NULL)
   {
      stress = this->externalLoadController->stress( VectorDim::Zero());
   }
   else
   {
      std::cout << "error: getResolvedShearStresses(), "
        << " this->externalLoadController == NULL" << std::endl;
      return rss;
   }

   for ( const auto& grain : this->poly.grains)
   {
      rss.clear();
      std::cout << "grain " << grainCount << std::endl;
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

   for ( const auto& itr : rss)
   {
      std::cout << "(grain,slipSystem,stress): " << "("
         << std::get<0>(itr) << ","
         << std::get<1>(itr) << ","
         << std::get<2>(itr) << ")" << std::endl;
   }

   return rss;
}

std::list<std::tuple< size_t, size_t, double>>
ddpy::DefectiveCrystalInterface::getResolvedShearStrains()
{
   size_t grainCount = 0;
   size_t slipSystemCount = 0;
   MatrixDim strain;
   std::list<std::tuple<size_t,size_t,double>> rss;
   if ( this->externalLoadController != NULL)
   {
      strain = this->externalLoadController->strain( VectorDim::Zero());
   }
   else
   {
      std::cout << "error: getResolvedShearStrains(), "
        << " this->externalLoadController == NULL" << std::endl;
      return rss;
   }
   for ( const auto& grain : this->poly.grains)
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
   for ( const auto& itr : rss)
   {
      std::cout << "(grain,slipSystem,strain): " << "("
         << std::get<0>(itr) << ","
         << std::get<1>(itr) << ","
         << std::get<2>(itr) << ")" << std::endl;
   }
   return rss;
}

std::list<std::tuple< size_t, size_t, double>>
ddpy::DefectiveCrystalInterface::getPlasticStrains()
{
   std::list<std::tuple< size_t, size_t, double>> plasticStrains;
   std::map<std::pair<int,int>,double> sspd(
         this->DN->slipSystemPlasticDistortion());
   // [[grain ID, slip system ID], strain value in that slip system]

   //std::cout << "sspd.size: " << sspd.size() << std::endl;

   // copy into list of tuples to be returned
   for (const auto& itr : sspd)
   {
       plasticStrains.push_back(
             std::tuple( itr.first.first, itr.first.second, itr.second)
             );
   }

   for ( const auto& itr : plasticStrains)
   {
      std::cout << "(grain,slipSystem,strain): " << "("
         << std::get<0>(itr) << ","
         << std::get<1>(itr) << ","
         << std::get<2>(itr) << ")" << std::endl;
   }

   return plasticStrains;
}

#endif
