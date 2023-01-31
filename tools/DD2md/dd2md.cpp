
#include <dd2md.h>

PYBIND11_MODULE( dd2md, m) {
   namespace py = pybind11;
   m.doc() = "TODO: revise m.doc() in dd2md.cpp";
   py::class_<dd2md::DefectiveCrystalInterface>( m, "DefectiveCrystalInterface")
      .def( py::init([](const std::string& folderName){ // lambda function that returns an instantiation
               return std::unique_ptr< dd2md::DefectiveCrystalInterface>(
                     new dd2md::DefectiveCrystalInterface( folderName)
                     );
               }), py::arg("folderPath"))
      .def("readLmpStream",
            &dd2md::DisplacementFieldInterface::readLmpStream)
      .def("computeDisplacements",
            &dd2md::DisplacementFieldInterface::computeDisplacements)
         ;
}
/**********************************************************************/
VectorDim DisplacementFieldInterface::computeDisplacements()//const std::string& filename)//, FieldPointType& fieldPoints, const double& b_SI)
{
    std::cout<<"Computing DislocationDisplacement at field points..."<<std::flush;
   // TODO:
   // TODO:
   // TODO:
   // TODO:
   // TODO:
};

void DisplacementFieldInterface::readLmpStream(const std::string& filename)//, FieldPointType& fieldPoints, const double& b_SI)
{
    
    std::ifstream infile(filename.c_str(),std::ifstream::in);
    
    if(infile.is_open())
    {
        std::cout<<"Reading file "<<filename<<"..."<<std::flush;
        const auto t0= std::chrono::system_clock::now();
        size_t sizeA;
        std::string line;
        std::stringstream ss;
        
        std::getline(infile, line);
        std::getline(infile, line);
        std::getline(infile, line);
        // Get the total atom number
        ss<<line;
        ss>>sizeA;
        ss.clear(); ss.str("");
        
        // Skip lines
        for(int k=0; k<12; k++)
            std::getline(infile, line);
        
        // Read atom positions
        for(int k=1; k<=sizeA; k++)
        {
            std::getline(infile, line);

            ss<<line;
            int sID, aType;
            ss>>sID;
            ss>>aType;
            Eigen::Matrix<double,3,1> P;
            Eigen::Matrix<double,3,1> S(0.0,0.0,0.0); // S vector at the current
            ss>>P[0];
            ss>>P[1];
            ss>>P[2];
            ss.clear();
//            fieldPoints.emplace_back(sID,aType,P,S);
            this->fieldPoints.emplace_back(sID,P*1e-10/(this->b_SI));
        }
        infile.close();
        model::cout<<"  "<<sizeA<<" Atoms...";
        model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
    }
    return;
}
