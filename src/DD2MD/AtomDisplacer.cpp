/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_AtomDisplacer_cpp_
#define model_AtomDisplacer_cpp_

#include <AtomDisplacer.h>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace model
{

AtomDisplacer::AtomDisplacer(const std::string& folderName):
/* init */ ddBase(folderName)
/* init */,microstructureGenerator(ddBase)
/* init */,configFields(ddBase,microstructureGenerator.config())
{
    
}

const DDconfigIO<3>& AtomDisplacer::config() const
{
    return microstructureGenerator.config();
}

DDconfigIO<3>& AtomDisplacer::config()
{
    return microstructureGenerator.config();
}

void AtomDisplacer::readConfiguration(const size_t& runID)
{
    config().read(runID);
    configFields.updateConfiguration();
}

void AtomDisplacer::readMicrostructure()
{
    microstructureGenerator.readMicrostructureFile();
    config().print();
    configFields.updateConfiguration();
}

void AtomDisplacer::writeConfiguration(const size_t& runID)
{
    microstructureGenerator.writeConfigFiles(runID);
}

double AtomDisplacer::solidAngle(const VectorDim& x) const
{
    return configFields.solidAngle(x);
}

double AtomDisplacer::solidAngle(const double& x,const double& y,const double& z) const
{
    return solidAngle((VectorDim()<<x,y,z).finished());
}

typename AtomDisplacer::VectorDim AtomDisplacer::dislocationPlasticDisplacement(const VectorDim& x) const
{
    return configFields.dislocationPlasticDisplacement(x);
}

typename AtomDisplacer::VectorDim AtomDisplacer::dislocationPlasticDisplacement(const double& x,const double& y,const double& z) const
{
    return dislocationPlasticDisplacement((VectorDim()<<x,y,z).finished());
}

//pybind11::array_t<double, pybind11::array::c_style> AtomDisplacer::dislocationPlasticDisplacement(const double& x,const double& y,const double& z) const
//{
//    const VectorDim temp(dislocationPlasticDisplacement ((VectorDim()<<x,y,z).finished()));
//    pybind11::array_t<double, pybind11::array::c_style> outArray;
//    outArray[0] = temp[0];
//    outArray[1] = temp[1];
//    outArray[2] = temp[2];
//    return outArray;
//}

PYBIND11_MODULE(DD2MD,m)
{
    namespace py=pybind11;
    py::class_<model::AtomDisplacer>(m,"AtomDisplacer")
        .def(py::init<const std::string&>())
        .def("readMicrostructure", &model::AtomDisplacer::readMicrostructure)
        .def("readConfiguration", &model::AtomDisplacer::readConfiguration)
        .def("writeConfiguration", &model::AtomDisplacer::writeConfiguration)
        .def("solidAngle", static_cast<double (model::AtomDisplacer::*)(const double&,const double&,const double&) const>(&model::AtomDisplacer::solidAngle))
//        .def("dislocationPlasticDisplacement", static_cast<pybind11::array_t<double, pybind11::array::c_style> (model::AtomDisplacer::*)(const double&,const double&,const double&) const>(&model::AtomDisplacer::dislocationPlasticDisplacement))
        .def("dislocationPlasticDisplacement", static_cast<typename model::AtomDisplacer::VectorDim (model::AtomDisplacer::*)(const double&,const double&,const double&) const>(&model::AtomDisplacer::dislocationPlasticDisplacement))
    ;
}


//    PYBIND11_MODULE(DD2MD,m)
//    {
//        namespace py=pybind11;
//        py::class_<model::AtomDisplacer>(m,"AtomDisplacer")
//            .def(py::init([](const std::string& folderName)
//                          {
//                            return std::unique_ptr<AtomDisplacer>(new AtomDisplacer(folderName));
//
//                            }
//                          ),
//                 py::arg("folderName")
//                 )
//    //        .def(py::init([](const std::string& folderName){return DD2MD(folderName);}),
//    //             py::arg("folderName",none(false))
//    //             )
//        ;
//    }

}
#endif
