
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

namespace dd2md
{
template<int dim>
struct DisplacementFieldPoint
{
      typedef Eigen::Matrix<double,3,1> VectorDim;
      typedef Eigen::Matrix<double,3,3> MatrixDim;
      typedef model::DefectiveCrystal<3,0> DefectiveCrystalType;
    const int aID; // Atom ID
    const int aType; // Atom Type
    const Eigen::Matrix<double,dim,1> P; // Atom position in A
    Eigen::Matrix<double,dim,1> S; // Atom displacement in A

    DisplacementFieldPoint(const int& aIDin,
                           const int& aTypein,
                           const Eigen::Matrix<double,dim,1>& Pin,
                           Eigen::Matrix<double,dim,1>& Sin) :
    /* init */ aID(aIDin),
    /* init */ aType(aTypein),
    /* init */ P(Pin),
    /* init */ S(Sin)
    {}
}; // DisplacementFieldPoint


class DisplacementFieldInterface : public model::DefectiveCrystal<3,0>
{

   std::vector<FEMnodeEvaluation<LagrangeElement<3,2>,3,1>> fieldPoints;

   public:
   DisplacementFieldInterface(
         const std::string& lammpsData
         ) : DefectiveCrystalType( lammpsData)
   {}
   //std::string lammpsData( "A/Zr_pristine.lmp");

   template<typename FieldPointType>
   void readLmpStream(const std::string& filename);
   //, FieldPointType& fieldPoints, const double& b_SI);
   VectorDim computeDisplacements();//const std::string& filename);
}; // DisplacementFieldInterface 

/**********************************************************************/
} // d2md namespace
