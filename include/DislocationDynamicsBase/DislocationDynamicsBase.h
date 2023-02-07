/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * Additional contributors:
 *   Nicholas Huebner Julian <njulian@lanl.gov>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DislocationDynamicsBase_H_
#define model_DislocationDynamicsBase_H_

#include <string>
#include <Eigen/Dense>

#include <DefectiveCrystalParameters.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
//#include <CrackSystem.h>
#include <DislocationDynamicsModule.h>
//#include <FEMnodeEvaluation.h>
//#include <DislocationNetwork.h>


namespace model
{
   template <int _dim> // corder isn't used until there's a DislocationNode
   class DislocationDynamicsBase
   {
      // desired members:
      //  those common to MicrostructureGenerator and DefectiveCrystal
      //  so that their members can be references to those in this class

      public:

         // constructor
         DislocationDynamicsBase( const std::string& folderName);

         static constexpr int dim=_dim;
         typedef DislocationDynamicsBase<dim> DislocationDynamicsBaseType;

         const std::string ddFolderName;
         DefectiveCrystalParameters simulationParameters;
         const std::set<int> periodicFaceIDs;
         const SimplicialMesh<dim> mesh;
         const Polycrystal<dim> poly;

   }; // class DislocationDynamicsBase

} // namespace model
#endif
