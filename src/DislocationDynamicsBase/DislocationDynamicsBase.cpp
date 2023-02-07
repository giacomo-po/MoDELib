/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * Additional contributors:
 *  Nicholas Huebner Julian <njulian@lanl.gov>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationDynamicsBase_CPP_
#define model_DislocationDynamicsBase_CPP_

#include <DislocationDynamicsBase.h>

namespace model
{
   template <int _dim>
   DislocationDynamicsBase<_dim>::DislocationDynamicsBase( const std::string& folderName) :
      /* init */ ddFolderName( folderName)
      /* init */,simulationParameters( ddFolderName)
      /* init */,periodicFaceIDs(
            TextFileParser(
               simulationParameters.traitsIO.polyFile
               ).template readSet<int>("periodicFaceIDs",true))
      /* init */,mesh(
            simulationParameters.traitsIO.meshFile,
            TextFileParser(
               simulationParameters.traitsIO.polyFile
               ).readMatrix<double>("A",3,3,true),
            TextFileParser(
               simulationParameters.traitsIO.polyFile
               ).readMatrix<double>("x0",1,3,true).transpose(),
            periodicFaceIDs
            )
        /* init */,poly(simulationParameters.traitsIO.polyFile,mesh)
   {
      if(!mesh.simplices().size())
      {
         throw std::runtime_error("Mesh is empty");
      }
   }

template class DislocationDynamicsBase <3>;

} // namespace model
#endif
