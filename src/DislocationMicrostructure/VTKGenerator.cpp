/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VTKGenerator_cpp_
#define model_VTKGenerator_cpp_

#include <numbers>
#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationLoopNodeIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGenerator.h>
#include <VTKGenerator.h>
#include <PlaneLineIntersection.h>

namespace model
{

    VTKGenerator::VTKGenerator(const std::string& fileName) :
    /* init */ MicrostructureGeneratorBase(fileName)
    {
        
    }

    void VTKGenerator::generateDensity(MicrostructureGenerator& mg)
    {
        throw std::runtime_error("VTKGenerator only works for single entities.");
    }

    void VTKGenerator::generateIndividual(MicrostructureGenerator& mg)
    {
        const std::string vtkRelFile(this->parser.readString("vtkFile",false));
        const std::string vtkFile(std::filesystem::path(this->microstructureFileName).parent_path().string()+"/"+vtkRelFile);
        
        std::cout<<magentaBoldColor<<"Generating VTK microstructure from "<<vtkFile<<defaultColor<<std::endl;
        
        const std::string dislocationType(this->parser.readString("dislocationType",false));
        
        if(dislocationType=="loop")
        {
            
            
//            mg.insertJunctionLoop(loopNodePos,glidePlane,
//                               slipSystem.s.cartesian(),glidePlane->referencePlane->unitNormal,
//                                  glidePlane->referencePlane->snapToPlane(P0),grainID,DislocationLoopIO<dim>::GLISSILELOOP);

            
        }
//        else if(dislocationType=="dipole")
//        {
//            
//        }
        else
        {
            throw std::runtime_error("Unknown type "+dislocationType);
        }
    
    }

}
#endif
