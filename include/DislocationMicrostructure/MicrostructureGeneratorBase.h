/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureGeneratorBase_H_
#define model_MicrostructureGeneratorBase_H_


#include <string>
#include <Eigen/Dense>

#include <TextFileParser.h>
#include <LatticeModule.h>
#include <GlidePlaneModule.h>
#include <DislocationNodeIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationLoopNodeIO.h>



namespace model
{


    class MicrostructureGenerator;
    

    struct PolyPoint
    {
        
        std::shared_ptr<PeriodicPlanePatch<3>> periodicPlanePatch() const;
    

    
    };

    struct MicrostructureGeneratorBase
    {
        
        constexpr static int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef LatticeDirection<dim> LatticeDirectionType;
        typedef DislocationLoopIO<dim>::DislocationLoopType DislocationLoopType;
        typedef Eigen::Matrix<double,dim,dim>    MatrixDimD;
        typedef Eigen::Matrix<long int,dim,dim>    MatrixDimI;

        const std::string microstructureFileName;
        TextFileParser parser;
        const std::string type;
        const std::string style;
        const std::string tag;
//        const std::string microstructureType;
        
        
        MicrostructureGeneratorBase(const std::string&);
        
        virtual void generateIndividual(MicrostructureGenerator&) =0;
        virtual void generateDensity(MicrostructureGenerator&) =0;
        
    };





}
#endif
