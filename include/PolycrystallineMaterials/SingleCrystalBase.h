/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SingleCrystalBase_H_
#define model_SingleCrystalBase_H_

#include <cmath>
#include <string>
#include <vector>
#include <tuple>

#include <TerminalColors.h> // defines mode::cout
#include <LatticeModule.h>
#include <SlipSystem.h>
#include <SecondPhase.h>


#include <TextFileParser.h>


namespace model
{
    
    template<int dim>
    struct SingleCrystalBase : public Lattice<dim>
    {
        
        typedef Lattice<dim> LatticeType;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef std::vector<std::shared_ptr<LatticePlaneBase>> PlaneNormalContainerType;
        typedef std::vector<std::shared_ptr<SlipSystem>> SlipSystemContainerType;
        typedef std::vector<std::shared_ptr<SecondPhase<dim>>> SecondPhaseContainerType;
        

        SingleCrystalBase(const MatrixDim& A,
                          const MatrixDim& C2G);
        
        virtual ~SingleCrystalBase();
        virtual const PlaneNormalContainerType& planeNormals() const =0;
        virtual const SlipSystemContainerType& slipSystems() const =0;
        virtual const SecondPhaseContainerType& secondPhases() const =0;
        
    };
    
    
}
#endif
