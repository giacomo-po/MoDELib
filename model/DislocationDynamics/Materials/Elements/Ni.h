/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Element_Ni_H_
#define model_Element_Ni_H_

#include <deque>
#include <Eigen/Dense>
#include <model/DislocationDynamics/Materials/MaterialSymmetry.h>
#include <model/DislocationDynamics/Materials/FCCcrystal.h>
#include <model/DislocationDynamics/Materials/BCCcrystal.h>
#include <model/DislocationDynamics/MobilityLaws/DislocationMobility.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundaryType.h>

namespace model
{
    using Eigen::Vector3d;
    
    /**************************************************************************/
    /**************************************************************************/
    template <>
    struct PeriodicElement<28,Isotropic>
    {
        typedef FCC CrystalStructure;
        static constexpr int Z=28;
        static constexpr auto name="Nickel";
        static constexpr double nu=0.31;                // Poisson ratio [-]
        static constexpr double mu=76e9;                // Shear modulus [Pa]
        static constexpr double b=0.2489e-9;            // Burgers vector[m]
        static constexpr double rho=8908.0;             // Mass density [kg/m^3]
        static constexpr double cs=sqrt(mu/rho);        // Shear wave speed [m/s]
        
        //! FCC mobility law with data from Olmsted MSMSE 13(3), 2005.
        static constexpr DislocationMobility<FCC> dm=DislocationMobility<FCC>(b,mu,cs,5.0e-08,6.4e-08);
    };
    
} // namespace model
#endif
