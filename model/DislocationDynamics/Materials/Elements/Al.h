/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Element_Al_H_
#define model_Element_Al_H_

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
    struct PeriodicElement<13,Isotropic>
    {
        typedef FCC CrystalStructure;
        static constexpr int Z=13;
        static constexpr auto name="Aluminum";
        static constexpr double nu=0.347;               // Poisson ratio [-]
        static constexpr double mu=26.0e9;              // Shear modulus [Pa]
        static constexpr double b=0.2851e-9;            // Burgers vector[m]
        static constexpr double rho=2700.0;             // Mass density [kg/m^3]
        static constexpr double cs=sqrt(mu/rho);        // Shear wave speed [m/s]
        
        //! FCC-mobility law with data from Olmsted MSMSE 13(3), 2005.
        static constexpr DislocationMobility<FCC> dm=DislocationMobility<FCC>(b,mu,cs,3.9e-08,7.5e-08);
        
        static const std::deque<GrainBoundaryType<3>> grainBoundaryTypes;

    };
    
    const std::deque<GrainBoundaryType<3>> PeriodicElement<13,Isotropic>::grainBoundaryTypes=
    {
        
        
        //! GB energy values obtained with LAMMPs data from Ercolessi and Adams, MSMSE 12(665-670), 2004
        
        
        
        
        //These are the "Special" grain boundaries - low angle GBs and high angle GBs : Defined by a total number of 2 SUM units per nominal periodic dimension.
        
        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LOW ANGLE 0-18.4 degrees   START           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        GrainBoundaryType<3>("symmetric-tilt [100](20 1 0) sigma 401",
                             Vector3d(1,0,0),Vector3d(20,1,0),Vector3d(-20,1,0), // axis, normal1, normal2
                             335.0, // GB energy
                             14.15980226,1.4142135623731 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](14 1 0) sigma 197",
                             Vector3d(1,0,0),Vector3d(14,1,0),Vector3d(-14,1,0), // axis, normal1, normal2
                             393.0, // GB energy
                             9.924723808,	sqrt(2)// GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](10 1 0) sigma 101",
                             Vector3d(1,0,0),Vector3d(10,1,0),Vector3d(-10,1,0), // axis, normal1, normal2
                             448.0, // GB energy
                             7.106335854,	sqrt(2) // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](8 1 0) sigma 65",
                             Vector3d(1,0,0),Vector3d(8,1,0),Vector3d(-8,1,0), // axis, normal1, normal2
                             481.0, // GB energy
                             5.700869464,	sqrt(2) // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](6 1 0) sigma 37",
                             Vector3d(1,0,0),Vector3d(6,1,0),Vector3d(-6,1,0), // axis, normal1, normal2
                             514.0, // GB energy
                             4.301164685,	sqrt(2) // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LOW ANGLE 0-18.4 degrees     END           ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
        
        /*~~~~~~~~~~~~~~~~~~~~~~ 'Special' High angle 18.4-73 degrees START ~~~~~~~~*/
        GrainBoundaryType<3>("symmetric-tilt [100](5 1 0) sigma 13",
                             Vector3d(1,0,0),Vector3d(5,1,0),Vector3d(-5,1,0), // axis, normal1, normal2
                             521.0, // GB energy
                             3.605546207,1.4142135623731 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](4 1 0) sigma 17",
                             Vector3d(1,0,0),Vector3d(4,1,0),Vector3d(-4,1,0), // axis, normal1, normal2
                             517.0, // GB energy
                             2.915475947,1.4142135623731 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](3 1 0) sigma 5",
                             Vector3d(1,0,0),Vector3d(3,1,0),Vector3d(-3,1,0), // axis, normal1, normal2
                             491.0, // GB energy
                             2.236067977,1.4142135623731 // GB dislocation spacing and Burgers vector
                             ),
        
        
        //@@@@@@@@@@   Angles larger than 45 degrees below: Burgers vector from sqrt(2) to 1.0  @@@@@
        GrainBoundaryType<3>("symmetric-tilt [100](2 1 0) sigma 5",
                             Vector3d(1,0,0),Vector3d(2,1,0),Vector3d(-2,1,0), // axis, normal1, normal2
                             459.0, // GB energy
                             1.58113883,1 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](3 2 0) sigma 13",
                             Vector3d(1,0,0),Vector3d(3,2,0),Vector3d(-3,2,0), // axis, normal1, normal2
                             465.0, // GB energy
                             2.549509757,1// GB dislocation spacing and Burgers vector
                             ),
        
        /*~~~~~~~~~~~~~~~~~~~~~~ 'Special' High angle 18.4-73 degrees END ~~~~~~~~*/
        /*~~~~~~~~~ LOW ANGLE 73-90 degrees START   ~~~~~~~~~~~~~~~~~~~~*/
        GrainBoundaryType<3>("symmetric-tilt [100](4 3 0) sigma 25",
                             Vector3d(1,0,0),Vector3d(4,3,0),Vector3d(-4,3,0), // axis, normal1, normal2
                             430.0, // GB energy
                             3.535533906,1  // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](5 4 0) sigma 41",
                             Vector3d(1,0,0),Vector3d(5,4,0),Vector3d(-5,4,0), // axis, normal1, normal2
                             395.0, // GB energy
                             4.527692569,1 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](6 5 0) sigma 61",
                             Vector3d(1,0,0),Vector3d(6,5,0),Vector3d(-6,5,0), // axis, normal1, normal2
                             364.0, // GB energy
                             5.522680509,1  // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](7 6 0) sigma 85",
                             Vector3d(1,0,0),Vector3d(7,6,0),Vector3d(-7,6,0), // axis, normal1, normal2
                             337.0, // GB energy
                             6.519202405,1 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](8 7 0) sigma 113",
                             Vector3d(1,0,0),Vector3d(8,7,0),Vector3d(-8,7,0), // axis, normal1, normal2
                             315.0, // GB energy
                             7.516648189,1  // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](9 8 0) sigma 145",
                             Vector3d(1,0,0),Vector3d(9,8,0),Vector3d(-9,8,0), // axis, normal1, normal2
                             296.0, // GB energy
                             8.514693183,1 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](10 9 0) sigma 181",
                             Vector3d(1,0,0),Vector3d(10,9,0),Vector3d(-10,9,0), // axis, normal1, normal2
                             280.0, // GB energy
                             9.513148795,1 // GB dislocation spacing and Burgers vector
                             ),
        /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ LOW ANGLE 73-90 degrees END ~~~*/
        
        
        
        //These are the "General" high angle symmetric tilt <001> grain boundaries
        GrainBoundaryType<3>("symmetric-tilt [100](9 2 0) sigma 85",
                             Vector3d(1,0,0),Vector3d(9,2,0),Vector3d(-9,2,0), // axis, normal1, normal2
                             535.0, // GB energy
                             3.259605126,1.41439574 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](7 2 0) sigma 53",
                             Vector3d(1,0,0),Vector3d(7,2,0),Vector3d(-7,2,0), // axis, normal1, normal2
                             531.0, // GB energy
                             2.57390753,1.41439574// GB dislocation spacing and Burgers vector
                             ),
        
        GrainBoundaryType<3>("symmetric-tilt [100](13 4 0) sigma 185",
                             Vector3d(1,0,0),Vector3d(13,4,0),Vector3d(-13,4,0), // axis, normal1, normal2
                             534.0, // GB energy
                             2.40442301,1.41439574 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](11 4 0) sigma 137",
                             Vector3d(1,0,0),Vector3d(11,4,0),Vector3d(-11,4,0), // axis, normal1, normal2
                             532.0, // GB energy
                             2.0691181,1.41439574 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](5 2 0) sigma 53",
                             Vector3d(1,0,0),Vector3d(5,2,0),Vector3d(-5,2,0), // axis, normal1, normal2
                             526.0, // GB energy
                             1.90394327,1.41439574 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](9 4 0) sigma 53",
                             Vector3d(1,0,0),Vector3d(9,4,0),Vector3d(-9,4,0), // axis, normal1, normal2
                             537.0, // GB energy
                             1.3928388,1 // GB dislocation spacing and Burgers vector
                             ),
        
        
        GrainBoundaryType<3>("symmetric-tilt [100](17 8 0) sigma 353.",
                             Vector3d(1,0,0),Vector3d(17,8,0),Vector3d(-17,8,0), // axis, normal1, normal2
                             533.0, // GB energy
                             1.3202348151,1 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](15 8 0) sigma 289",
                             Vector3d(1,0,0),Vector3d(15,8,0),Vector3d(-15,8,0), // axis, normal1, normal2
                             522.0, // GB energy
                             1.717259326,1// GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](7 4 0) sigma 65",
                             Vector3d(1,0,0),Vector3d(7,4,0),Vector3d(-7,4,0), // axis, normal1, normal2
                             513.0, // GB energy
                             1.900292375,1 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](5 3 0) sigma 17. Spacing is HALVED",
                             Vector3d(1,0,0),Vector3d(5,3,0),Vector3d(-5,3,0), // axis, normal1, normal2
                             498.0, // GB energy
                             2.061552813,1 // GB dislocation spacing and Burgers vector
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](13 8 0) sigma 233",
                             Vector3d(1,0,0),Vector3d(13,8,0),Vector3d(-13,8,0), // axis, normal1, normal2
                             511.0, // GB energy
                             2.158703314,1 // GB dislocation spacing and Burgers vector
                             ),
        
        GrainBoundaryType<3>("symmetric-tilt [100](7 5 0) sigma 37",
                             Vector3d(1,0,0),Vector3d(7,5,0),Vector3d(-7,5,0), // axis, normal1, normal2
                             461.0, // GB energy
                             3.041381265,1 // GB dislocation spacing and Burgers vector
                             ),
        
        //THESE ARE EXPERIMENTAL BENCHMARK "general" asymmetric tilt GBS
        
        //High angle case:
        GrainBoundaryType<3>("Experiment asymmetric tilt-high angle",
                             Vector3d(0.3144011287936295806, -0.7820602669121764494,  0.5380833291698851051),Vector3d(-8.485281374238571317,  7.071067811865476394,  7.071067811865476394),Vector3d(-0.707106781186548683, -3.535533905932737753, -2.121320343559642385), // axis, normal1, normal2
                             9999999.0, // GB energy
                             9999999.0, 1 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("Experiment asymmetric tilt-low angle",
                             Vector3d(0.707107,   -0.707107, 2.14345e-15),Vector3d(8,8,7),Vector3d(8,8,9), // axis, normal1, normal2
                             9999999.0, // GB energy
                             9999999.0, 1 // GB dislocation spacing and Burgers vector MAGNITUDE
                             )
        
        
    };
    
} // namespace model
#endif
