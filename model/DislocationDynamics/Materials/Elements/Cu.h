/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Element_Cu_H_
#define model_Element_Cu_H_

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
    struct PeriodicElement<29,Isotropic>
    {
        typedef FCC CrystalStructure;
        static constexpr int Z=29;
        static constexpr auto name="Copper";
        static constexpr double nu=0.34;                // Poisson ratio [-]
        static constexpr double mu=48e9;                // Shear modulus [Pa]
        static constexpr double b=0.2556e-9;            // Burgers vector[m]
        static constexpr double rho=8940.0;             // Mass density [kg/m^3]
        static constexpr double cs=sqrt(mu/rho);        // Shear wave speed [m/s]
        
        static constexpr DislocationMobility<FCC> dm=DislocationMobility<FCC>(b,mu,cs,3.3333e-07,3.3333e-07);
        
        static const std::deque<GrainBoundaryType<3>> grainBoundaryTypes;
    };
    
    const std::deque<GrainBoundaryType<3>> PeriodicElement<29,Isotropic>::grainBoundaryTypes=
    {
        GrainBoundaryType<3>("symmetric-tilt [100](20 1 0) sigma 401",
                             Vector3d(1,0,0),Vector3d(20,1,0),Vector3d(-20,1,0), // axis, normal1, normal2
                             473.0, // GB energy
                             20.041811162,	2 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](14 1 0) sigma 197",
                             Vector3d(1,0,0),Vector3d(14,1,0),Vector3d(-14,1,0), // axis, normal1, normal2
                             579.0, // GB energy
                             14.0377844517,	2// GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](10 1 0) sigma 101",
                             Vector3d(1,0,0),Vector3d(10,1,0),Vector3d(-10,1,0), // axis, normal1, normal2
                             686.0, // GB energy
                             10.0509161128,	2 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](8 1 0) sigma 65",
                             Vector3d(1,0,0),Vector3d(8,1,0),Vector3d(-8,1,0), // axis, normal1, normal2
                             757.0, // GB energy
                             8.0622761524,	2 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](6 1 0) sigma 37",
                             Vector3d(1,0,0),Vector3d(6,1,0),Vector3d(-6,1,0), // axis, normal1, normal2
                             838.0, // GB energy
                             6.0842421077,	2 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](5 1 0) sigma 13",
                             Vector3d(1,0,0),Vector3d(5,1,0),Vector3d(-5,1,0), // axis, normal1, normal2
                             878.0, // GB energy
                             3.3557157572,	1.316227766 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](9 2 0) sigma 85",
                             Vector3d(1,0,0),Vector3d(9,2,0),Vector3d(-9,2,0), // axis, normal1, normal2
                             910.0, // GB energy
                             3.0334710406,	1.316227766 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](4 1 0) sigma 17",
                             Vector3d(1,0,0),Vector3d(4,1,0),Vector3d(-4,1,0), // axis, normal1, normal2
                             914.0, // GB energy
                             2.7137086322,	1.316227766 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](7 2 0) sigma 53",
                             Vector3d(1,0,0),Vector3d(7,1,0),Vector3d(-7,1,0), // axis, normal1, normal2
                             939.0, // GB energy
                             2.3956286304,	1.316227766// GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](13 4 0) sigma 185",
                             Vector3d(1,0,0),Vector3d(13,4,0),Vector3d(-13,4,0), // axis, normal1, normal2
                             940.0, // GB energy
                             2.2375409042,	1.316227766 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](3 1 0) sigma 5",
                             Vector3d(1,0,0),Vector3d(3,1,0),Vector3d(-3,1,0), // axis, normal1, normal2
                             905.0, // GB energy
                             0.9999973204,	0.632455532 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](11 4 0) sigma 137",
                             Vector3d(1,0,0),Vector3d(11,4,0),Vector3d(-11,4,0), // axis, normal1, normal2
                             964.0, // GB energy
                             1.3085063159,	0.894427191 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](5 2 0) sigma 29. Spacing is HALVED",
                             Vector3d(1,0,0),Vector3d(5,2,0),Vector3d(-5,2,0), // axis, normal1, normal2
                             983.0, // GB energy
                             1.2042335191,	0.894427191 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](9 4 0) sigma 97",
                             Vector3d(1,0,0),Vector3d(9,4,0),Vector3d(-9,4,0), // axis, normal1, normal2
                             992.0, // GB energy
                             1.2456522298,	0.894427191 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](17 8 0) sigma 353. SPACING DIVIDED BY EIGHTEEN, AS PER ATOMISTIC STUDIES",
                             Vector3d(1,0,0),Vector3d(17,8,0),Vector3d(-17,8,0), // axis, normal1, normal2
                             984.0, // GB energy
                             1.3202348151,	0.894427191 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        
        
        
        
    
        GrainBoundaryType<3>("symmetric-tilt [100](2 1 0) sigma 5",
                             Vector3d(1,0,0),Vector3d(2,1,0),Vector3d(-2,1,0), // axis, normal1, normal2
                             951.0, // GB energy
                             1.4142097728,	0.894427191 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](15 8 0) sigma 289",
                             Vector3d(1,0,0),Vector3d(15,8,0),Vector3d(-15,8,0), // axis, normal1, normal2
                             941.0, // GB energy
                             1.9819847686,	1.1543203767// GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](7 4 0) sigma 65",
                             Vector3d(1,0,0),Vector3d(7,4,0),Vector3d(-7,4,0), // axis, normal1, normal2
                             901.0, // GB energy
                             2.1935628736,	1.1543203767 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](5 3 0) sigma 17. Spacing is HALVED",
                             Vector3d(1,0,0),Vector3d(5,3,0),Vector3d(-5,3,0), // axis, normal1, normal2
                             856.0, // GB energy
                             2.3798990201,	1.1543203767 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](13 8 0) sigma 233",
                             Vector3d(1,0,0),Vector3d(13,8,0),Vector3d(-13,8,0), // axis, normal1, normal2
                             853.0, // GB energy
                             2.4922915514,	1.1543203767 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](3 2 0) sigma 13",
                             Vector3d(1,0,0),Vector3d(3,2,0),Vector3d(-3,2,0), // axis, normal1, normal2
                             790.0, // GB energy
                             2.9429337209,	1.1543203767// GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](7 5 0) sigma 37",
                             Vector3d(1,0,0),Vector3d(7,5,0),Vector3d(-7,5,0), // axis, normal1, normal2
                             732.0, // GB energy
                             4.3022088527,	1.4142135624 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](4 3 0) sigma 25",
                             Vector3d(1,0,0),Vector3d(4,3,0),Vector3d(-4,3,0), // axis, normal1, normal2
                             677.0, // GB energy
                             5.0000625254,	1.4142135624  // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](5 4 0) sigma 41",
                             Vector3d(1,0,0),Vector3d(5,4,0),Vector3d(-5,4,0), // axis, normal1, normal2
                             595.0, // GB energy
                             6.4033171014,	1.4142135624 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](6 5 0) sigma 61",
                             Vector3d(1,0,0),Vector3d(6,5,0),Vector3d(-6,5,0), // axis, normal1, normal2
                             533.0, // GB energy
                             7.80939344,	1.4142135624  // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](7 6 0) sigma 85",
                             Vector3d(1,0,0),Vector3d(7,6,0),Vector3d(-7,6,0), // axis, normal1, normal2
                             484.0, // GB energy
                             9.216837054,	1.4142135624 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        
        GrainBoundaryType<3>("symmetric-tilt [100](8 7 0) sigma 113",
                             Vector3d(1,0,0),Vector3d(8,7,0),Vector3d(-8,7,0), // axis, normal1, normal2
                             677.0, // GB energy
                             10.6275717332,	1.4142135624  // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](9 8 0) sigma 145",
                             Vector3d(1,0,0),Vector3d(9,8,0),Vector3d(-9,8,0), // axis, normal1, normal2
                             595.0, // GB energy
                             12.0468155227,	1.4142135624 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("symmetric-tilt [100](10 9 0) sigma 181",
                             Vector3d(1,0,0),Vector3d(10,9,0),Vector3d(-10,9,0), // axis, normal1, normal2
                             533.0, // GB energy
                             13.4437604857,	1.4142135624 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        
        
        
        
        
        //GrainBoundaryType<3>("symmetric-tilt [100](210) sigma 5",
        //                     Vector3d(1,0,0),Vector3d(2,1,0),Vector3d(2,1,0), // axis, normal1, normal2
        //                     0.0, // GB energy
        //                     10.0,1 // GB dislocation spacing and Burgers vector MAGNITUDE
        //                     ),
        //GrainBoundaryType<3>("symmetric-tilt [100](310) sigma 5",
        //                     Vector3d(1,0,0),Vector3d(3,1,0),Vector3d(3,1,0), // axis, normal1, normal2
        //                     0.0, // GB energy
        //                     sqrt(5), 1 // GB dislocation spacing and Burgers vector MAGNITUDE
        //                     ),
        GrainBoundaryType<3>("symmetric-tilt [110](111) sigma 3",
                             Vector3d(1,1,0),Vector3d(1,1,1),Vector3d(1,1,1), // axis, normal1, normal2
                             0.0, // GB energy
                             10.0, 0 // GB dislocation spacing and Burgers vector MAGNITUDE
                             ),
        GrainBoundaryType<3>("asymmetric-tilt [110](1,1,1)(11,11,1) sigma ???",
                             Vector3d(1,1,0),Vector3d(1,1,1),Vector3d(11,11,1), // axis, normal1, normal2
                             0.0, // GB energy
                             30.0, 0 // GB dislocation spacing and Burgers vector MAGNITUDE
                             )
        
    };
    
} // namespace model
#endif
