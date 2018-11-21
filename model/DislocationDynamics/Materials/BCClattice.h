/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_BCClattice_H_
#define model_BCClattice_H_

#include <vector>
#include <Eigen/Dense>

#include <LatticeMath.h>
#include <SlipSystem.h>

namespace model
{

    template<int dim>
    struct BCClattice
    {
    
    };

    template<>
    struct BCClattice<3> : public Lattice<3>
    {
        
        static constexpr int dim=3;
        static constexpr auto name="BCC";
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;

        
        BCClattice(const MatrixDim& Q) :
        /* init */ Lattice<dim>(getLatticeBasis(),Q)
        {
        
        }
    

        /**********************************************************************/
//        template <int dim>
        static Eigen::Matrix<double,dim,dim> getLatticeBasis()
        {/*!\returns The matrix of lattice vectors (in columns), in units of the
          * crystallographic Burgers vector.
          */
            
            Eigen::Matrix<double,dim,dim> temp;
            temp << -1.0,  1.0,  1.0,
            /*   */  1.0, -1.0,  1.0,
            /*   */  1.0,  1.0, -1.0;
            
            return temp/sqrt(3.0);
        }
        
        /**********************************************************************/
//        template <int dim>
        static std::vector<LatticePlaneBase> reciprocalPlaneNormals(const Lattice<dim>& lat)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            typedef Eigen::Matrix<long int,dim,1> VectorDimI;
            
            typedef LatticeVector<dim> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0),lat);
            LatticeVectorType a2(VectorDimI(0,1,0),lat);
            LatticeVectorType a3(VectorDimI(0,0,1),lat);
            LatticeVectorType  y(VectorDimI(1,1,1),lat);
            
            std::vector<LatticePlaneBase> temp;
            temp.emplace_back(a3,a1); // is ( 1, 0, 1) in cartesian
            temp.emplace_back( y,a2); // is ( 1, 0,-1) in cartesian
            temp.emplace_back(a2,a3); // is ( 0, 1, 1) in cartesian
            temp.emplace_back( y,a1); // is ( 0,-1, 1) in cartesian
            temp.emplace_back(a1,a2); // is ( 1, 1, 0) in cartesian
            temp.emplace_back( y,a3); // is (-1, 1, 0) in cartesian
            
            return temp;
        }

        /**********************************************************************/
//        template <int dim>
        static std::vector<SlipSystem> slipSystems(const Lattice<dim>& lat)
        {/*!\returns a std::vector of ReciprocalLatticeDirection(s) corresponding
          * the slip plane normals of the FCC lattice
          */
            typedef Eigen::Matrix<long int,3,1> VectorDimI;
            
            typedef LatticeVector<3> LatticeVectorType;
            LatticeVectorType a1(VectorDimI(1,0,0),lat);
            LatticeVectorType a2(VectorDimI(0,1,0),lat);
            LatticeVectorType a3(VectorDimI(0,0,1),lat);
            LatticeVectorType  y(VectorDimI(1,1,1),lat);
            
            std::vector<SlipSystem> temp;
            
            temp.emplace_back(a3,a1, a3); // is ( 1, 0, 1) in cartesian
            temp.emplace_back(a3,a1,a3*(-1)); // is ( 1, 0, 1) in cartesian
            temp.emplace_back(a3,a1, a1); // is ( 1, 0, 1) in cartesian
            temp.emplace_back(a3,a1,a1*(-1)); // is ( 1, 0, 1) in cartesian
            
            temp.emplace_back( y,a2, y); // is ( 1, 0,-1) in cartesian
            temp.emplace_back( y,a2,y*(-1)); // is ( 1, 0,-1) in cartesian
            temp.emplace_back( y,a2, a2); // is ( 1, 0,-1) in cartesian
            temp.emplace_back( y,a2,a2*(-1)); // is ( 1, 0,-1) in cartesian
            
            temp.emplace_back(a2,a3, a2); // is ( 0, 1, 1) in cartesian
            temp.emplace_back(a2,a3,a2*(-1)); // is ( 0, 1, 1) in cartesian
            temp.emplace_back(a2,a3, a3); // is ( 0, 1, 1) in cartesian
            temp.emplace_back(a2,a3,a3*(-1)); // is ( 0, 1, 1) in cartesian
            
            temp.emplace_back( y,a1, y); // is ( 0,-1, 1) in cartesian
            temp.emplace_back( y,a1,y*(-1)); // is ( 0,-1, 1) in cartesian
            temp.emplace_back( y,a1, a1); // is ( 0,-1, 1) in cartesian
            temp.emplace_back( y,a1,a1*(-1)); // is ( 0,-1, 1) in cartesian
            
            temp.emplace_back(a1,a2, a1); // is ( 1, 1, 0) in cartesian
            temp.emplace_back(a1,a2,a1*(-1)); // is ( 1, 1, 0) in cartesian
            temp.emplace_back(a1,a2, a2); // is ( 1, 1, 0) in cartesian
            temp.emplace_back(a1,a2,a2*(-1)); // is ( 1, 1, 0) in cartesian
            
            temp.emplace_back( y,a3, y); // is (-1, 1, 0) in cartesian
            temp.emplace_back( y,a3,y*(-1)); // is (-1, 1, 0) in cartesian
            temp.emplace_back( y,a3, a3); // is (-1, 1, 0) in cartesian
            temp.emplace_back( y,a3,a3*(-1)); // is (-1, 1, 0) in cartesian
            
            return temp;
        }

    };    

} // namespace model
#endif

