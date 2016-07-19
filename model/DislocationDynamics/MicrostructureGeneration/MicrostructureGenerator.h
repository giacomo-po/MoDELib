/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureGenerator_H_
#define model_MicrostructureGenerator_H_

#include <chrono>
#include <random>
#include <assert.h>
#include <model/Utilities/EigenDataReader.h>
#include <model/Mesh/SimplicialMesh.h> // defines mode::cout
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/DislocationDynamics/Polycrystals/Polycrystal.h> // defines mode::cout

namespace model
{
    
    class MicrostructureGenerator
    {
        constexpr static int dim=3;
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        
        typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
        typedef Eigen::Matrix<long int,dim,dim>	MatrixDimI;
        
        std::mt19937 generator;
        std::uniform_real_distribution<double> distribution;
        std::uniform_real_distribution<double> sizeDistribution;
        double _minSize;
        
        /**********************************************************************/
        VectorDimD randomPoint()
        {
            VectorDimD P0;
            
            P0 << mesh.xMin(0)+distribution(generator)*(mesh.xMax(0)-mesh.xMin(0)),
            /* */ mesh.xMin(1)+distribution(generator)*(mesh.xMax(1)-mesh.xMin(1)),
            /* */ mesh.xMin(2)+distribution(generator)*(mesh.xMax(2)-mesh.xMin(2));
            
            return P0;
            
        }
        
        /**********************************************************************/
        static double min(const double& a,const double& b)
        {
            return a<b? a : b;
        }
        
    public:
        
        SimplicialMesh<dim> mesh;
        Polycrystal<dim> poly;
        
        MicrostructureGenerator() :
        /* init list */ generator(std::chrono::system_clock::now().time_since_epoch().count()),
        /* init list */ distribution(0.0,1.0),
        /* init list */ sizeDistribution(0.1,0.5),
        /* init list */ _minSize(0),
        /* init list */ poly(mesh)
        {
            int meshID(0);
            EigenDataReader EDR;
            bool use_boundary=false;
            EDR.readScalarInFile("./DDinput.txt","use_boundary",use_boundary);
            
            
            
            
            if (use_boundary)
            {
                
                EDR.readScalarInFile("./DDinput.txt","meshID",meshID);
                mesh.readMesh(meshID);
                _minSize=min(mesh.xMax(0)-mesh.xMin(0),min(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2)));
                
                poly.init("./DDinput.txt");
                
            }
            else
            {
                assert(0 && "MICROSTRUCTURE GENERATION IN INFINITE SPACE NOT SUPPORTED YET.");
            }
            
            
            unsigned int materialZ;
            EDR.readScalarInFile("./DDinput.txt","material",materialZ); // material by atomic number Z
            Material<Isotropic>::select(materialZ);
//            MatrixDimD C2Gtemp;
//            EDR.readMatrixInFile("./DDinput.txt","C2G",C2Gtemp); // crystal-to-global orientation
//            Material<Isotropic>::rotateCrystal(C2Gtemp);
            
            
        }
        
        /**********************************************************************/
        std::pair<LatticeVector<dim>,int> randomPointInMesh()
        {
            VectorDimD P0=randomPoint();
            auto searchResult=mesh.search(P0);
            if(searchResult.first)
            {// point inside
                LatticeVector<dim> L0 = poly.grain(searchResult.second->region->regionID).snapToLattice(P0);
                searchResult=mesh.searchRegionWithGuess(L0.cartesian(),searchResult.second);
                if(searchResult.first)
                {// point inside
                    return std::make_pair(L0,searchResult.second->region->regionID);
                }
                else
                {
                    return randomPointInMesh();
                }
            }
            else
            {
                return randomPointInMesh();
            }
            
            //            LatticeVector<dim> L0=LatticeBase<dim>::snapToLattice(randomPoint());
            //            bool inside=mesh.search(L0.cartesian()).first;
            //            while(!inside)
            //            {
            //                L0=LatticeBase<dim>::snapToLattice(randomPoint());
            //                inside=mesh.search(L0.cartesian()).first;
            //            }
            //            return L0;
        }
        
        /**********************************************************************/
        double randomSize()
        {
            return _minSize*sizeDistribution(generator);
        }
        
    };
    
}
#endif
