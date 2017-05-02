/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *                        Yinan Cui <cuiyinan@ucla.edu>
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

        int sizeDistributeType=0;
        double paramA;
        double paramB;
        double lmin;
        enum{Uniform=0,Gaussian=1,exponential=2,weibull=3};  
   

      //  std::uniform_real_distribution<double> sizeDistribution;
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

        
        MicrostructureGenerator() :
        /* init list */ generator(std::chrono::system_clock::now().time_since_epoch().count()),
        /* init list */ distribution(0.0,1.0),
        /* init list */ sizeDistributeType(0),
        /* init list */ paramA(0.0),
        /* init list */ paramB(1.0),
        /* init list */ lmin(1.0),
       // /* init list */ sizeDistribution(0.1,0.5),
        /* init list */ _minSize(0)
        {
            int meshID(0);
            EigenDataReader EDR;
            bool use_boundary=false;
            EDR.readScalarInFile("./DDinput.txt","use_boundary",use_boundary);
            
            EDR.readScalarInFile("./microstructureInput.txt","sizeDistributeType",sizeDistributeType);
            
             switch (sizeDistributeType)
		    {
		       case Uniform:
		           {
		           //paraA is minimum value, paraB is maximum value. ranges:[0 1]
            EDR.readScalarInFile("./microstructureInput.txt","unifA",paramA);
            EDR.readScalarInFile("./microstructureInput.txt","unifB",paramB);
		           }
		            break;
		       case Gaussian:
		            //For normal distribution, paraA=mean value (normlized by sample size along x direction) range: (0 1)
		            //                         paraB=standard deviation (normlized by mean value) range: >0
		           {
            EDR.readScalarInFile("./microstructureInput.txt","gaussA",paramA);
            EDR.readScalarInFile("./microstructureInput.txt","gaussB",paramB);
		           }
		            break;
		       case exponential:
		            //P(x|paraA)=paraA*exp(-paraA*x)  generally generate lots of small source
		           {
            EDR.readScalarInFile("./microstructureInput.txt","expA",paramA);
		           }
		            break;
		       case weibull:
		            //P(x|paraA,paraB)=(paraA/paraB)*(x/paraB)^(paraA-1)*exp(-[x/paraB]^paraA)
		            //paraA=shape parameter, paraB=scale parameter; generally generate lots of small source
		           {
            EDR.readScalarInFile("./microstructureInput.txt","weibA",paramA);
            EDR.readScalarInFile("./microstructureInput.txt","weibB",paramB);
		           }
		            break;
		    }
            
            
            if (use_boundary)
            {
                
                EDR.readScalarInFile("./DDinput.txt","meshID",meshID);
                mesh.readMesh(meshID);
                _minSize=min(mesh.xMax(0)-mesh.xMin(0),min(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2)));
                
            }
            else
            {
                assert(0 && "MICROSTRUCTURE GENERATION IN INFINITE SPACE NOT SUPPORTED YET.");
            }
            
            
            unsigned int materialZ;
            EDR.readScalarInFile("./DDinput.txt","material",materialZ); // material by atomic number Z
            Material<Isotropic>::select(materialZ);
            MatrixDimD C2Gtemp;
            EDR.readMatrixInFile("./DDinput.txt","C2G",C2Gtemp); // crystal-to-global orientation
            Material<Isotropic>::rotateCrystal(C2Gtemp);
            EDR.readScalarInFile("./DDinput.txt","Lmin",lmin); // material by atomic number Z
            

            
        }

        /**********************************************************************/
        LatticeVector<dim> randomPointInMesh()
        {
            LatticeVector<dim> L0=LatticeBase<dim>::snapToLattice(randomPoint());
            bool inside=mesh.search(L0.cartesian()).first;
            while(!inside)
            {
                L0=LatticeBase<dim>::snapToLattice(randomPoint());
                inside=mesh.search(L0.cartesian()).first;
            }
            return L0;
        }
    
        /**********************************************************************/
        double randomSize()
        {
            double randomvalue=0.0;
            while (_minSize*randomvalue<2.1*lmin)
            {
				switch (sizeDistributeType)
				{
				   case Uniform:
					   {
					   //paraA is minimum value, paraB is maximum value. ranges:[0 1]
						std::uniform_real_distribution<double> sizeDistribution(paramA,paramB);
						randomvalue=sizeDistribution(generator);
					   }
						break;
				   case Gaussian:
						//For normal distribution, paraA=mean value (normlized by sample size along x direction) range: (0 1)
						//                         paraB=standard deviation (normlized by mean value) range: >0
					   {
						std::normal_distribution<double> sizeDistribution(paramA,paramB*paramA);
						randomvalue=sizeDistribution(generator);
					   }
						break;
				   case exponential:
						//P(x|paraA)=paraA*exp(-paraA*x)  generally generate lots of small source
					   {
						std::exponential_distribution<double> sizeDistribution(paramA);
						randomvalue=sizeDistribution(generator);
					   }
						break;
				   case weibull:
						//P(x|paraA,paraB)=(paraA/paraB)*(x/paraB)^(paraA-1)*exp(-[x/paraB]^paraA)
						//paraA=shape parameter, paraB=scale parameter; generally generate lots of small source
					   {
						std::weibull_distribution<double> sizeDistribution(paramA,paramB);
						randomvalue=sizeDistribution(generator);
					   }
						break;
				}
            }
 
            return _minSize*randomvalue;
 
            
        }
        
    };

}
#endif
