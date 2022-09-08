/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SLIPSYSTEM_cpp_
#define model_SLIPSYSTEM_cpp_

#include <memory>
#include <assert.h>
#include <LatticeModule.h>
//#include <LatticePlaneBase.h>
//#include <LatticeVector.h>
//#include <RationalLatticeDirection.h>
#include <DislocationMobilityBase.h>
#include <SlipSystem.h>

namespace model
{

SlipSystem::SlipSystem(const LatticePlaneBase& n_in,
                       const LatticeVector<3>& slip_in,
                       const std::shared_ptr<DislocationMobilityBase>& mobility_in,
                       const std::shared_ptr<GammaSurface>& gammaSurface_in,
                       const std::shared_ptr<GlidePlaneNoise>& planeNoise_in):
    /* init */ n(n_in)
    /* init */,s(slip_in)
    /* init */,unitNormal(n.cartesian().normalized())
    /* init */,unitSlip(s.cartesian().normalized())
    /* init */,unitSlipFull(n.primitiveVectors.first.cartesian().normalized())
    /* init */,G2Lfull((MatrixDim()<<unitSlipFull,unitNormal.cross(unitSlipFull),unitNormal).finished())
    /* init */,mobility(mobility_in)
    /* init */,gammaSurface(gammaSurface_in)
    /* init */,planeNoise(planeNoise_in)
    {
        
        std::cout<<greenColor<<"Creating SlipSystem "<<this->sID<<defaultColor<<std::endl;
        std::cout<<"  s= "<<std::setprecision(15)<<std::scientific<<s.cartesian().transpose()<<std::endl;
        std::cout<<"  n= "<<std::setprecision(15)<<std::scientific<<n.cartesian().transpose()<<std::endl;
        std::cout<<"  mobility= "<<mobility->name<<std::endl;
        
        
        if(s.dot(n)!=0)
        {
            throw std::runtime_error("SlipSystem: PLANE NORMAL AND SLIP DIRECTION ARE NOT ORTHOGONAL.");
        }
        if(!mobility)
        {
            throw std::runtime_error("SlipSystem: MOBILITY CANNOT BE A NULLPTR.");
        }
    }

SlipSystem::SlipSystem(const LatticePlaneBase& n_in,
                       const RationalLatticeDirection<3>& slip_in,
                       const std::shared_ptr<DislocationMobilityBase>& mobility_in,
                       const std::shared_ptr<GammaSurface>& gammaSurface_in,
                       const std::shared_ptr<GlidePlaneNoise>& planeNoise_in):
    /* init */ n(n_in)
    /* init */,s(slip_in)
    /* init */,unitNormal(n.cartesian().normalized())
    /* init */,unitSlip(s.cartesian().normalized())
    /* init */,unitSlipFull(n.primitiveVectors.first.cartesian().normalized())
    /* init */,G2Lfull((MatrixDim()<<unitSlipFull,unitNormal.cross(unitSlipFull),unitNormal).finished())
    /* init */,mobility(mobility_in)
    /* init */,gammaSurface(gammaSurface_in)
    /* init */,planeNoise(planeNoise_in)
    {
        
        std::cout<<greenColor<<"Creating partial SlipSystem "<<this->sID<<defaultColor<<std::endl;
        std::cout<<"  s= "<<std::setprecision(15)<<std::scientific<<s.cartesian().transpose()<<std::endl;
        std::cout<<"  n= "<<std::setprecision(15)<<std::scientific<<n.cartesian().transpose()<<std::endl;
        std::cout<<"  mobility= "<<mobility->name<<std::endl;
        
        
        if(s.dot(n)!=0)
        {
            throw std::runtime_error("SlipSystem: PLANE NORMAL AND SLIP DIRECTION ARE NOT ORTHOGONAL.");
        }
        if(!mobility)
        {
            throw std::runtime_error("SlipSystem: MOBILITY CANNOT BE A NULLPTR.");
        }
    }

    bool SlipSystem::isPartial() const
    {
        return abs(s.rat.d)!=1;
    }

    bool SlipSystem::isSameAs(const RationalLatticeDirection<3>& s1,const ReciprocalLatticeDirection<3>& n1)
    {
        if(   ((s-s1).squaredNorm()==0 && (n-n1).squaredNorm()==0)
           || ((s+s1).squaredNorm()==0 && (n+n1).squaredNorm()==0)
           )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    double SlipSystem::misfitEnergy(const Eigen::Matrix<double,3,1>& b)
    {
        return gammaSurface? gammaSurface->operator()(b) : 0.0;
    }

Eigen::Matrix<double,2,1> SlipSystem::globalToLocal(const VectorDim& x) const
{
    return (G2Lfull*x).template segment<2>(0);
}

Eigen::Matrix<double,3,1> SlipSystem::localToGlobal(const Eigen::Matrix<double,2,1>& x) const
{
    return (G2Lfull.transpose()).template block<3,2>(0,0)*x;
}

std::tuple<Eigen::Matrix<double,3,3>,double,double> SlipSystem::gridInterp(const VectorDim& x, const bool Debugflag)
{   // Added by Hyunsoo (hyunsol@g.clemson.edu)
    
    if(planeNoise)
    {
        const std::tuple<double,double,double> gridNoise(planeNoise->gridInterp(globalToLocal(x),Debugflag));
        
        const Eigen::Matrix<double,3,3> ssStress( std::get<0>(gridNoise)*(G2Lfull.row(2).transpose()*G2Lfull.row(0)+G2Lfull.row(0).transpose()*G2Lfull.row(2))
                                                 +std::get<1>(gridNoise)*(G2Lfull.row(2).transpose()*G2Lfull.row(1)+G2Lfull.row(1).transpose()*G2Lfull.row(2)));
        const double ssRSS((ssStress*unitSlip).dot(unitNormal));
        return std::make_tuple(ssStress,ssRSS,std::get<2>(gridNoise));
    }
    else
    {
        return std::make_tuple(Eigen::Matrix<double,3,3>::Zero(),0.0,0.0);
    }
}

//std::tuple<Eigen::Matrix<double,3,3>,double,double> SlipSystem::gridInterp(const Eigen::Matrix<double,2,1>& localPos, const bool Debugflag)
//{   // Added by Hyunsoo (hyunsol@g.clemson.edu)
//
//    if(planeNoise)
//    {
//
//        /* ************ used parameters ****************
//         <x0,y1>(w2),Ind2     <x1,y1>(w3),Ind3
//         *------------------*
//         |   s2   |   s3     |
//         |        |(localPos)|
//         |________+_________ |
//         |        |          |
//         |        |          |
//         |   s0   |   s1     |
//         *------------------*
//         <x0,y0>(w0),Ind0     <x1,y0>(w1),Ind1
//         ______________________________________________
//         localPos(0), localPos(1) : gauss points (x_g, y_g)
//         x0,x1,y0,y1 : node coordinates
//         Ind0,Ind1,Ind2,Ind3 : local node indices
//         s0,s1,s2,s3 : fraction of areas
//         w0,w1,w2,w3 : weight
//         */
//
//        // Get the indices of the grid that contains the gauss point
//        // grid_num_x, grid_num_y always find "Ind3" indices due to the usage of "ceil" function)
////        const int global_num_x = ceil( localPos(0)/planeNoise->gridSpacing(0) ); // Mesh (grid) index, x/dx, grid_num_x is equivalent to "i" in GlidePlaneActor.cpp
////        const int global_num_y = ceil( localPos(1)/planeNoise->gridSpacing(1) ); // Mesh (grid) index, y/dx, grid_num_y is equivalent to "j" in GlidePlaneActor.cpp
//
////        const int grid_num_x = global_num_x-std::floor((global_num_x*1.0) /planeNoise->gridSize(0)) * planeNoise->gridSize(0); // periodic condition (0 != 255) & (0 = 256), ex) 256 - (256.0/256 * 256) = 0
////        const int grid_num_y = global_num_y-std::floor((global_num_y*1.0) /planeNoise->gridSize(1)) * planeNoise->gridSize(1); // periodic condition (0 != 255) & (0 = 256)
//
//        const auto global_num(planeNoise->gridIndex(localPos));
//        const int& global_num_x(global_num(0));
//        const int& global_num_y(global_num(1));
//        const auto grid_num(planeNoise->periodicGridIndex(global_num));
//        const int& grid_num_x(grid_num(0));
//        const int& grid_num_y(grid_num(1));
//
//        // if the periodic condition algorithm fails and finds the grid points outside of the plane, it stops the program
//        if (grid_num_x > planeNoise->gridSize(0)-1 || grid_num_y > planeNoise->gridSize(1)-1 || grid_num_x < 0 || grid_num_y < 0)
//        {
//            throw std::runtime_error("grid index is out of range!");
//        }
//        // By default, the previous_x is the grid point located right before the grid point that is found.
//        const int previous_x = (grid_num_x-1) >= 0 ? grid_num_x-1 : planeNoise->gridSize(0) - 1; // periodic condition (0 != 255) & (-1 = 255), if previous_x = -1, then previous_x = 255
//        const int previous_y = (grid_num_y-1) >= 0 ? grid_num_y-1 : planeNoise->gridSize(1) - 1;
//
//
//        // Convert the local indices to global indices;
////        const int Ind3 = planeNoise->gridSize(1)*planeNoise->gridSize(2)*grid_num_x+planeNoise->gridSize(2)*grid_num_y; // (i = 0 ~ 255, planeNoise->gridSize = 1 ~ 256), so i = grid_num_x =
////        const int Ind2 = planeNoise->gridSize(1)*planeNoise->gridSize(2)*previous_x+planeNoise->gridSize(2)*grid_num_y;
////        const int Ind1 = planeNoise->gridSize(1)*planeNoise->gridSize(2)*grid_num_x+planeNoise->gridSize(2)*previous_y;
////        const int Ind0 = planeNoise->gridSize(1)*planeNoise->gridSize(2)*previous_x+planeNoise->gridSize(2)*previous_y;     // gridsSize(1)*planeNoise->gridSize(2) -> decalre variable
//
//        const int Ind3 = planeNoise->linearIndex(grid_num_x,grid_num_y); // (i = 0 ~ 255, planeNoise->gridSize = 1 ~ 256), so i = grid_num_x =
//        const int Ind2 = planeNoise->linearIndex(previous_x,grid_num_y);
//        const int Ind1 = planeNoise->linearIndex(grid_num_x,previous_y);
//        const int Ind0 = planeNoise->linearIndex(previous_x,previous_y);     // gridsSize(1)*planeNoise->gridSize(2) -> decalre variable
//
//
//        const Eigen::Array<double,2,1> x1y1( (Eigen::Array<double,2,1>()<<global_num_x*planeNoise->gridSpacing(0), global_num_y*planeNoise->gridSpacing(1)).finished() ); // copy the matrix
//        const Eigen::Array<double,2,1> x1y0( (Eigen::Array<double,2,1>()<<global_num_x*planeNoise->gridSpacing(0), (global_num_y-1)*planeNoise->gridSpacing(1)).finished() );  // since it is periodic, it will give us correct area
//        const Eigen::Array<double,2,1> x0y1( (Eigen::Array<double,2,1>()<<(global_num_x-1)*planeNoise->gridSpacing(0), global_num_y*planeNoise->gridSpacing(1)).finished() );
//        const Eigen::Array<double,2,1> x0y0( (Eigen::Array<double,2,1>()<<(global_num_x-1)*planeNoise->gridSpacing(0), (global_num_y-1)*planeNoise->gridSpacing(1)).finished() );
//
//        // calculate the fraction areas and store the area values in a vector (s0, s1, s2, s3)
//        const Eigen::Array<double,4,1> fracArea( (Eigen::Array<double,4,1>()<<
//                                                  std::abs((localPos(0)-x0y0(0))*(localPos(1)-x0y0(1))), std::abs((localPos(0)-x1y0(0))*(localPos(1)-x1y0(1))),  // (localpos.array() - x0y0).prod()
//                                                  std::abs((localPos(0)-x0y1(0))*(localPos(1)-x0y1(1))), std::abs((localPos(0)-x1y1(0))*(localPos(1)-x1y1(1)))).finished() );
//
//        // get the weight; w0, w1, w2, w3
//        const double totalArea = planeNoise->gridSpacing(0)*planeNoise->gridSpacing(1);
//        const Eigen::Matrix<double,4,1> weight( (Eigen::Matrix<double,4,1>()<<fracArea(3)/totalArea, fracArea(2)/totalArea, fracArea(1)/totalArea, fracArea(0)/totalArea).finished() );
//
//        const Eigen::Matrix<double,4,1> sfNoise(planeNoise->stackingFault? (Eigen::Matrix<double,4,1>()<< planeNoise->stackingFault->operator[](Ind0), planeNoise->stackingFault->operator[](Ind1),
//                                                                            planeNoise->stackingFault->operator[](Ind2), planeNoise->stackingFault->operator[](Ind3)).finished() : Eigen::Matrix<double,4,1>::Zero());
//        const double effsfNoise(weight.dot(sfNoise));  //weight(0)*sfNoise(0) + weight(1)*sfNoise(1) + weight(2)*sfNoise(2) + weight(3)*sfNoise(3);
//
//        const Eigen::Matrix<double,4,1> solNoiseXZ(planeNoise->solidSolution? (Eigen::Matrix<double,4,1>()<< planeNoise->solidSolution->operator[](Ind0)(0), planeNoise->solidSolution->operator[](Ind1)(0),
//                                                                               planeNoise->solidSolution->operator[](Ind2)(0), planeNoise->solidSolution->operator[](Ind3)(0)).finished() : Eigen::Matrix<double,4,1>::Zero());
//        const Eigen::Matrix<double,4,1> solNoiseYZ(planeNoise->solidSolution? (Eigen::Matrix<double,4,1>()<< planeNoise->solidSolution->operator[](Ind0)(1), planeNoise->solidSolution->operator[](Ind1)(1),
//                                                                               planeNoise->solidSolution->operator[](Ind2)(1), planeNoise->solidSolution->operator[](Ind3)(1)).finished() : Eigen::Matrix<double,4,1>::Zero());
//        const double effsolNoiseXZ(weight.dot(solNoiseXZ));
//        const double effsolNoiseYZ(weight.dot(solNoiseYZ));
//
//
//        const Eigen::Matrix<double,3,3> ssStress( effsolNoiseXZ*(G2Lfull.row(2).transpose()*G2Lfull.row(0)+G2Lfull.row(0).transpose()*G2Lfull.row(2))
//                                                 +effsolNoiseYZ*(G2Lfull.row(2).transpose()*G2Lfull.row(1)+G2Lfull.row(1).transpose()*G2Lfull.row(2)));
//        const double ssRSS((ssStress*unitSlip).dot(unitNormal));
//
//
//        if (Debugflag)
//        {
//            std::cout << "grid_num_x = " << grid_num_x << std::endl
//            << "grid_num_y = " << grid_num_y << std::endl
//            << "previous_x = "  << previous_x << std::endl
//            << "previous_y = "  << previous_y << std::endl
//            << "Ind0 = " << Ind0 << std::endl
//            << "Ind1 = " << Ind1 << std::endl
//            << "Ind2 = " << Ind2 << std::endl
//            << "Ind3 = " << Ind3 << std::endl
//            << "x0y0 = " << x0y0.transpose() << std::endl
//            << "x1y0 = " << x1y0.transpose() << std::endl
//            << "x0y1 = " << x0y1.transpose() << std::endl
//            << "x1y1 = " << x1y1.transpose() << std::endl
//            << "fracArea = " << fracArea.transpose() << std::endl
//            << "weight = " << weight.transpose() << std::endl;
//            //                              << "sfNoise = " << sfNoise.transpose() << std::endl;
//        }
//        return std::make_tuple(ssStress,ssRSS,effsfNoise);
//    }
//    else
//    {
//        return std::make_tuple(Eigen::Matrix<double,3,3>::Zero(),0.0,0.0);
//    }
//
//
//
//}

}
#endif
