/* First define the non-singluar method used for calculations
 * _MODEL_NON_SINGULAR_DD_ = 0 classical theory
 * _MODEL_NON_SINGULAR_DD_ = 1 Cai's regularization method
 * _MODEL_NON_SINGULAR_DD_ = 2 Lazar's regularization method
 */
#define _MODEL_NON_SINGULAR_DD_ 0

#include <stdlib.h> // atoi
#include <chrono>
#include <model/LatticeMath/LatticeMath.h>
#include <model/DislocationDynamics/DislocationNetwork.h>
#include <model/Utilities/SequentialOutputFile.h>
#include <model/ParticleInteraction/SingleFieldPoint.h>
#include <model/DislocationDynamics/StressStraight.h>

using namespace model;

/******************************************************************************/
Eigen::Matrix<double,3,3> stressStraightEdge(const Eigen::Matrix<double,3,1>& P)
{/*!@param[in] P the field point
  * \returns the stress field of a straight edge dislocation at P
  */
    Eigen::Matrix<double,3,3> stress(Eigen::Matrix<double,3,3>::Zero()); // 3x3 matrix of zeros
    
    // From Hirth and Lothe ch 3.4, p76.
    // Assume unitary b and mu
    
    const double& x(P(0));
    const double& y(P(1));
    const double x2(x*x);
    const double y2(y*y);
    const double R4(std::pow(x2+y2,2));
    
    stress(0,0)=-y*(3.0*x2+y2)/R4;  // sigma_xx
    stress(1,1)= y*(    x2-y2)/R4;  // sigma_yy
    stress(2,2)=-2.0*Material<Isotropic>::nu*y/(x2+y2);  // sigma_zz
    stress(0,1)= x*(    x2-y2)/R4;  // sigma_xy
    stress(1,0)=stress(0,1);        // symmetry
    
    return stress/(2.0*M_PI*(1.0-Material<Isotropic>::nu));
}

/******************************************************************************/
int main(int argc, char * argv[])
{
    // Some definitions
    typedef  DislocationNetwork<3,1,CatmullRom,UniformOpen> DislocationNetworkType;
    typedef typename DislocationNetworkType::StressField StressField;
    
    // Create DislocationNetwork
    DislocationNetworkType DN(argc,argv);
    
    // The (unitary) Burgers vector of the dislocation. Length units are in b.
    Eigen::Matrix<double,3,1> b(1.0,0.0,0.0); // b along x1
    
    // Generate straight dislocation along z, centered at origin
    LatticeDirection<3> lz((Eigen::Matrix<double,3,1>()<<0.0,0.0,1.0).finished()); // shortest lattice vector along z
    const double dz=3.0*lz.cartesian().norm();  // spacing between vertices
    const int nz=50;        // dislocation has 2*nz+1 vertices
    const double Lz=dz*nz;  // dislocation extends along x3 from -Lz to Lz
    std::cout<<"Dislocation extends from z="<<-Lz<<" to z="<<Lz<<std::endl;
    
    size_t oldID(0);
    for (int k=-nz;k<nz+1;++k)
    {
        Eigen::Matrix<double,3,1> P(0.0,0.0,k*dz); // position of current vertex
        const size_t newID(DN.insertVertex(P).first->first); // insert vertex ind DislocationNetwork
        if (k>-nz)
        {
            DN.connect(oldID,newID,b); // connect two vertices with a DislocationSegment
        }
        oldID=newID;
    }
    DN.output(1);
    
    /**************************************************************************/
    // Construct a grid of field points on the plane z=0
    std::cout<<"Creating grid of FieldPoints..."<<std::flush;
    const auto t0= std::chrono::system_clock::now();
    typedef SingleFieldPoint<StressField> FieldPointType;
    std::deque<FieldPointType> fieldPoints; // the container of field points
    
    int nx=100;              // the number of grid points in x is 2*nx+1
    int ny=100;              // the number of grid points in y is 2*ny+1
    bool outputToFile(false);
    
    if (argc>1) // overwrite nx and ny using input arguments to program
    {
        ny=nx=atoi(argv[1]);
        if (argc>2)
        {
            ny=atoi(argv[2]);
            if (argc>3)
            {
                outputToFile=atoi(argv[3]);
            }
        }
    }
    
    double Lx( 1.0e+01);    // grid extends in [-Lx,Lx] along x direction (units of b)
    
    double Ly( 1.0e+01);    // grid extends in [-Ly,Ly] along y direction (units of b)
    
    for (int i=-nx;i<nx+1;++i)
    {
        for (int j=-ny;j<ny+1;++j)
        {
            Eigen::Matrix<double,3,1> P(i*Lx/nx,j*Ly/ny,0.0); // position of current field point
            fieldPoints.emplace_back(P,true);
        }
    }
    std::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
    
    
    /**************************************************************************/
    // Let the DislocationNetwork compute the stress at the field points
    std::cout<<"Computing DislocationStress at field points..."<<std::flush;
    const auto t1= std::chrono::system_clock::now();
    DN.updateQuadraturePoints();
    DN.computeField<FieldPointType,StressField>(fieldPoints);
    std::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<std::endl;
    
    /**************************************************************************/
    // Ouput results to file
    if(outputToFile)
    {
        std::cout<<"Outputing to file..."<<std::flush;
        const auto t2= std::chrono::system_clock::now();
        model::SequentialOutputFile<'S',1>  numericalFile;  // this is file S/S_0.txt
        model::SequentialOutputFile<'S',1>  analyticalFile; // this is file S/S_1.txt
        model::SequentialOutputFile<'S',1>  straightFile; // this is file S/S_2.txt
        
        Eigen::Matrix<double,3,1> P0(0.0,0.0,-Lz); // position of current vertex
        Eigen::Matrix<double,3,1> P1(0.0,0.0,+Lz); // position of current vertex

        StressStraight<3> stressStraight(P0,P1,b);
        
        for (unsigned int k=0;k<fieldPoints.size();++k)
        {
//            Eigen::Matrix<double,3,3> temp(Material<Isotropic>::C2*(fieldPoints[k].field<StressField>()+fieldPoints[k].field<StressField>().transpose()));
            Eigen::Matrix<double,3,3> temp(fieldPoints[k].field());
            
            
            numericalFile << fieldPoints[k].P.transpose()<<" "<<temp.row(0)
            /*                                        */ <<" "<<temp.row(1)
            /*                                        */ <<" "<<temp.row(2)<<"\n";
            
            Eigen::Matrix<double,3,3> s(stressStraightEdge(fieldPoints[k].P));
            analyticalFile << fieldPoints[k].P.transpose()<<" "<<s.row(0)
            /*                                         */ <<" "<<s.row(1)
            /*                                         */ <<" "<<s.row(2)<<"\n";
            
            Eigen::Matrix<double,3,3> sStraight(stressStraight.stress(fieldPoints[k].P));
            straightFile << fieldPoints[k].P.transpose()<<" "<<sStraight.row(0)
            /*                                       */ <<" "<<sStraight.row(1)
            /*                                       */ <<" "<<sStraight.row(2)<<"\n";
            
        }
        std::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]"<<std::endl;
    }
    return 0;
}

