#include <stdlib.h> // atoi
#include <model/DislocationDynamics/DislocationNetwork.h>
#include <model/Utilities/SequentialOutputFile.h>

using namespace model;

/******************************************************************************/
Eigen::Matrix<double,3,3> stressStraightEdge(const Eigen::Matrix<double,3,1>& P)
{
    Eigen::Matrix<double,3,3> stress(Eigen::Matrix<double,3,3>::Zero()); // 3x3 matrix of zeros
    
    // Assume unitary b and mu

    //  stress(0,0)=...; // sigma_xx
    //  stress(0,1)=...; // sigma_xy
    //  stress(1,1)=...; // sigma_yy

    //  stress(1,0)=stress(0,1); // sigma_yx is same as sigma_xy
    
    return stress;
}


/******************************************************************************/
struct SimpleFieldPoint :
/* inheritance */ public FieldPoint<SimpleFieldPoint,3,DislocationStress<3> >
{    
    const Eigen::Matrix<double,3,1> P;
    
    SimpleFieldPoint(const Eigen::Matrix<double,3,1>& Pin) : P(Pin){}
    
};


/******************************************************************************/
int main(int argc, char * argv[])
{
    // Some definitions
    typedef  DislocationNetwork<3,1,CatmullRom,16,UniformOpen> DislocationNetworkType;
    typedef typename DislocationNetworkType::StressField StressField;

    // Create DislocationNetwork
    DislocationNetworkType DN(argc,argv);

    // The (unitary) Burgers vector of the dislocation. Length units are in b.
    Eigen::Matrix<double,3,1> b(1.0,0.0,0.0);

    // Generate straight dislocation along z, centered at origin
    const double Lz=200.0;  // dislocation extends along x3 from -Lz to Lz
    const int nz=50;        // dislocation has 2*nz+1 vertices
    const double dz=Lz/nz;  // spacing between vertices
    
    size_t oldID(0);
    for (int k=0;k<2*nz+1;++k)
    {
        Eigen::Matrix<double,3,1> P(0.0,0.0,k*dz-Lz); // position of current vertex
        const size_t newID(DN.insertVertex(P)); // insert vertex ind DislocationNetwork
        if (k>0)
        {
           DN.connect(oldID,newID,b); // connect two vertices with a DislocationSegment 
        }
        oldID=newID;
    }
    
    // Construct a grid of field points on the plane z=0
    std::deque<SimpleFieldPoint> fieldPoints; // the container of field points
    
    int nx=500;              // the number of grid points in x is 2*nx+1
    int ny=500;              // the number of grid points in y is 2*ny+1

    if (argc>1) // overwrite nx and ny using input arguments to program
    {
        ny=nx=atoi(argv[1]);
        if (argc>2)
        {
            ny=atoi(argv[2]);
        }
    }
    
    double Lx( 1.0e+02);    // grid extends in [-Lx,Lx] along x direction (units of b)
    
    double Ly( 1.0e+02);    // grid extends in [-Ly,Ly] along y direction (units of b)

    for (int i=-nx;i<nx+1;++i)
    {
        for (int j=-ny;j<ny+1;++j)
        {
            Eigen::Matrix<double,3,1> P(i*Lx/nx,j*Ly/ny,0.0); // position of current field point
            fieldPoints.emplace_back(P);
        }
    }
    
    // Let the DislocationNetwork compute the stress at the field points
    DN.updateQuadraturePoints();
    DN.computeField<SimpleFieldPoint,StressField>(fieldPoints);

//    // Ouput results to file
//    model::SequentialOutputFile<'S',1>  numericalFile; // this is file S/S_0.txt
//    model::SequentialOutputFile<'S',1> analyticalFile; // this is file S/S_1.txt
//    
//    for (unsigned int k=0;k<fieldPoints.size();++k)
//    {
//        numericalFile << fieldPoints[k].P.transpose()<<" "<<fieldPoints[k].field<StressField>().row(0)
//        /*                                        */ <<" "<<fieldPoints[k].field<StressField>().row(1)
//        /*                                        */ <<" "<<fieldPoints[k].field<StressField>().row(2)<<"\n";
//        
//        Eigen::Matrix<double,3,3> temp(stressStraightEdge(fieldPoints[k].P));
//        analyticalFile << fieldPoints[k].P.transpose()<<" "<<temp.row(0)
//        /*                                         */ <<" "<<temp.row(1)
//        /*                                         */ <<" "<<temp.row(2)<<"\n";
//    }
    
    
    return 0;
}
