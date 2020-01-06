/* First define the non-singluar method used for calculations
 * _MODEL_NON_SINGULAR_DD_ = 0 classical theory
 * _MODEL_NON_SINGULAR_DD_ = 1 Cai's regularization method
 * _MODEL_NON_SINGULAR_DD_ = 2 Lazar's regularization method
 */
#define _MODEL_NON_SINGULAR_DD_ 0

#include <stdlib.h> // atoi
#include <math.h>
#include <chrono>
#include <model/DislocationDynamics/DislocationNetwork.h>
#include <model/Utilities/SequentialOutputFile.h>
#include <model/ParticleInteraction/FieldPoint.h>
#include <model/DislocationDynamics/StressStraight.h>

using namespace model;


/******************************************************************************/
template<typename Field>
struct DisplacementFieldPoint :
/* inheritance */ public FieldPoint<DisplacementFieldPoint<Field>,Field::dim,Field>
{
    
    typedef FieldPoint<DisplacementFieldPoint<Field>,Field::dim,Field> FieldPointType;
    const Eigen::Matrix<double,Field::dim,1> P;
    const Eigen::Matrix<double,Field::dim,1> S;
    
    DisplacementFieldPoint(const Eigen::Matrix<double,Field::dim,1>& Pin,
                           const Eigen::Matrix<double,Field::dim,1>& Sin) :
    /* init */ FieldPointType(true),
    /* init */ P(Pin),
    /* init */ S(Sin)
    {}
    
};

/******************************************************************************/
Eigen::Matrix<double,3,1> displacementStraightEdge(const Eigen::Matrix<double,3,1>& P)
{/*!@param[in] P the field point
  * \returns the displacement field at P
  *
  *
  */
    Eigen::Matrix<double,3,1> disp(Eigen::Matrix<double,3,1>::Zero()); // 3x3 matrix of zeros
    
    const double& x(P(0));
    const double& y(P(1));
    const double  x2(x*x);
    const double  y2(y*y);
    const double  R2(x2+y2);
    const double& nu(Material<Isotropic>::nu);
    
    double at2yx=atan2(y,x); // atan2 output results in [-pi;pi]
    if(at2yx<0.0)
    {// transform at2yx in [0;2*pi]
        at2yx+=2.0*M_PI;
    }
    at2yx-=M_PI; // divide equally jump across slip plane
    
    // From Hirth and Lothe, Therory of Dislocations 2nd Ed, p 78. Eq. (3-45) and (3.56)
    //    disp(0)= at2yx + x*y/(2.0*(1.0-nu)*R2);
    //    disp(1)= -((1.0-2.0*nu)/(4.0*(1.0-nu))*log(R2)+(x2-y2)/(4.0*(1.0-nu)*R2));
    
    // Eshelby (Brit J Appl Phys 1966 v.17)
    disp(0)= at2yx - x*y/(2.0*(1.0-nu)*R2);
    disp(1)= -(1.0-2.0*nu)/(4.0*(1.0-nu))*log(R2)+y2/(2.0*(1.0-nu)*R2);
    //    disp(1)= -(1.0-2.0*nu)/(4.0*(1.0-nu))*log(x2+y2)+(x2-y2)/(4.0*(1.0-nu)*R2);
    disp *= 1.0/(2.0*M_PI);
    
    //From Nabarro, p.57ttttt
    return disp;
}

/******************************************************************************/
int main(int argc, char * argv[])
{
    // Some definitions
    typedef  DislocationNetwork<3,1,CatmullRom,UniformOpen> DislocationNetworkType;
    typedef typename DislocationNetworkType::DisplacementField DisplacementField;
    
    // Create DislocationNetwork
    DislocationNetworkType DN(argc,argv);
    
    // The (unitary) Burgers vector of the dislocation. Length units are in b.
    Eigen::Matrix<double,3,1> b(1.0,0.0,0.0);
    
    // Generate straight dislocation along z, centered at origin
    const double Lz=200.0;  // dislocation extends along x3 from -Lz to Lz
    const double Lx=200.0;  // dislocation extends along x3 from -Lz to Lz
    
    Eigen::Matrix<double,3,1> P0(0.0,0.0,-Lz/2.0);
    Eigen::Matrix<double,3,1> P1(0.0,0.0, Lz/2.0);
    Eigen::Matrix<double,3,1> P2( Lx,0.0, Lz/2.0);
    Eigen::Matrix<double,3,1> P3( Lx,0.0,-Lz/2.0);
    
    // Connect four straight segments (do no create a loop to avoid curvature)
    const size_t n0(DN.insertVertex(P0).first->first); // insert vertex ind DislocationNetwork
    const size_t n1(DN.insertVertex(P1).first->first); // insert vertex ind DislocationNetwork
    DN.connect(n0,n1,b); // connect two vertices with a DislocationSegment
    
    const size_t n2(DN.insertVertex(P1).first->first); // insert vertex ind DislocationNetwork
    const size_t n3(DN.insertVertex(P2).first->first); // insert vertex ind DislocationNetwork
    DN.connect(n2,n3,b); // connect two vertices with a DislocationSegment
    
    const size_t n4(DN.insertVertex(P2).first->first); // insert vertex ind DislocationNetwork
    const size_t n5(DN.insertVertex(P3).first->first); // insert vertex ind DislocationNetwork
    DN.connect(n4,n5,b); // connect two vertices with a DislocationSegment
    
    const size_t n6(DN.insertVertex(P3).first->first); // insert vertex ind DislocationNetwork
    const size_t n7(DN.insertVertex(P0).first->first); // insert vertex ind DislocationNetwork
    DN.connect(n6,n7,b); // connect two vertices with a DislocationSegment
    
    // Reach a stable number of nodes
    size_t nodeOrder=0;
    while(nodeOrder!=DN.nodeOrder())
    {
        nodeOrder=DN.nodeOrder();
        DN.remesh();
    }
    
    //    DN.runSteps();
    
    /**************************************************************************/
    // Construct a grid of field points on the plane z=0
    std::cout<<"Creating grid of FieldPoints..."<<std::flush;
    const auto t0= std::chrono::system_clock::now();
    typedef DisplacementFieldPoint<DisplacementField> FieldPointType;
    std::deque<FieldPointType> fieldPoints; // the container of field points
    int nx=10;              // the number of grid points in x is 2*nx+1
    int ny=10;              // the number of grid points in y is 2*ny+1
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
    
    const double Gx( nx);    // grid extends in [-Lx,Lx] along x direction (units of b)
    const double Gy( ny);    // grid extends in [-Ly,Ly] along y direction (units of b)
    for (int i=-nx;i<nx+1;++i)
    {
        for (int j=-ny;j<ny+1;++j)
        {
            Eigen::Matrix<double,3,1> P(i*Gx/nx,j*Gy/ny,0.0); // position of current field point
            Eigen::Matrix<double,3,1> S(0.0,1.0,0.0); // S vector at the current field point, constructed so that it does not intersect the slip surface
            if(P(1)<0.0)
            {
                S*=-1.0;
            }
            fieldPoints.emplace_back(P,S);
        }
    }
    std::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
    
    /**************************************************************************/
    // Let the DislocationNetwork compute the stress at the field points
    DN.updateQuadraturePoints();
    std::cout<<"Computing DislocationDisplacement at field points..."<<std::flush;
    const auto t1= std::chrono::system_clock::now();
    DN.computeField<FieldPointType,DisplacementField>(fieldPoints);
    std::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<std::endl;
    
    /**************************************************************************/
    // Ouput results to file
    if(outputToFile)
    {
        std::cout<<"Outputing to file..."<<std::flush;
        const auto t2= std::chrono::system_clock::now();
        model::SequentialOutputFile<'D',1>  numericalFile;  // this is file D/D_0.txt
        model::SequentialOutputFile<'D',1>  analyticalFile; // this is file D/D_1.txt
        
        for (unsigned int k=0;k<fieldPoints.size();++k)
        {
            numericalFile << fieldPoints[k].P.transpose()<<" "<<fieldPoints[k].field<DisplacementField>().transpose()<<"\n";
            analyticalFile << fieldPoints[k].P.transpose()<<" "<<displacementStraightEdge(fieldPoints[k].P).transpose()<<"\n";
        }
        std::cout<<" done.["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]"<<std::endl;
    }
    return 0;
}

