#include <model/DislocationDynamics/DislocationNetwork.h>
#include <model/Utilities/SequentialOutputFile.h>

int main()
{
    
    DislocationNetwork<3,1,CatmullRom,16,UniformOpen> DN(argc,argv);

    // The Burgers vector of the dislocation
    Eigen::Matrix<double,3,1> B;
    B<<1.0,0.0,0.0;
    
    const double Lz=200.0;
    const int nz=50;
    const double dz=Lz/nz;
    const int np=2*nz+1;
    
    
    // Generate straight dislocation along z
    size_t oldID(0);
    for (int k=0;k<2*nz+1;++k)
    {
        Eigen::Matrix<double,3,np> P;
        P<<0.0,0.0,k*dz-Lz;
        const size_t newID(DN.insertVertex(P));
        if (k>0)
        {
           DN.connect(oldID,newID,B);
        }
        oldID=newID;
    }
    

    
    
    double Lx( 1.0e+03);
    double Ly( 1.0e+03);
    
    
    
    model::SequentialOutputFile<'A',1> analyticalFile;
    model::SequentialOutputFile<'N',1>  numericalFile;
    
    for (int i=-nL;i<nL+1;++i)
    {
        for (int j=-nL;j<nL+1;++j)
        {
            Eigen::Matrix<double,3,np> P;
            analyticalFile << P.transpose()<<"\n";
             numericalFile << P.transpose()<<"\n";
        }

    }
    
    

    
    return 0;
}
