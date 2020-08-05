#include <iostream>
#include <DislocatedMaterial.h>
#include <StraightDislocationSegment.h>
#include <TextFileParser.h>
#include <vector>


using namespace model;


int main(int argc, char** argv)
{
    
    const bool computeSelfEnergiesOnly= argc > 1 ? atoi(argv[1]) : false;
    
    if(computeSelfEnergiesOnly)
    {
        const int nGP(TextFileParser("input.txt").readScalar<int>("nGP",false));
        
        const Eigen::MatrixXd segments(TextFileParser("input.txt").readMatrixCols<double>("S",9,false));
        std::ofstream pointsFile("points.txt");
        for(int r=0;r<segments.rows();++r)
        {
            const Eigen::Vector3d P0(segments.template block<1,3>(r,0));
            const Eigen::Vector3d P1(segments.template block<1,3>(r,3));
            const Eigen::Vector3d b(segments.template block<1,3>(r,3));
            const double length((P1-P0).norm());
            const Eigen::Vector3d t((P1-P0)/length);
            
            StraightDislocationSegment sds(P0,P1,b,length,t);
            
            const double Wns(DislocatedMaterial<3,Isotropic>::C2);
            
            const double dL(length/nGP);
            for(int k=0;k<nGP;++k)
            {
                const Eigen::Vector3d x(P0+(0.5+k)*dL*t);
                const double e(sds.elasticInteractionEnergy(x,t,b));
                pointsFile<<r<<" "<<k<<" "<<x.transpose()<<" "<<e<<" "<<dL<<" "<<Wns<<"\n";
            }
            
        }
    }
    else
    {
        
    }
    
    
    return 0;
}
