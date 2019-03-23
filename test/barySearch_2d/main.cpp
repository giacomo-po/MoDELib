
#include <model/Utilities/SequentialOutputFile.h>

#define _MODEL_BENCH_BARYSEARCH_
model::SequentialOutputFile<'S',1> searchFile;

#include <model/Mesh/SimplicialMesh.h>


using namespace model;

int main(int argc, char * argv[])
{
    Eigen::Matrix<double,2,1> P2(Eigen::Matrix<double,2,1>::Zero());
    
    int meshID(0);
    
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }
    
    if (argc>3)
    {
        P2<<atof(argv[2]),atof(argv[3]);
    }
    
    SimplicialMesh<2> mesh2;
    mesh2.readMesh(meshID);
    
    searchFile<<P2.transpose()<<" "<<Eigen::Matrix<double,1,3>::Zero()<<"\n";    //    auto p=mesh2.search(P2);
    auto p=mesh2.searchWithGuess(P2,&(mesh2.rbegin()->second));
    
    std::cout<<"Found? "<<p.first<<std::endl;
    if(!p.first)
    {// output boundary intersection
        
        const Eigen::Matrix<double,3,1> bary1(Eigen::Matrix<double,3,1>::Ones()/3.0);
        const Eigen::Matrix<double,3,1> bary2(p.second->pos2bary(P2));
        
        int faceID;
        bary2.minCoeff(&faceID);
        
        searchFile<<p.second->bary2pos(p.second->faceLineIntersection(bary1,bary2,faceID)).transpose()<<" "
        /*      */<<p.second->xID<<std::endl;
    }
    
    SequentialOutputFile<'P',true> pFile;
    for (typename SimplexObserver<2,2>::const_iterator sIter=SimplexObserver<2,2>::simplexBegin();
         sIter!=SimplexObserver<2,2>::simplexEnd();++sIter)
    {
        pFile<<sIter->second->vertexPositionMatrix()<<"\n";
    }
    
    return 0;
}
