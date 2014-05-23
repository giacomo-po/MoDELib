
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
    
    searchFile<<P2.transpose()<<"\n";
    auto p=mesh2.search(P2);
    std::cout<<"Found? "<<p.first<<std::endl;
    if(!p.first)
    {// output boundary intersection
        searchFile<<p.second->bary2pos(p.second->baryFaceIntersection(Eigen::Matrix<double,3,1>::Ones()/3.0,p.second->pos2bary(P2))).transpose()<<std::endl;
    }
    
    SequentialOutputFile<'P',true> pFile;
    for (typename SimplexObserver<2,2>::const_iterator sIter=SimplexObserver<2,2>::simplexBegin();
         sIter!=SimplexObserver<2,2>::simplexEnd();++sIter)
    {
        pFile<<sIter->second->vertexPositionMatrix()<<"\n";
    }
    
    return 0;
}
