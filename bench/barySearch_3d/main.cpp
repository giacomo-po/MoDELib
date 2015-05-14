
#include <model/Utilities/SequentialOutputFile.h>

#define _MODEL_BENCH_BARYSEARCH_
model::SequentialOutputFile<'S',true> searchFile;

#include <model/Mesh/SimplicialMesh.h>


using namespace model;

int main(int argc, char * argv[])
{
    Eigen::Matrix<double,3,1> P(Eigen::Matrix<double,3,1>::Zero());
    Eigen::Matrix<size_t,1,4> xID(Eigen::Matrix<size_t,1,4>::Zero());
    
    int meshID(0);
    
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }
    
    if (argc>3)
    {
        P<<atof(argv[2]),atof(argv[3]),atof(argv[4]);
    }
    
    xID<<atoi(argv[5]),atoi(argv[6]),atoi(argv[7]),atoi(argv[8]);
    
    SimplicialMesh<3> mesh;
    mesh.readMesh(meshID);
    
    searchFile<<P.transpose()<<" "<<Eigen::Matrix<double,1,4>::Zero()<<"\n";
    //    auto p=mesh2.search(P2);
    
    const Simplex<3,3>* guess=&(mesh.rbegin()->second);
    SimplicialMesh<3>::iterator sIter=mesh.find(xID);
    if(sIter!=mesh.end())
    {
        guess=&(sIter->second);
        std::cout<<"guess found"<<std::endl;
    }
    
    auto p=mesh.searchWithGuess(P,guess);

        std::cout<<"lastSearch= "<<p.second->xID<<std::endl;
    std::cout<<"Found? "<<p.first<<std::endl;
    std::cout<<"Boundary Simplex? "<<p.second->isBoundarySimplex()<<std::endl;

    if(!p.first)
    {// output boundary intersection
        
        const Eigen::Matrix<double,4,1> bary1(Eigen::Matrix<double,4,1>::Ones()/4.0);
        const Eigen::Matrix<double,4,1> bary2(p.second->pos2bary(P));
        
        int faceID;
        bary2.minCoeff(&faceID);
        
        searchFile<<p.second->bary2pos(p.second->faceLineIntersection(bary1,bary2,faceID)).transpose()<<" "
        /*      */<<p.second->xID<<std::endl;
    }
//
    SequentialOutputFile<'P',true> pFile;
    for (typename SimplexObserver<3,3>::const_iterator sIter=SimplexObserver<3,3>::simplexBegin();
         sIter!=SimplexObserver<3,3>::simplexEnd();++sIter)
    {
        pFile<<sIter->second->vertexPositionMatrix()<<"\n";
    }
    
    return 0;
}
