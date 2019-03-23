
#include <model/Utilities/SequentialOutputFile.h>
#include <model/Mesh/SimplicialMesh.h>


using namespace model;

int main(int argc, char * argv[])
{
    
    int meshID(0);
    
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }
    
    SimplicialMesh<3> mesh3;
    mesh3.readMesh(meshID);
    
    SequentialOutputFile<'P',true> pFile;
    SequentialOutputFile<'D',true> dFile;
    
    for (typename SimplexObserver<3,3>::const_iterator sIter=SimplexObserver<3,3>::simplexBegin();
         sIter!=SimplexObserver<3,3>::simplexEnd();++sIter)
    {
        pFile<<sIter->second->vertexPositionMatrix()<<"\n";
        dFile<<sIter->second->nda<<"\n";
    }
    
    return 0;
}
