
#include <iostream>
#include <UniformPeriodicGrid.h>



using namespace model;




int main(int argc, char * argv[])
{
    
    typedef typename UniformPeriodicGrid<2>::ArrayDimI Array2I;
    typedef typename UniformPeriodicGrid<2>::ArrayDimD Array2D;
    typedef typename UniformPeriodicGrid<3>::ArrayDimI Array3I;
    typedef typename UniformPeriodicGrid<3>::ArrayDimD Array3D;

    const Array2I gridSize(10,10);
    const Array2D gridSpacing(1.0,1.0);
    
    
    
    UniformPeriodicGrid<2> upg(gridSize,gridSpacing);

    const Array2D pos(10,5);
    const auto idx(upg.posToPeriodicIdx(pos));
    
    std::cout<<pos.transpose()<<", "<<idx.first.transpose()<<", "<<idx.second.transpose()<<std::endl;
    
    const auto idxAndWeights(upg.posToPeriodicCornerIdxAndWeights(pos));
    for(int p=0;p<idxAndWeights.first.size();++p)
    {
        std::cout<<idxAndWeights.first[p].transpose()<<"->"<<idxAndWeights.second[p]<<std::endl;
    }
    
//    const auto cornerIdx2(UniformPeriodicGrid<2>::cornerIdx(std::make_pair(Array2I(9,5),Array2I(0,6))));
//    for(const auto& x : cornerIdx2)
//    {
//        std::cout<<x.transpose()<<std::endl;
//    }
//    std::cout<<std::endl;
//
//    const auto cornerIdx3(UniformPeriodicGrid<3>::cornerIdx(std::make_pair(Array3I(9,5,4),Array3I(0,6,7))));
//    for(const auto& x : cornerIdx3)
//    {
//        std::cout<<x.transpose()<<std::endl;
//    }
    
//    std::cout<<std::endl;
    
//    int dim=2;
//    for(int p=0;p<std::pow(2,dim);++p)
//    {
//        for(int i=0;i<dim;++i)
//        {
//            const int ithIndex= (p%int(std::pow(2,dim-i)))/std::pow(2,dim-1-i);
//            std::cout<<ithIndex<<" ";
//        }
//        std::cout<<std::endl;
//    }
    
//    upg.cornerIdx(idx);
    
//    const auto idxAndWeights(upg.idxAndWeights(pos));
    


    return 0;
}

