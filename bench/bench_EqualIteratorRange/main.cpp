
#include <iostream>
#include <iomanip>
#include <map>
#include <math.h>
#include <chrono>
#include <omp.h>
#include <Model/Threads/EqualIteratorRange.h>



using namespace model;

int main (int argc, char * const argv[])
{
    
    typedef std::map<int,int> MapType;
    
    int N=25;
    if (argc>1)
    {
        N=atoi(argv[1]);
    }
    
    MapType v;
    for(int k=0;k<N;++k)
    {
        v.emplace(k,0);
    }
    
    
    auto t1= std::chrono::system_clock::now();
    
    const size_t nThreads = omp_get_max_threads();
    EqualIteratorRange<MapType > eir(v.begin(),v.end(),nThreads);
    
#pragma omp parallel for
    for (unsigned int k = 0; k < eir.size(); ++k)
    {
        for (MapType::iterator iter=eir[k].first; iter!=eir[k].second; ++iter)
        {
            iter->second=omp_get_thread_num();
        }
        
    }
    
    for (MapType::const_iterator iter=v.begin(); iter!=v.end(); ++iter)
    {
        std::cout<<iter->first<<" "<<iter->second<<std::endl;
    }
    
    return 0;
}

