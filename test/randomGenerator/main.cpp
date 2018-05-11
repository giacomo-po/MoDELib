#include <iostream>
#include <random>
#include <chrono>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif


int main()
{
    typedef std::default_random_engine T;
    std::vector<T> generators;
    std::vector<std::normal_distribution<double>> distributions;

    const int seed=1;
    generators.resize(omp_get_max_threads());
    distributions.resize(omp_get_max_threads());
    for (int t=0;t<omp_get_max_threads();++t)
    {
        generators[t].seed(seed);
    }
    //generator(std::chrono::system_clock::now().time_since_epoch().count());
    //generator.seed(3);
    std::normal_distribution<double> distribution;


    for(int k=0;k<3;++k)
    {
//#pragma omp parallel for
        for (int t=0;t<omp_get_max_threads();++t)
        {
            std::cout<<distributions[t](generators[t])<<" ";
        }
        std::cout<<std::endl;
    }

    
    return 0;
}
