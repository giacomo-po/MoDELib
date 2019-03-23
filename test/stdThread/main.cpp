
#include <iostream>
#include <iomanip>

#include <math.h>
#include <chrono>
//#include <ctime>
#include <thread>
#include <omp.h>
#include <Model/Threads/EqualIteratorRange.h>
//#include <Model/Threads/ParallelFor.h>



using namespace model;

int main (int argc, char * const argv[])
{
    
    int N=130000000;
    std::vector<int> v;
    v.reserve(N);
    for(int k=0;k<N;++k)
    {
        v.emplace_back(k);
    }
    
    
    //    auto t0 = std::chrono::system_clock::now();
    //    auto t0_sys = std::chrono::system_clock::to_time_t(t0);
    //    std::cout<<"Starting computation on "<<std::ctime(&t0_sys);
    
    auto t1= std::chrono::system_clock::now();
    //    double t1c=clock();
    
    const size_t nThreads = std::thread::hardware_concurrency();
    EqualIteratorRange<std::vector<int> > eir(v.begin(),v.end(),nThreads);
    
    
    
#pragma omp parallel for
    for (unsigned int k = 0; k < eir.size(); ++k)
    {
        for (std::vector<int>::const_iterator iter=eir[k].first; iter!=eir[k].second; ++iter)
        {
            double a = sin(*iter)+exp(cos(*iter))/exp(sin(*iter));
        }
        
    }
    
    std::cout<< "elapsed time= " <<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count() << "s\n";
    //    std::cout<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t1c)/CLOCKS_PER_SEC<<" sec]."<<std::endl;
    
    
    auto t2= std::chrono::system_clock::now();
    //    double t2c=clock();
    
#pragma omp parallel for
    for (unsigned int k = 0; k < v.size(); ++k)
    {
        double a = sin(v[k])+exp(cos(v[k]))/exp(sin(v[k]));
    }
    std::cout<< "elapsed time= " <<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count() << "s\n";
    //    std::cout<<std::setprecision(3)<<std::scientific<<" ["<<(clock()-t2c)/CLOCKS_PER_SEC<<" sec]."<<std::endl;
    
    
    
    return 0;
}




//#include <cmath>
//#include <omp.h>
//
//
//#include <thread>
//#include <vector>
//#include <assert.h>
//#include <iostream>


//void foo(long int load)
//{
////#pragma omp parallel for // ENABLE for openmp version
//    for (unsigned int i = 0; i < load; ++i)
//    {
//        double v = std::sin(i)+std::exp(std::cos(i))/std::exp(std::sin(i));
//    }
//}
//
//struct ParallelFor :
///* inheritance  */ private std::vector<std::thread>
//{
//
//    const long int nThreads;
//
//
//    ParallelFor(const long int& n = std::thread::hardware_concurrency()) :
//    /* init list */ nThreads(n)
//    {
//        std::cout<<"Creating ParallelFor with "<<nThreads<<" threads."<<std::endl;
//
//        assert(nThreads>0 && "NUMBER OF THREADS MUST BE >0");
//
//    }
//
//    template< class Function, class... Args >
//    void run(Function&& foo, Args&&... args )
//    {
//        for (auto i=0;i<nThreads;++i)
//        {
//            this->emplace_back(foo,args...);
//        }
//
//        for (auto i=0;i<this->size();++i)
//        {
//            this->operator[](i).join();
//        }
//
//        this->clear();
//
//    }
//
//};


//struct MyClass
//{
//
//    void bar(long int load)
//    {
//        //#pragma omp parallel for // ENABLE for openmp version
//        for (unsigned int i = 0; i < load; ++i)
//        {
//            double v = std::sin(i)+std::exp(std::cos(i))/std::exp(std::sin(i));
//        }
//    }
//
//
//
//    MyClass()
//    {
//        long int load = 800000000;
//
//        // std::thread version
//        ParallelFor pf;
//        pf.run(&MyClass::bar,this,load/pf.nThreads);
//
//
//    }
//
//};

