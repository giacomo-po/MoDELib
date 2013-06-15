

/* Compile
 g++ bench_LongestProcessingTime.cpp -o bench_LongestProcessingTime -O3 -std=c++11 -I../ -I/usr/local/include
*/

#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <model/MPI/LongestProcessingTime.h>
#include <time.h>       /* time */

using namespace model;

struct Dummy{};


int main (int argc, char * const argv[]) {
    
    srand (time(NULL));

    
    Dummy d;
    
    LongestProcessingTime<Dummy> lpt;

    // Insert may "light" Dummies with weight in [0,1]
    for(int k=0;k<1000;++k)
    {
        const double w(rand());
        lpt.insert(w/RAND_MAX,&d);
    }

    // Insert a few "heavy" Dummies with weight in [0,100]
    for(int k=0;k<100;++k)
    {
        const double w(rand());
        lpt.insert(w/RAND_MAX*10.0,&d);
    }
    
    int nP=(rand()*19.0)/RAND_MAX+1; // a random number of processors in [1 20]
    
    lpt.partition(nP);
    

    lpt.show();
    
    return 0;
}
