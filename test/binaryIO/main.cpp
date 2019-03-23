#include <iostream>
#include <fstream>
#include <vector>


/**********************************************************************/
template<typename T>
static void binWrite(std::ofstream& file,const T& o)
{
    file.write((char *) &o, sizeof (o));
}

/**********************************************************************/
template<typename T>
static void binRead(std::ifstream& file,
                    T*& memblock,
                    const size_t& arraySize)
{
    memblock = new T [arraySize];
    file.read (reinterpret_cast<char*>(memblock), arraySize*sizeof (T));
}

struct A
{
    
    int a;
    double b;
    
    A(){}
    A(int a1,double b1) :a(a1), b(b1){}
    
};

struct B
{
    size_t t;
    
    B(){}
    B(size_t t1) :t(t1){}
    
};


int main()
{
    // Create data
    std::vector<A> vA;
    vA.emplace_back(1,1.5);
    vA.emplace_back(2,2.5);
    vA.emplace_back(3,3.5);
    
    std::vector<B> vB;
    vB.emplace_back(1);
    vB.emplace_back(2);

    
    // Writing binary file
    const std::string filename("out.bin");
    std::ofstream outfile(filename.c_str(), std::ios::out  | std::ios::binary);
    
    const size_t sA(vA.size());
    binWrite(outfile,sA);

    const size_t sB(vB.size());
    binWrite(outfile,sB);

    
    for(const auto& a : vA)
    {
        binWrite(outfile,a);
    }
    
    for(const auto& b : vB)
    {
        binWrite(outfile,b);
    }
    
    outfile.close();

    
    // Reading binary file
    std::ifstream infile (filename.c_str(), std::ios::in|std::ios::binary);
    
    size_t* sizeA;
    binRead(infile,sizeA,1);
    std::cout<<"sizeA="<<sizeA[0]<<std::endl;

    size_t* sizeB;
    binRead(infile,sizeB,1);
    std::cout<<"sizeB="<<sizeB[0]<<std::endl;

    
        A* memblockA;
        binRead(infile,memblockA,sizeA[0]);


    
    B* memblockB;
    binRead(infile,memblockB,sizeB[0]);

    
    for(int i=0;i<sizeA[0];++i)
    {
        std::cout<<memblockA[i].a<<","<<memblockA[i].b<<std::endl;
    }
    
    for(int i=0;i<sizeB[0];++i)
    {
        std::cout<<memblockB[i].t<<std::endl;
    }
    
    return 0;
}
