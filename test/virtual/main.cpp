#include <iostream>
#include <string>
#include <sstream>      // std::istringstream
#include <vector>      // std::istringstream
#include <memory>      // std::istringstream
#include <model/DislocationDynamics/Materials/BCClattice.h>      // std::istringstream
#include <model/DislocationDynamics/Materials/Material.h>      // std::istringstream
#include <model/DislocationDynamics/Materials/SingleCrystal.h>      // std::istringstream
#include <model/DislocationDynamics/Polycrystals/PolyCrystal.h>      // std::istringstream

#include <model/IO/TextFileParser.h>

//#include <istringstream>

struct Base
{

    virtual int f() const = 0;
    
};

struct Derived0 : public Base
{
    
    int f() const { return 0;}
    
};


struct Derived1 : public Base
{
    
    int f() const { return 1;}
    
};

using namespace model;

int main()
{
    
    std::vector<std::unique_ptr<Base>> v;
    v.emplace_back(new Derived0());
    v.emplace_back(new Derived1());
    
    std::cout<<v[0]->f()<<std::endl;
    std::cout<<v[1]->f()<<std::endl;
    
    BCClattice<3> bcc;
    std::cout<<bcc.latticeBasis<<std::endl;
    
    Lattice<3> lat=bcc;
    bcc.rotate( Eigen::AngleAxisd(0.25*M_PI, Eigen::Vector3d::UnitX()).toRotationMatrix());
    
    std::cout<<bcc.latticeBasis<<std::endl;

    std::cout<<lat.latticeBasis<<std::endl;

    TextFileParser parser("input.txt");
    std::string material=parser.readString("material");
    
//    long a;
//    a=parser.readScalar<decltype(a)>("a");
    
    SingleCrystal<3,Isotropic> cr(material);

    
//    cr.rotate( Eigen::AngleAxisd(0.0*M_PI, Eigen::Vector3d::UnitX()).toRotationMatrix());
//    std::cout<<cr.latticeBasis()<<std::endl;

//    std::cout<<cr2.latticeBasis()<<std::endl;
    
    return 0;
}
