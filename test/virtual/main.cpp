#include <iostream>
#include <string>
#include <sstream>      // std::istringstream
#include <vector>      // std::istringstream
#include <memory>      // std::istringstream

#include <model/IO/TextFileParser.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/DislocationDynamics/MobilityLaws/DislocationMobility.h>

using namespace model;


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


static std::unique_ptr<Base> selectDerived(const int& k)
{
    if(k==0)
    {
        return std::make_unique<Derived0>();
    }
    else if(k==1)
    {
        return std::make_unique<Derived1>();
    }
    else
    {
        std::cout<<"no implementation"<<std::endl;
        exit(1);
    }
    
    
}


static std::unique_ptr<DislocationMobilityBase> getMobility(const  Material<3,Isotropic>& material)
{
    if(material.crystalStructure=="FCC")
    {
        TextFileParser parser(material.materialFile);
        const double B1e=parser.readScalar<double>("B1e_SI",true);
        const double B1s=parser.readScalar<double>("B1s_SI",true);
        return std::make_unique<DislocationMobility<FCClattice<3>>>(material.b_SI,
                                                                      material.mu_SI,
                                                                      material.cs_SI,
                                                                      B1e,
                                                                      B1s);
    }
    else
    {
        std::cout<<"Unknown mobility for crystal structure '"<<material.crystalStructure<<"'. Exiting."<<std::endl;
        exit(EXIT_FAILURE);
    }
}


int main()
{
    std::unique_ptr<Base> p0=selectDerived(0);
    std::unique_ptr<Base> p1=selectDerived(1);
    
    std::cout<<p0->f()<<std::endl;
    std::cout<<p1->f()<<std::endl;
    
    Material<3,Isotropic> material("../../tutorials/DislocationDynamics/MaterialsLibrary/Cu.txt");
    
    TextFileParser parser(material.materialFile);
//    const double B1e=parser.readScalar<double>("B1e_SI",true);
//    const double B1s=parser.readScalar<double>("B1s_SI",true);
//    std::unique_ptr<DislocationMobilityBase> mobility=std::make_unique<DislocationMobility<FCClattice<3>>>(material.b_SI,
//                                                                                                           material.mu_SI,
//                                                                                                           material.cs_SI,
//                                                                                                           B1e,
//                                                                                                           B1s);
    
    std::unique_ptr<DislocationMobilityBase> mobility=getMobility(material);
    
    double v=mobility->velocity(Eigen::Matrix<double,3,3>(),
                                Eigen::Matrix<double,3,1>(),
                                Eigen::Matrix<double,3,1>(),
                                Eigen::Matrix<double,3,1>(),
                                300.0,
                                10.0,
                                1.0,
                                false);
    std::cout<<v<<std::endl;
    
    return 0;
}
