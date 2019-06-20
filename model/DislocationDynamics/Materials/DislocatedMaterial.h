/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DISLOCATEDMATERIAL_H_
#define model_DISLOCATEDMATERIAL_H_


#include <cmath>
#include <string>
//#include <PeriodicElement.h>
//#include <CrystalOrientation.h>
#include <MPIcout.h> // defines mode::cout
#include <TerminalColors.h> // defines mode::cout
//#include <EigenDataReader.h> // defines mode::cout
#include <MaterialSymmetry.h>
//#include <MaterialBase.h>
//#include <BCClattice.h>
//#include <FCClattice.h>
#include <TextFileParser.h>
#include <MaterialBase.h>
#include <DislocationMobilityBase.h>
#include <DislocationMobilityBCC.h>
#include <DislocationMobilityFCC.h>
#include <DislocationMobilityHEXbasal.h>
#include <DislocationMobilityHEXprismatic.h>
#include <DislocationMobilityHEXpyramidal.h>



namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template<int dim,typename SymmetryType>
    class DislocatedMaterial { };
    
    /**************************************************************************/
    /**************************************************************************/
    template<int dim>
    class DislocatedMaterial<dim,Isotropic> : public MaterialBase
    {
        
        
        /**********************************************************************/
        static std::map<std::string,std::shared_ptr<DislocationMobilityBase>> getMobilities(const MaterialBase& materialBase)
        {
            
            std::map<std::string,std::shared_ptr<DislocationMobilityBase>> temp;
            
            if(materialBase.crystalStructure=="BCC")
            {
                temp.emplace("bcc",new DislocationMobilityBCC(materialBase));
            }
            else if(materialBase.crystalStructure=="FCC")
            {
                temp.emplace("fcc",new DislocationMobilityFCC(materialBase));
            }
            else if(materialBase.crystalStructure=="HEX")
            {
                temp.emplace("hexBasal",new DislocationMobilityHEXbasal(materialBase));
                temp.emplace("hexPrismatic",new DislocationMobilityHEXprismatic(materialBase));
                temp.emplace("hexPyramidal",new DislocationMobilityHEXpyramidal(materialBase));
            }
            else
            {
                std::cout<<"Unknown mobility for crystal structure '"<<materialBase.crystalStructure<<"'. Exiting."<<std::endl;
                exit(EXIT_FAILURE);
            }
            
            return temp;
        }
        
        
    public:
        
        
        
        // Universal constants

//        const std::string crystalStructure;
//
//        // Material constants in SI units
//        const double T;         // Absolute temparature [K]
//        const double Tm;        // Melting temparature [K]
//        const double mu0_SI;    // [Pa]
//        const double mu1_SI;    // [Pa/K]
//        const double mu_SI;     // temperature-dependent shear modulus mu=m0+m1*T [Pa]
//        const double nu;        // Poisson's ratio
//        const double rho_SI;    // mass density [Kg/m^3]
//        const double cs_SI;     // shear wave speed [m/s]
//        const double b_SI;      // Burgers vector [m]
//
//        // Material constants in code units
//        const double kB;        // Boltzmann constant [-]
//        const double mu;        // shear modulus [-]
//        const double b;         // Burgers vector [-]
//        const double cs;        // shear wave speed [-]
        const std::map<std::string,std::shared_ptr<DislocationMobilityBase>> dislocationMobilities;
        
        static double C1;        // 1-nu
        static double C2;        // 1.0/(4.0*M_PI*C1)
        static double C3;        // 1.0-2.0*nu;
        static double C4;        // 0.5*C2;
        
        
        /**********************************************************************/
        DislocatedMaterial(const std::string& fileName) :
        /* init */ MaterialBase(fileName)
//        /* init */,crystalStructure(TextFileParser(materialFile).readString("crystalStructure",true))
//        /* init */,T(TextFileParser(materialFile).readScalar<double>("T",true))
//        /* init */,Tm(TextFileParser(materialFile).readScalar<double>("Tm",true))
//        /* init */,mu0_SI(TextFileParser(materialFile).readScalar<double>("mu0_SI",true))
//        /* init */,mu1_SI(TextFileParser(materialFile).readScalar<double>("mu1_SI",true))
//        /* init */,mu_SI(mu0_SI+mu1_SI*T)
//        /* init */,nu(TextFileParser(materialFile).readScalar<double>("nu",true))
//        /* init */,rho_SI(TextFileParser(materialFile).readScalar<double>("rho_SI",true))
//        /* init */,cs_SI(sqrt(mu_SI/rho_SI))
//        /* init */,b_SI(TextFileParser(materialFile).readScalar<double>("b_SI",true))
//        /* init */,kB(kB_SI/mu_SI/std::pow(b_SI,3))
//        /* init */,mu(1.0)
//        /* init */,b(1.0)
//        /* init */,cs(1.0)
        /* init */,dislocationMobilities(getMobilities(*this))
        {
            model::cout<<magentaColor<<"  units of stress (shear modulus): mu="<<mu_SI<<" [Pa]"<<std::endl;
            model::cout<<magentaColor<<"  units of length (Burgers vector): b="<<b_SI<<" [m]"<<std::endl;
            model::cout<<magentaColor<<"  units of speed (shear-wave speed): cs="<<cs_SI<<" [m/s]"<<std::endl;
            model::cout<<magentaColor<<"  units of time: b/cs="<<b_SI/cs_SI<<" [sec]"<<defaultColor<<std::endl;

            C1=1.0-nu;
            C2=1.0/(4.0*M_PI*C1);
            C3=1.0-2.0*nu;
            C4=0.5*C2;

            
        }
        
//        /**********************************************************************/
//        static double velocity(const Eigen::Matrix<double,3,3>& S,
//                               const Eigen::Matrix<double,3,1>& b,
//                               const Eigen::Matrix<double,3,1>& xi,
//                               const Eigen::Matrix<double,3,1>& n,
//                               const double& dL,
//                               const double& dt,
//                               const bool& use_stochasticForce)
//        {
//            
//            assert(0 && "FINISH VELOCITY");
//        }
        
    };
    
    template<int dim>
    double DislocatedMaterial<dim,Isotropic>::C1=1.0-0.34;    // 1-nu
    
    template<int dim>
    double DislocatedMaterial<dim,Isotropic>::C2=1.0/(4.0*M_PI*(1.0-0.34));  // 1/(4*pi*(1-nu))

    template<int dim>
    double DislocatedMaterial<dim,Isotropic>::C3=1.0-2.0*0.34; // 1-2*nu
    
    template<int dim>
    double DislocatedMaterial<dim,Isotropic>::C4=1.0/(8.0*M_PI*(1.0-0.34));  // 1/(4*pi*(1-nu))
    
}
#endif

//        /**********************************************************************/
//        template<int Z>
//        static void select()
//        {
//            typedef PeriodicElement<Z,Isotropic> IM;
//
//            mu=1.0;     // mu is used to normalized stress
//            b =1.0;     // b is used to normalized length
//            cs=1.0;     // sc*=sqrt(mu*/rho*) is used to normalize speed
//
//            b_real=PeriodicElement<Z,Isotropic>::b;
//            //            B =1.0;     // B=A*T [Pa*sec] is effectively used to normalize time
//            //            Binv=1.0/B;
//            //            rho=IM::mu*std::pow(IM::b/(IM::Ae*T),2)*IM::rho;  //! rho* = mu*(b/B)^2 * rho, B=A*T
//            //            //        cs=IM::B/IM::b*std::pow(IM::rho*IM::mu,-0.5);
//            //            cs=sqrt(mu/rho); // sc*=sqrt(mu*/rho*)
//            ////            cs =1.0;    // shear wave speed is used to normalized velocity
//            nu=IM::nu;
//            lambda=2.0*mu*nu/(1.0-2.0*nu);
//            C1=1.0-nu;
//            C2=1.0/(4.0*M_PI*C1);
//            C3=1.0-2.0*nu;
//            C4=0.5*C2;
//
//            /* still unused, to be used in BCC mobility */
//            kB=1.38e-23/IM::mu/std::pow(IM::b,3); // [-]
//
//            model::cout<<magentaColor<<"Material is now: "<<IM::name<<defaultColor<<std::endl;
//            model::cout<<greenColor<<"  units of stress (shear modulus): mu="<<IM::mu<<" [Pa]"<<std::endl;
//            model::cout<<greenColor<<"  units of length (Burgers vector): b="<<IM::b<<" [m]"<<std::endl;
//            model::cout<<greenColor<<"  units of speed (shear-wave speed): cs="<<IM::cs<<" [m/s]"<<std::endl;
//            model::cout<<greenColor<<"  units of time: b/cs="<<IM::b/IM::cs<<" [sec]"<<defaultColor<<std::endl;
//        }


//    public:

//        Material(const std::string& materialFileName)
//        {
//
////            EigenDataReader EDR;
//
//
//
//
//        }

//        static double b;
//        static double b_real;
//        static double mu;
//        static double cs;
//
//        static double T;
//        //        static double rho;
//        static double nu;
//        static double lambda;
//        static double C1;
//        static double C2;
//        static double C3;
//        static double C4;
//        static double kB;



//        static Eigen::Matrix<double,Eigen::Dynamic,2> dH0;

//        /**********************************************************************/
//        static void select(const unsigned int& Z)
//        {/*!\param[in] Z the atomic number of the element to be selected
//          *
//          * Selects the element Z as the current material. All material properties
//          * are updated accordingly.
//          */
//            switch (Z)
//            {
//                case Al.Z:
//                    selectedMaterial=Al.Z;
//                    select<Al.Z>();
//                    break;
//                case Ni.Z:
//                    selectedMaterial=Ni.Z;
//                    select<Ni.Z>();
//                    break;
//                case Cu.Z:
//                    selectedMaterial=Cu.Z;
//                    select<Cu.Z>();
//                    break;
//                    //                case Fe:
//                    //                    selectedMaterial=Fe;
//                    //                    select<Fe>();
//                    //                    break;
//                case W.Z:
//                    selectedMaterial=W.Z;
//                    select<W.Z>();
//                    break;
//                default:
//                    assert(0 && "Material not implemented.");
//                    break;
//            }
//
//        }

//        /**********************************************************************/
//        static double velocity(const Eigen::Matrix<double,3,3>& S,
//                               const Eigen::Matrix<double,3,1>& b,
//                               const Eigen::Matrix<double,3,1>& xi,
//                               const Eigen::Matrix<double,3,1>& n,
//                               const double& dL,
//                               const double& dt,
//                               const bool& use_stochasticForce)
//        {
//            switch (selectedMaterial)
//            {
//                case Al.Z:
//                    return Al.dm.velocity(S,b,xi,n,T,dL,dt,use_stochasticForce);
//                case Ni.Z:
//                    return Ni.dm.velocity(S,b,xi,n,T,dL,dt,use_stochasticForce);
//                case Cu.Z:
//                    return Cu.dm.velocity(S,b,xi,n,T,dL,dt,use_stochasticForce);
//                    //              case Fe:
//                    //                  return PeriodicElement<Cu,Isotropic>::dm.velocity(S,b,xi,n,T);
//                case W.Z:
//                    return W.dm.velocity(S,b,xi,n,T,dL,dt,use_stochasticForce);
//
//                default:
//                    assert(0 && "velocity function not implemented.");
//                    return 0.0;
//            }
//        }

//    template<int dim,typename SymmetryType>
//    class Crystal : public Lattice<dim>
////                    public Material<dim,SymmetryType>
//
//    {
//
//
//        static Lattice<dim> getLattice(const std::string& inputFileName)
//        {
//
//            std::string crystalStructure="FCC";
//            if(crystalStructure=="BCC")
//            {
//                return BCClattice<dim>();
//            }
//            else if(crystalStructure=="FCC")
//            {
//                return FCClattice<dim>();
//            }
//            else
//            {
//                exit(EXIT_FAILURE);
//            }
//
//
//        }
//
//    public:
//
//        Crystal(const std::string& inputFileName) :
//        Lattice<dim>(getLattice(inputFileName))
//        {
//
//        }
//
//
//    };

//    // Static data
//    int Material<Isotropic>::selectedMaterial=29;
//    double Material<Isotropic>::b=1.0;    // dimensionless Burgers vector
//    double Material<Isotropic>::b_real=0.2556e-9;    // dimensionless Burgers vector
//    double Material<Isotropic>::mu=1.0;
//    double Material<Isotropic>::cs=1.0;
//    double Material<Isotropic>::T=300.0;  // Temperature [K]
//    double Material<Isotropic>::nu=0.34;
//    double Material<Isotropic>::lambda=2.0*1.0*0.34/(1.0-2.0*0.34);
//    double Material<Isotropic>::C1=1.0-0.34;    // 1-nu
//    double Material<Isotropic>::C2=1.0/(4.0*M_PI*(1.0-0.34));  // 1/(4*pi*(1-nu))
//    double Material<Isotropic>::C3=1.0-2.0*0.34; // 1-2*nu
//    double Material<Isotropic>::C4=1.0/(8.0*M_PI*(1.0-0.34));  // 1/(4*pi*(1-nu))
//    double Material<Isotropic>::kB=1.38e-23/48.0e9/std::pow(0.2556e-9,3);  // Boltzmann constant
