/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpF@ucla.edu>.
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneNoise_H
#define model_GlidePlaneNoise_H

#include <cmath>
#include <random>

#include <Eigen/Dense>


#ifdef _MODEL_GLIDE_PLANE_NOISE_GENERATOR_
//#include <complex.h>
#include <fftw3.h>
#include <boost/math/special_functions/bessel.hpp>
#endif

//#include <DDtraitsIO.h>
#include <PolycrystallineMaterialBase.h>
#include <TerminalColors.h>
//#include <SlipSystem.h>
//#include <GlidePlane.h>

namespace model
{

    struct NoiseTraitsBase
    {
        typedef double REAL_SCALAR;
        typedef std::complex<double> COMPLEX;
        typedef Eigen::Array<int,2,1> GridSizeType;
        typedef Eigen::Array<double,2,1> GridSpacingType;
    };

    template <int N>
    struct NoiseTraits
    {
        typedef typename NoiseTraitsBase::REAL_SCALAR REAL_SCALAR;
        typedef typename NoiseTraitsBase::COMPLEX COMPLEX;
        typedef typename NoiseTraitsBase::GridSizeType GridSizeType;
        typedef Eigen::Matrix<REAL_SCALAR,1,N> NoiseType;
        typedef std::vector<NoiseType> NoiseContainerType;
    };

    template<>
    struct NoiseTraits<1>
    {
        typedef typename NoiseTraitsBase::REAL_SCALAR REAL_SCALAR;
        typedef typename NoiseTraitsBase::COMPLEX COMPLEX;
        typedef typename NoiseTraitsBase::GridSizeType GridSizeType;
        typedef REAL_SCALAR NoiseType;
        typedef std::vector<NoiseType> NoiseContainerType;
    };

struct SolidSolutionNoiseReader : public NoiseTraits<2>::NoiseContainerType
{
    
    typedef typename NoiseTraitsBase::REAL_SCALAR REAL_SCALAR;
    typedef typename NoiseTraitsBase::COMPLEX COMPLEX;
    typedef typename NoiseTraitsBase::GridSizeType GridSizeType;
    typedef typename NoiseTraitsBase::GridSpacingType GridSpacingType;
    typedef typename NoiseTraits<2>::NoiseType NoiseType;
    typedef typename NoiseTraits<2>::NoiseContainerType NoiseContainerType;
    
    
    
    static int LittleEndian()
    {
        int num = 1;
        if(*(char *)&num == 1)
        {
            return 1;       //little endian
        }
        else
        {
            return 0;       // big endian
        }
    }
    
    static float ReverseFloat( const float inFloat )
    {
        float retVal;
        char *FloatToConvert = ( char* ) & inFloat;
        char *returnFloat = ( char* ) & retVal;
        
        // swap the bytes into a temporary buffer
        returnFloat[0] = FloatToConvert[3];
        returnFloat[1] = FloatToConvert[2];
        returnFloat[2] = FloatToConvert[1];
        returnFloat[3] = FloatToConvert[0];
        
        return retVal;
    }
    
    static double ReverseDouble( const double inDouble )
    {
        double retVal;
        char *DoubleToConvert = ( char* ) & inDouble;
        char *returnDouble = ( char* ) & retVal;
        
        // swap the bytes into a temporary buffer
        returnDouble[0] = DoubleToConvert[7];
        returnDouble[1] = DoubleToConvert[6];
        returnDouble[2] = DoubleToConvert[5];
        returnDouble[3] = DoubleToConvert[4];
        returnDouble[4] = DoubleToConvert[3];
        returnDouble[5] = DoubleToConvert[2];
        returnDouble[6] = DoubleToConvert[1];
        returnDouble[7] = DoubleToConvert[0];
        
        return retVal;
    }
    
    static std::pair<GridSizeType,GridSpacingType> Read_dimensions(const char *fname)
    {
        int NX, NY, NZ;
        double DX, DY, DZ;
        char line[200];
//        double temp;
//        int Nr_in;
//        char *dum;
//        int flag;
        FILE *InFile=fopen(fname,"r");
        
        if (InFile == NULL)
        {
            fprintf(stderr, "Can't open noise file %s\n",fname);
            exit(1);
        }
        
        for(int i=0;i<5;i++)
        {
//            char *dum=
            fgets(line, 200, InFile);
        }
//        int flag=
        fscanf(InFile, "%s %lf %lf %lf\n", line, &(DX), &(DY), &(DZ));
//        dum=fgets(line, 200, InFile);
//        flag=
        fscanf(InFile, "%s %d %d %d\n", line, &(NX), &(NY), &(NZ));
//        return std::make_pair((GridSizeType()<<NX,NY,NZ).finished(),(GridSpacingType()<<DX,DY,DZ).finished());
        return std::make_pair((GridSizeType()<<NX,NY).finished(),(GridSpacingType()<<DX,DY).finished());

    }
    
    
    
    static void Read_noise_vtk(const char *fname, REAL_SCALAR *Noise, int Nr,const double& MSS)
    {
        char line[200];
        double temp;
//        char *dum;
//        int flag;
        FILE *InFile=fopen(fname,"r");
        
        if (InFile == NULL)
        {
            fprintf(stderr, "Can't open noise file %s\n",fname);
            exit(1);
        }
        
        for(int i=0;i<10;i++)
        {
//            char *dum=
            fgets(line, 200, InFile);
        }
        
        if(LittleEndian()) // if machine works with LittleEndian
        {
            for(int ind=0;ind<Nr;ind++)
            {
//                int flag =
                fread(&temp, sizeof(double), 1, InFile);
                Noise[ind] = REAL_SCALAR(ReverseDouble(temp));
            }
        }
        else // if machine works with BigEndian
        {
            for(int ind=0;ind<Nr;ind++)
            {
//                int flag =
                fread(&temp, sizeof(double), 1, InFile);
                Noise[ind] = MSS*REAL_SCALAR(temp);
            }
        }
    }
    
    SolidSolutionNoiseReader(const std::string& noiseFile,const PolycrystallineMaterialBase& mat,
                             const GridSizeType& _gridSize, const GridSpacingType& _gridSpacing_A)
    {
        std::cout<<"Reading SolidSolutionNoise files"<<std::endl;
//        NoiseContainerType Noise;
        
        const std::string fileName_xz(TextFileParser(noiseFile).readString("solidSolutionNoiseFile_xz",true));
        const auto gridSize_xz(Read_dimensions(fileName_xz.c_str()));
        const std::string fileName_yz(TextFileParser(noiseFile).readString("solidSolutionNoiseFile_yz",true));
        const auto gridSize_yz(Read_dimensions(fileName_yz.c_str()));
        const double MSSS_SI(TextFileParser(mat.materialFile).readScalar<double>("MSSS_SI",true));
        const double MSS(std::sqrt(MSSS_SI)/mat.mu);

        if((gridSize_xz.first-gridSize_yz.first).matrix().squaredNorm()==0 && (gridSize_xz.first-_gridSize).matrix().squaredNorm()==0)
        {
            if((gridSize_xz.second-gridSize_yz.second).matrix().squaredNorm()==0.0 && (gridSize_xz.second-_gridSpacing_A).matrix().squaredNorm()==0.0)
            {
                const size_t Nr(gridSize_xz.first.array().prod());
                // allocate noises
                REAL_SCALAR *Noise_xz = (REAL_SCALAR *) malloc(sizeof(REAL_SCALAR)*Nr);
                // read noise vtk file
                Read_noise_vtk(fileName_xz.c_str(), Noise_xz, Nr,MSS);
                
                // allocate noises
                REAL_SCALAR *Noise_yz = (REAL_SCALAR *) malloc(sizeof(REAL_SCALAR)*Nr);
                // read noise vtk file
                Read_noise_vtk(fileName_yz.c_str(), Noise_yz, Nr,MSS);
                
                this->resize(Nr);
                for(size_t k=0;k<Nr;++k)
                {
                    this->operator[](k)<<Noise_xz[k],Noise_yz[k];
                }
            }
            else
            {
                std::cout<<"gridSpacing in "<<fileName_xz<<std::endl;
                std::cout<<gridSize_xz.second<<std::endl;
                std::cout<<"gridSpacing in "<<fileName_yz<<std::endl;
                std::cout<<gridSize_yz.second<<std::endl;
                std::cout<<"input gridSpacing_A = "<<_gridSpacing_A<<std::endl;
                throw std::runtime_error("gridSpacing mismatch.");
            }
        }
        else
        {
            std::cout<<"gridSize in "<<fileName_xz<<std::endl;
            std::cout<<gridSize_xz.first<<std::endl;
            std::cout<<"gridSize in "<<fileName_yz<<std::endl;
            std::cout<<gridSize_yz.first<<std::endl;
            std::cout<<"input gridSize = "<<_gridSize<<std::endl;
            throw std::runtime_error("gridSize mismatch.");
        }
    }
    
    
};

    class SolidSolutionNoise : public NoiseTraits<2>::NoiseContainerType
    {
        typedef typename NoiseTraitsBase::REAL_SCALAR REAL_SCALAR;
        typedef typename NoiseTraitsBase::COMPLEX COMPLEX;
        typedef typename NoiseTraitsBase::GridSizeType GridSizeType;
        typedef typename NoiseTraitsBase::GridSpacingType GridSpacingType;
        typedef typename NoiseTraits<2>::NoiseType NoiseType;
        typedef typename NoiseTraits<2>::NoiseContainerType NoiseContainerType;
        
        
        const NoiseContainerType& noiseVector() const
        {
            return *this;
        }
        
        NoiseContainerType& noiseVector()
        {
            return *this;
        }
        
        const GridSizeType gridSize;
        const GridSpacingType gridSpacing_A;

    public:
        
        SolidSolutionNoise(const std::string& noiseFile,const PolycrystallineMaterialBase& mat,
                           const GridSizeType& _gridSize, const GridSpacingType& _gridSpacing_A, const int& solidSolutionNoiseMode) :
        /* init */ gridSize(_gridSize)
        /* init */,gridSpacing_A(_gridSpacing_A)
        {
            
            switch (solidSolutionNoiseMode)
            {
                case 1:
                {// read noise
                    std::cout<<greenBoldColor<<"Reading SolidSolutionNoise"<<defaultColor<<std::endl;
                    noiseVector()=(SolidSolutionNoiseReader(noiseFile,mat,gridSize,gridSpacing_A));
                    break;
                }
                    
                case 2:
                {// compute noise
                    std::cout<<greenBoldColor<<"Computing SolidSolutionNoise"<<defaultColor<<std::endl;
                    break;
                }
                    
                default:
                    break;
            }
            
            
            NoiseType ave(NoiseType::Zero());
            for(const auto& valArr: noiseVector())
            {
                ave+=valArr;
            }
            ave/=noiseVector().size();
            
            NoiseType var(NoiseType::Zero());
            for(const auto& valArr: noiseVector())
            {
                var+= ((valArr-ave).array()*(valArr-ave).array()).matrix();
            }
            var/=noiseVector().size();
            
            std::cout<<"gridSize= "<<gridSize<<std::endl;
            std::cout<<"gridSpacing_A= "<<gridSpacing_A<<std::endl;
            std::cout<<"noiseAverage="<<ave<<std::endl;
            std::cout<<"noiseVariance="<<var<<std::endl;
        }
        
    };

//    struct SolidSolutionNoiseGenerator
//    {
//
//        //    int mod(int a, int b)
//        //    {
//        //        int r = a % b;
//        //        return r < 0 ? r + b : r;
//        //    }
//
//        typedef typename NoiseTraits<2>::REAL_SCALAR REAL_SCALAR;
//        typedef typename NoiseTraits<2>::COMPLEX COMPLEX;
//        typedef typename NoiseTraits<2>::GridSizeType GridSizeType;
//        typedef typename NoiseTraits<2>::NoiseType NoiseType;
//        typedef typename NoiseTraits<2>::NoiseContainerType NoiseContainerType;
//
//        int NX, NY, NZ;
//        REAL_SCALAR DX, DY, DZ;
//        REAL_SCALAR a;
//        REAL_SCALAR a_cai;
//        int seed;
//        //    int flag;
//        REAL_SCALAR LX, LY, LZ;
//        REAL_SCALAR DV;
//        int NR;
//        int NK;
//        REAL_SCALAR Norm;
//
//        SolidSolutionNoiseGenerator() :
//        /*init*/ NX(256)     // dimension along x
//        /*init*/,NY(256)     // dimension along y
//        /*init*/,NZ(64)      // dimension along z
//        /*init*/,DX(1.0)     // grid spacing [AA]
//        /*init*/,DY(1.0)     // grid spacing [AA]
//        /*init*/,DZ(1.0)     // grid spacing [AA]
//        /*init*/,a(1.0)      // spreading length for stresses [AA]
//        /*init*/,a_cai(3.0)  // spreading length for non-singular dislocaion theory [AA]
//        /*init*/,seed(1234)  // random seed
//        /*init*/,LX(NX*DX)
//        /*init*/,LY(NY*DY)
//        /*init*/,LZ(NZ*DZ)
//        /*init*/,DV(DX*DY*DZ)
//        /*init*/,NR(NX*NY*NZ)
//        /*init*/,NK(NX*NY*(NZ/2+1))
//        /*init*/,Norm(1./REAL_SCALAR(NR))
//        {
//
//        }
//
//    #ifdef _MODEL_GLIDE_PLANE_NOISE_GENERATOR_
//
//        // Cai doubly-convoluted spreading function in Fourier space
//        REAL_SCALAR Wk_Cai(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz, REAL_SCALAR a) const
//        {
//            REAL_SCALAR k = sqrt(kx*kx + ky*ky + kz*kz);
//            if(k>0)
//            {
//                return a*k*sqrt(0.5*boost::math::cyl_bessel_k(2,a*k));
//            }
//            else
//            {
//                return 1.;
//            }
//
//        }
//
//
//        // Cai spreading function
//        REAL_SCALAR W_Cai(REAL_SCALAR r2, REAL_SCALAR a) const
//        {
//            return 15.*a*a*a*a/(8.*M_PI*pow(r2+a*a,7./2.));
//        }
//
//        REAL_SCALAR W_t_Cai(REAL_SCALAR r2, REAL_SCALAR a) const
//        {
//            return 0.3425*W_Cai(r2,0.9038*a) + 0.6575*W_Cai(r2,0.5451*a);
//        }
//
//        // normalized auto-correlation function in Fourier space for sigma_xy
//        REAL_SCALAR S_xy_k(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz) const
//        {
//            REAL_SCALAR k2 = kx*kx + ky*ky + kz*kz;
//            return 120.*M_PI*sqrt(M_PI)*a*a*a/(LX*LY*LZ)*(kx*kx*ky*ky)/(k2*k2)*exp(-a*a*k2);
//        }
//
//        // normalized auto-correlation function in Fourier space for sigma_xz
//        REAL_SCALAR S_xz_k(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz) const
//        {
//            REAL_SCALAR k2 = kx*kx + ky*ky + kz*kz;
//            return 120.*M_PI*sqrt(M_PI)*a*a*a/(LX*LY*LZ)*(kx*kx*kz*kz)/(k2*k2)*exp(-a*a*k2);
//        }
//
//        // normalized auto-correlation function in Fourier space for sigma_yz
//        REAL_SCALAR S_yz_k(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz) const
//        {
//            REAL_SCALAR k2 = kx*kx + ky*ky + kz*kz;
//            return 120.*M_PI*sqrt(M_PI)*a*a*a/(LX*LY*LZ)*(ky*ky*kz*kz)/(k2*k2)*exp(-a*a*k2);
//        }
//
//        SolidSolutionNoise getNoiseBase() const
//        {
//            GridSizeType gridSize;
//            NoiseContainerType noise;
//
//            std::cout<<"Computing SolidSolutionNoise"<<std::endl;
//            int ind = 0;
//            REAL_SCALAR kx,ky,kz;
//
//            // read input command line
//
//            //        char* fname_out_yz = argv[10];  // name of the output vtk
//            //        char* fname_out_xz = argv[11];  // name of the output vtk
//
//            REAL_SCALAR var, std;
//
//
//
//            REAL_SCALAR *Rr_yz;          // pointer of the REAL_SCALAR space data
//            COMPLEX *Rk_yz;       // pointer of the fourier space data
//
//            REAL_SCALAR *Rr_xz;          // pointer of the REAL_SCALAR space data
//            COMPLEX *Rk_xz;       // pointer of the fourier space data
//
//            fftw_plan plan_R_yz_r2c, plan_R_yz_c2r;  // fft plan
//            fftw_plan plan_R_xz_r2c, plan_R_xz_c2r;  // fft plan
//
//            //    REAL_SCALAR *Sig_xz;
//            //    REAL_SCALAR *Sig_yz;
//
//            // allocate
//            Rr_yz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR);
//            Rk_yz = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*NK);
//            Rr_xz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR);
//            Rk_xz = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*NK);
//
//
//            //    Sig_xy = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR);
//            //    Sig_xz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR);
//            //    Sig_yz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR);
//
//            // prepare plans
//            //    plan_W_r2c = fftw_plan_dft_r2c_3d(NX, NY, NZ, Wr, reinterpret_cast<fftw_complex*>(Wk), FFTW_ESTIMATE);
//            //    plan_W_c2r = fftw_plan_dft_c2r_3d(NX, NY, NZ, reinterpret_cast<fftw_complex*>(Wk), Wr, FFTW_ESTIMATE);
//            plan_R_yz_r2c = fftw_plan_dft_r2c_3d(NX, NY, NZ, Rr_yz, reinterpret_cast<fftw_complex*>(Rk_yz), FFTW_ESTIMATE);
//            plan_R_yz_c2r = fftw_plan_dft_c2r_3d(NX, NY, NZ, reinterpret_cast<fftw_complex*>(Rk_yz), Rr_yz, FFTW_ESTIMATE);
//            plan_R_xz_r2c = fftw_plan_dft_r2c_3d(NX, NY, NZ, Rr_xz, reinterpret_cast<fftw_complex*>(Rk_xz), FFTW_ESTIMATE);
//            plan_R_xz_c2r = fftw_plan_dft_c2r_3d(NX, NY, NZ, reinterpret_cast<fftw_complex*>(Rk_xz), Rr_xz, FFTW_ESTIMATE);
//
//
//
//            // define the noise in Fourier space based on the auto-correlation function [Geslin et al. JPMS 2021]
//            REAL_SCALAR Nk_yz, Mk_yz;
//            REAL_SCALAR Nk_xz, Mk_xz;
//
//            std::default_random_engine generator(seed);
//            std::normal_distribution<REAL_SCALAR> distribution(0.0,1.0);
//
//
//            /////////////////////////////////////////////////////
//            // generate yz and xz correlated component //
//            /////////////////////////////////////////////////////
//            for(int i=0; i<NX; i++)
//            {
//                for(int j=0; j<NY; j++)
//                {
//                    for(int k=0; k<(NZ/2+1); k++)
//                    {
//                        ind = NY*(NZ/2+1)*i + j*(NZ/2+1) + k;
//
//                        kx = 2.*M_PI/LX*REAL_SCALAR(i);
//                        if(i>NX/2)
//                        {
//                            kx = 2.*M_PI/LX*REAL_SCALAR(i-NX);
//                        }
//
//                        ky = 2*M_PI/LY*REAL_SCALAR(j);
//                        if(j>NY/2)
//                        {
//                            ky = 2.*M_PI/LY*REAL_SCALAR(j-NY);
//                        }
//
//                        kz = 2.*M_PI/LZ*REAL_SCALAR(k);
//
//                        // random numbers
//                        Nk_yz = distribution(generator);
//                        Mk_yz = distribution(generator);
//                        if(kx*ky>=0)
//                        {
//                            Nk_xz = Nk_yz;
//                            Mk_xz = Mk_yz;
//                        }
//                        else
//                        {
//                            Nk_xz = -Nk_yz;
//                            Mk_xz = -Mk_yz;
//                        }
//
//                        //                if(strcmp(comp,"corr")==0)
//                        //                {
//                        if(k==0) // /!\ special case for k=0 and k==NZ/2 because of folding of C2R Fourier transform
//                        {
//                            Rk_yz[ind] = sqrt(S_yz_k(kx,ky,kz))*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
//                            Rk_xz[ind] = sqrt(S_xz_k(kx,ky,kz))*(Nk_xz+Mk_xz*COMPLEX(0.0,1.0));
//                        }
//                        else if(k==NZ/2)
//                        {
//                            Rk_yz[ind] = sqrt(S_yz_k(kx,ky,kz))*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
//                            Rk_xz[ind] = sqrt(S_xz_k(kx,ky,kz))*(Nk_xz+Mk_xz*COMPLEX(0.0,1.0));
//                        }
//                        else
//                        {
//                            Rk_yz[ind] = sqrt(S_yz_k(kx,ky,kz)/2.)*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
//                            Rk_xz[ind] = sqrt(S_xz_k(kx,ky,kz)/2.)*(Nk_xz+Mk_xz*COMPLEX(0.0,1.0));
//                        }
//
//                        if(a_cai>0)
//                        {
//                            Rk_yz[ind] = Rk_yz[ind]*Wk_Cai(kx, ky, kz, a_cai);
//                            Rk_xz[ind] = Rk_xz[ind]*Wk_Cai(kx, ky, kz, a_cai);
//                        }
//                        //                }
//                    }
//                }
//            }
//            Rk_yz[0] = 0;
//            Rk_xz[0] = 0;
//
//
//            // FFT back to REAL_SCALAR space
//            fftw_execute(plan_R_yz_c2r);
//            fftw_execute(plan_R_xz_c2r);
//
//            gridSize[0]=NX;
//            gridSize[1]=NY;
//            gridSize[2]=NZ;
//
//            noise.reserve(NR);
//            for(int i=0;i<NX;i++)
//            {
//                for(int j=0;j<NY;j++)
//                {
//                    for(int k=0;k<NZ;k++)
//                    {
//                        ind = NY*NZ*i + j*NZ + k;
//                        //                    noise.push_back(std::array<double,2>{SolidSolutionNoise::ReverseDouble(double(Rr_xz[ind])),SolidSolutionNoise::ReverseDouble(double(Rr_yz[ind]))});
//                        noise.push_back((NoiseType()<<Rr_xz[ind],Rr_yz[ind]).finished());
//
//                        //                    std::cout<<noise.back()[0]<<" "<<noise.back()[1]<<std::endl;
//                        //                    temp=ReverseDouble(double(F[ind]));
//                        //                    fwrite(&temp, sizeof(double), 1, OutFile);
//                    }
//                }
//            }
//            //        for(size_t k=0;k<NR;++k)
//            //        {
//            //            noise.push_back(std::array<double,2>{Rr_xz[k],Rr_yz[k]});
//            //        }
//
//            return SolidSolutionNoise(gridSize,noise);
//        }
//    #else
//
//        SolidSolutionNoise getNoiseBase() const
//        {
//            GridSizeType gridSize;
//            NoiseContainerType  noise;
//            return SolidSolutionNoise(gridSize,noise);
//        }
//    #endif
//
//
//
//
//    };




    struct StackingFaultNoise : public NoiseTraits<1>::NoiseContainerType
    {
        typedef typename NoiseTraits<1>::REAL_SCALAR REAL_SCALAR;
//        typedef typename NoiseTraits<2>::COMPLEX COMPLEX;
        typedef typename NoiseTraits<1>::GridSizeType GridSizeType;
        typedef typename NoiseTraits<1>::NoiseType NoiseType;
        typedef typename NoiseTraits<1>::NoiseContainerType NoiseContainerType;
        
        std::default_random_engine generator;


        
        StackingFaultNoise(const std::string&, // noiseFile
                           const PolycrystallineMaterialBase& mat,
                           const NoiseTraitsBase::GridSizeType& gridSize,
                           const NoiseTraitsBase::GridSpacingType& gridSpacing_SI)
        {
            
            std::cout<<greenBoldColor<<"Creating StackingFaultNoise"<<defaultColor<<std::endl;
            
//            const std::string noiseFileName(traitsIO.inputFilesFolder+"/"+TextFileParser(traitsIO.noiseFile).readString("stackingFaultNoiseFile"));
            const double isfEnergyDensityMEAN(TextFileParser(mat.materialFile).readScalar<double>("isfEnergyDensityMEAN_SI",true)/(mat.mu_SI*mat.b_SI));
            const double isfEnergyDensitySTD(TextFileParser(mat.materialFile).readScalar<double>("isfEnergyDensitySTD_SI",true)/std::sqrt(gridSpacing_SI(0)*gridSpacing_SI(1))/(mat.mu_SI*mat.b_SI));
                        
            std::normal_distribution<double> distribution (isfEnergyDensityMEAN,isfEnergyDensitySTD);

            const size_t N(gridSize.array().prod());
            this->reserve(N);
            for(size_t k=0;k<N;++k)
            {// J/m^2 = N/m = Pa*m
//                this->push_back((distribution(generator)-isfEnergyDensityMEAN));
                this->push_back(distribution(generator));
            }
            
            NoiseType ave(0.0);
            for(const auto& valArr: *this)
            {
                ave+=valArr;
            }
            ave/=this->size();
            
            NoiseType var(0.0);
            for(const auto& valArr: *this)
            {
                var+= (valArr-ave)*(valArr-ave);
            }
            var/=this->size();
            
            std::cout<<"gridSize= "<<gridSize<<std::endl;
            std::cout<<"gridSpacing_SI= "<<gridSpacing_SI<<std::endl;
            std::cout<<"noiseAverage="<<ave<<std::endl;
            std::cout<<"noiseVariance="<<var<<std::endl;
            
        }
        
    };

    struct GlidePlaneNoise
    {
        typedef typename Eigen::Matrix<double,3,1> VectorDim;

        const NoiseTraitsBase::GridSizeType gridSize;
        const NoiseTraitsBase::GridSpacingType gridSpacing_SI;
        const NoiseTraitsBase::GridSpacingType gridSpacing;
        const int solidSolutionNoiseMode;
        const int stackingFaultNoiseMode;

        const std::shared_ptr<SolidSolutionNoise> solidSolution;
        const std::shared_ptr<StackingFaultNoise> stackingFault;
                
        
        GlidePlaneNoise(const std::string& noiseFile,const PolycrystallineMaterialBase& mat) :
        /* init */ gridSize(TextFileParser(noiseFile).readMatrix<int,1,2>("gridSize",true))
        /* init */,gridSpacing_SI(TextFileParser(noiseFile).readMatrix<double,1,2>("gridSpacing_SI",true))
        /* init */,gridSpacing(gridSpacing_SI/mat.b_SI)
        /* init */,solidSolutionNoiseMode(TextFileParser(noiseFile).readScalar<int>("solidSolutionNoiseMode"))
        /* init */,stackingFaultNoiseMode(TextFileParser(noiseFile).readScalar<int>("stackingFaultNoiseMode"))
        /* init */,solidSolution(solidSolutionNoiseMode? new SolidSolutionNoise(noiseFile,mat,gridSize,gridSpacing_SI*1.0e10,solidSolutionNoiseMode) : nullptr)
        /* init */,stackingFault(stackingFaultNoiseMode? new StackingFaultNoise(noiseFile,mat,gridSize,gridSpacing_SI) : nullptr)
        {
            
        }
        
        Eigen::Array<int,2,1> gridIndex(const Eigen::Matrix<double,2,1>& localPos) const
        {
            return (localPos.array()/gridSpacing).ceil().template cast<int>();
        }
        
        Eigen::Array<int,2,1> periodicGridIndex(const Eigen::Array<int,2,1>& gi) const
        {
            const Eigen::Array<  double,2,1> gd(gi.template cast<double>());
            return gi-(gd/gridSize.template cast<double>()).floor().template cast<int>()*gridSize;
        }
        
        Eigen::Array<int,2,1> periodicGridIndex(const Eigen::Matrix<double,2,1>& localPos) const
        {
            const Eigen::Array<int,2,1> gi(gridIndex(localPos));
            return periodicGridIndex(gi);
        }
        
        int linearIndex(const int& i,const int& j) const
        {
            return gridSize(1)*i+j;
        }
        
        Eigen::Matrix<double,2,1> gridPos(const Eigen::Array<int,2,1>& gridIdx) const
        {
            return (gridIdx.template cast<double>()*gridSpacing).matrix();
        }
        
        std::tuple<double,double,double> gridInterp(const Eigen::Matrix<double,2,1>& localPos, const bool Debugflag)
        {   // Added by Hyunsoo (hyunsol@g.clemson.edu)
                                    
                /* ************ used parameters ****************
                 <x0,y1>(w2),Ind2     <x1,y1>(w3),Ind3
                 *------------------*
                 |   s2   |   s3     |
                 |        |(localPos)|
                 |________+_________ |
                 |        |          |
                 |        |          |
                 |   s0   |   s1     |
                 *------------------*
                 <x0,y0>(w0),Ind0     <x1,y0>(w1),Ind1
                 ______________________________________________
                 localPos(0), localPos(1) : gauss points (x_g, y_g)
                 x0,x1,y0,y1 : node coordinates
                 Ind0,Ind1,Ind2,Ind3 : local node indices
                 s0,s1,s2,s3 : fraction of areas
                 w0,w1,w2,w3 : weight
                 */
                
                // Get the indices of the grid that contains the gauss point
                // grid_num_x, grid_num_y always find "Ind3" indices due to the usage of "ceil" function)
        //        const int global_num_x = ceil( localPos(0)/gridSpacing(0) ); // Mesh (grid) index, x/dx, grid_num_x is equivalent to "i" in GlidePlaneActor.cpp
        //        const int global_num_y = ceil( localPos(1)/gridSpacing(1) ); // Mesh (grid) index, y/dx, grid_num_y is equivalent to "j" in GlidePlaneActor.cpp
                
        //        const int grid_num_x = global_num_x-std::floor((global_num_x*1.0) /gridSize(0)) * gridSize(0); // periodic condition (0 != 255) & (0 = 256), ex) 256 - (256.0/256 * 256) = 0
        //        const int grid_num_y = global_num_y-std::floor((global_num_y*1.0) /gridSize(1)) * gridSize(1); // periodic condition (0 != 255) & (0 = 256)

                const auto global_num(gridIndex(localPos));
                const int& global_num_x(global_num(0));
                const int& global_num_y(global_num(1));
                const auto grid_num(periodicGridIndex(global_num));
                const int& grid_num_x(grid_num(0));
                const int& grid_num_y(grid_num(1));
                
                // if the periodic condition algorithm fails and finds the grid points outside of the plane, it stops the program
                if (grid_num_x > gridSize(0)-1 || grid_num_y > gridSize(1)-1 || grid_num_x < 0 || grid_num_y < 0)
                {
                    throw std::runtime_error("grid index is out of range!");
                }
                // By default, the previous_x is the grid point located right before the grid point that is found.
                const int previous_x = (grid_num_x-1) >= 0 ? grid_num_x-1 : gridSize(0) - 1; // periodic condition (0 != 255) & (-1 = 255), if previous_x = -1, then previous_x = 255
                const int previous_y = (grid_num_y-1) >= 0 ? grid_num_y-1 : gridSize(1) - 1;

                
                // Convert the local indices to global indices;
        //        const int Ind3 = gridSize(1)*gridSize(2)*grid_num_x+gridSize(2)*grid_num_y; // (i = 0 ~ 255, gridSize = 1 ~ 256), so i = grid_num_x =
        //        const int Ind2 = gridSize(1)*gridSize(2)*previous_x+gridSize(2)*grid_num_y;
        //        const int Ind1 = gridSize(1)*gridSize(2)*grid_num_x+gridSize(2)*previous_y;
        //        const int Ind0 = gridSize(1)*gridSize(2)*previous_x+gridSize(2)*previous_y;     // gridsSize(1)*gridSize(2) -> decalre variable
                
                const int Ind3 = linearIndex(grid_num_x,grid_num_y); // (i = 0 ~ 255, gridSize = 1 ~ 256), so i = grid_num_x =
                const int Ind2 = linearIndex(previous_x,grid_num_y);
                const int Ind1 = linearIndex(grid_num_x,previous_y);
                const int Ind0 = linearIndex(previous_x,previous_y);     // gridsSize(1)*gridSize(2) -> decalre variable


                const Eigen::Array<double,2,1> x1y1( (Eigen::Array<double,2,1>()<<global_num_x*gridSpacing(0), global_num_y*gridSpacing(1)).finished() ); // copy the matrix
                const Eigen::Array<double,2,1> x1y0( (Eigen::Array<double,2,1>()<<global_num_x*gridSpacing(0), (global_num_y-1)*gridSpacing(1)).finished() );  // since it is periodic, it will give us correct area
                const Eigen::Array<double,2,1> x0y1( (Eigen::Array<double,2,1>()<<(global_num_x-1)*gridSpacing(0), global_num_y*gridSpacing(1)).finished() );
                const Eigen::Array<double,2,1> x0y0( (Eigen::Array<double,2,1>()<<(global_num_x-1)*gridSpacing(0), (global_num_y-1)*gridSpacing(1)).finished() );
                
                // calculate the fraction areas and store the area values in a vector (s0, s1, s2, s3)
                const Eigen::Array<double,4,1> fracArea( (Eigen::Array<double,4,1>()<<
                                                          std::abs((localPos(0)-x0y0(0))*(localPos(1)-x0y0(1))), std::abs((localPos(0)-x1y0(0))*(localPos(1)-x1y0(1))),  // (localpos.array() - x0y0).prod()
                                                          std::abs((localPos(0)-x0y1(0))*(localPos(1)-x0y1(1))), std::abs((localPos(0)-x1y1(0))*(localPos(1)-x1y1(1)))).finished() );
                
                // get the weight; w0, w1, w2, w3
                const double totalArea = gridSpacing(0)*gridSpacing(1);
                const Eigen::Matrix<double,4,1> weight( (Eigen::Matrix<double,4,1>()<<fracArea(3)/totalArea, fracArea(2)/totalArea, fracArea(1)/totalArea, fracArea(0)/totalArea).finished() );
                
                const Eigen::Matrix<double,4,1> sfNoise(stackingFault? (Eigen::Matrix<double,4,1>()<< stackingFault->operator[](Ind0), stackingFault->operator[](Ind1),
                                                                                    stackingFault->operator[](Ind2), stackingFault->operator[](Ind3)).finished() : Eigen::Matrix<double,4,1>::Zero());
                const double effsfNoise(weight.dot(sfNoise));  //weight(0)*sfNoise(0) + weight(1)*sfNoise(1) + weight(2)*sfNoise(2) + weight(3)*sfNoise(3);
                
                const Eigen::Matrix<double,4,1> solNoiseXZ(solidSolution? (Eigen::Matrix<double,4,1>()<< solidSolution->operator[](Ind0)(0), solidSolution->operator[](Ind1)(0),
                                                                                       solidSolution->operator[](Ind2)(0), solidSolution->operator[](Ind3)(0)).finished() : Eigen::Matrix<double,4,1>::Zero());
                const Eigen::Matrix<double,4,1> solNoiseYZ(solidSolution? (Eigen::Matrix<double,4,1>()<< solidSolution->operator[](Ind0)(1), solidSolution->operator[](Ind1)(1),
                                                                                       solidSolution->operator[](Ind2)(1), solidSolution->operator[](Ind3)(1)).finished() : Eigen::Matrix<double,4,1>::Zero());
                const double effsolNoiseXZ(weight.dot(solNoiseXZ));
                const double effsolNoiseYZ(weight.dot(solNoiseYZ));

            if (Debugflag)
            {
                std::cout << "grid_num_x = " << grid_num_x << std::endl
                << "grid_num_y = " << grid_num_y << std::endl
                << "previous_x = "  << previous_x << std::endl
                << "previous_y = "  << previous_y << std::endl
                << "Ind0 = " << Ind0 << std::endl
                << "Ind1 = " << Ind1 << std::endl
                << "Ind2 = " << Ind2 << std::endl
                << "Ind3 = " << Ind3 << std::endl
                << "x0y0 = " << x0y0.transpose() << std::endl
                << "x1y0 = " << x1y0.transpose() << std::endl
                << "x0y1 = " << x0y1.transpose() << std::endl
                << "x1y1 = " << x1y1.transpose() << std::endl
                << "fracArea = " << fracArea.transpose() << std::endl
                << "weight = " << weight.transpose() << std::endl;
                //                              << "sfNoise = " << sfNoise.transpose() << std::endl;
            }
            
                return std::make_tuple(effsolNoiseXZ,effsolNoiseYZ,effsfNoise);
        }
        
        
    };


}
#endif

