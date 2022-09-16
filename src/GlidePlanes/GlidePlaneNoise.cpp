/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpF@ucla.edu>.
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneNoise_cpp
#define model_GlidePlaneNoise_cpp

#include <GlidePlaneNoise.h>

#include <cmath>
#include <random>
#include <Eigen/Dense>


#ifdef _MODEL_GLIDE_PLANE_NOISE_GENERATOR_
#include <fftw3.h>
#include <boost/math/special_functions/bessel.hpp>
#endif

#include <PolycrystallineMaterialBase.h>
#include <TerminalColors.h>
#include <UniformPeriodicGrid.h>
#include <filesystem>

namespace model
{

    int SolidSolutionNoiseReader::LittleEndian()
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

    float SolidSolutionNoiseReader::ReverseFloat( const float inFloat )
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

    double SolidSolutionNoiseReader::ReverseDouble( const double inDouble )
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

    std::pair<typename SolidSolutionNoiseReader::GridSizeType,typename SolidSolutionNoiseReader::GridSpacingType> SolidSolutionNoiseReader::Read_dimensions(const char *fname)
    {
        int NX, NY, NZ;
        double DX, DY, DZ;
        char line[200];
        FILE *InFile=fopen(fname,"r");
        
        if (InFile == NULL)
        {
            fprintf(stderr, "Can't open noise file %s\n",fname);
            exit(1);
        }
        
        for(int i=0;i<5;i++)
        {
            fgets(line, 200, InFile);
        }
        fscanf(InFile, "%s %lf %lf %lf\n", line, &(DX), &(DY), &(DZ));
        fscanf(InFile, "%s %d %d %d\n", line, &(NX), &(NY), &(NZ));
        return std::make_pair((GridSizeType()<<NX,NY).finished(),(GridSpacingType()<<DX,DY).finished());
    }



    void SolidSolutionNoiseReader::Read_noise_vtk(const char *fname, REAL_SCALAR *Noise, int Nr,const double& MSS)
    {
        char line[200];
        double temp;
        FILE *InFile=fopen(fname,"r");
        
        if (InFile == NULL)
        {
            fprintf(stderr, "Can't open noise file %s\n",fname);
            exit(1);
        }
        
        for(int i=0;i<10;i++)
        {
            fgets(line, 200, InFile);
        }
        
        if(LittleEndian()) // if machine works with LittleEndian
        {
            for(int ind=0;ind<Nr;ind++)
            {
                fread(&temp, sizeof(double), 1, InFile);
                Noise[ind] = MSS*REAL_SCALAR(ReverseDouble(temp));
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

    SolidSolutionNoiseReader::SolidSolutionNoiseReader(const std::string& noiseFile,const PolycrystallineMaterialBase& mat,
                                                       const typename SolidSolutionNoiseReader::GridSizeType& _gridSize, const typename SolidSolutionNoiseReader::GridSpacingType& _gridSpacing_A)
    {
        std::cout<<"Reading SolidSolutionNoise files"<<std::endl;

        const std::string fileName_xz(std::filesystem::path(noiseFile).parent_path().string()+"/"+TextFileParser(noiseFile).readString("solidSolutionNoiseFile_xz",true));
        const auto gridSize_xz(Read_dimensions(fileName_xz.c_str()));
        const std::string fileName_yz(std::filesystem::path(noiseFile).parent_path().string()+"/"+TextFileParser(noiseFile).readString("solidSolutionNoiseFile_yz",true));
        const auto gridSize_yz(Read_dimensions(fileName_yz.c_str()));
        const double MSSS_SI(TextFileParser(mat.materialFile).readScalar<double>("MSSS_SI",true));
        const double MSS(std::sqrt(MSSS_SI)/mat.mu_SI);
        
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

    const typename SolidSolutionNoise::NoiseContainerType& SolidSolutionNoise::noiseVector() const
    {
        return *this;
    }

    typename SolidSolutionNoise::NoiseContainerType& SolidSolutionNoise::noiseVector()
    {
        return *this;
    }

    SolidSolutionNoise::SolidSolutionNoise(const std::string& noiseFile,const PolycrystallineMaterialBase& mat,
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

    StackingFaultNoise::StackingFaultNoise(const std::string&, // noiseFile
                                           const PolycrystallineMaterialBase& mat,
                                           const NoiseTraitsBase::GridSizeType& gridSize,
                                           const NoiseTraitsBase::GridSpacingType& gridSpacing_SI)
    {
        
        std::cout<<greenBoldColor<<"Creating StackingFaultNoise"<<defaultColor<<std::endl;
        
        const double isfEnergyDensityMEAN(TextFileParser(mat.materialFile).readScalar<double>("isfEnergyDensityMEAN_SI",true)/(mat.mu_SI*mat.b_SI));
        const double isfEnergyDensitySTD(TextFileParser(mat.materialFile).readScalar<double>("isfEnergyDensitySTD_SI",true)/std::sqrt(gridSpacing_SI(0)*gridSpacing_SI(1))/(mat.mu_SI*mat.b_SI));
        
        std::normal_distribution<double> distribution (isfEnergyDensityMEAN,isfEnergyDensitySTD);
        
        const size_t N(gridSize.array().prod());
        this->reserve(N);
        for(size_t k=0;k<N;++k)
        {// J/m^2 = N/m = Pa*m
            this->push_back(distribution(generator)-isfEnergyDensityMEAN);
            //                this->push_back(distribution(generator));
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

    GlidePlaneNoise::GlidePlaneNoise(const std::string& noiseFile,const PolycrystallineMaterialBase& mat) :
    /* init */ UniformPeriodicGrid<2>(TextFileParser(noiseFile).readMatrix<int,1,2>("gridSize",true),TextFileParser(noiseFile).readMatrix<double,1,2>("gridSpacing_SI",true)/mat.b_SI)
    /* init */,solidSolutionNoiseMode(TextFileParser(noiseFile).readScalar<int>("solidSolutionNoiseMode"))
    /* init */,stackingFaultNoiseMode(TextFileParser(noiseFile).readScalar<int>("stackingFaultNoiseMode"))
    /* init */,solidSolution(solidSolutionNoiseMode? new SolidSolutionNoise(noiseFile,mat,gridSize,this->gridSpacing*mat.b_SI*1.0e10,solidSolutionNoiseMode) : nullptr)
    /* init */,stackingFault(stackingFaultNoiseMode? new StackingFaultNoise(noiseFile,mat,gridSize,this->gridSpacing*mat.b_SI) : nullptr)
    {
        
    }

    int GlidePlaneNoise::storageIndex(const int& i,const int& j) const
    {/*!\param[in] localPos the  position vector on the grid
      * \returns The grid index periodically wrapped within the gridSize bounds
      */
        return this->gridSize(1)*i+j;
    }

    std::tuple<double,double,double> GlidePlaneNoise::gridInterp(const Eigen::Matrix<double,2,1>& localPos) const
    {   // Added by Hyunsoo (hyunsol@g.clemson.edu)
        
        const auto idxAndWeights(this->posToPeriodicCornerIdxAndWeights(localPos));
        double effsolNoiseXZ(0.0);
        double effsolNoiseYZ(0.0);
        double effsfNoise(0.0);
        for(size_t p=0;p<idxAndWeights.first.size();++p)
        {
            const int storageID(storageIndex(idxAndWeights.first[p](0),idxAndWeights.first[p](1)));
            if(solidSolution)
            {
                effsolNoiseXZ+=solidSolution->operator[](storageID)(0)*idxAndWeights.second[p];
                effsolNoiseYZ+=solidSolution->operator[](storageID)(1)*idxAndWeights.second[p];
            }
            if(stackingFault)
            {
                effsfNoise+=stackingFault->operator[](storageID)*idxAndWeights.second[p];
            }
        }
        
        return std::make_tuple(effsolNoiseXZ,effsolNoiseYZ,effsfNoise);
    }

    std::tuple<double,double,double> GlidePlaneNoise::gridVal(const Eigen::Array<int,2,1>& idx) const
    {   // Added by Hyunsoo (hyunsol@g.clemson.edu)
        const Eigen::Array<int,2,1> pidx(this->idxToPeriodicIdx(idx));
        const int storageID(storageIndex(pidx(0),pidx(1)));
        double effsolNoiseXZ(0.0);
        double effsolNoiseYZ(0.0);
        double effsfNoise(0.0);
        if(solidSolution)
        {
            effsolNoiseXZ=solidSolution->operator[](storageID)(0);
            effsolNoiseYZ=solidSolution->operator[](storageID)(1);
        }
        if(stackingFault)
        {
            effsfNoise=stackingFault->operator[](storageID);
        }
        return std::make_tuple(effsolNoiseXZ,effsolNoiseYZ,effsfNoise);
    }

}
#endif

