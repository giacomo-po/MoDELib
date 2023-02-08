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

#include <DislocationFieldBase.h>
#include <PolycrystallineMaterialBase.h>
#include <UniformPeriodicGrid.h>

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
        
        
        
        static int LittleEndian();
        
        static float ReverseFloat( const float inFloat );
        
        static double ReverseDouble( const double inDouble );
        
        static std::pair<GridSizeType,GridSpacingType> Read_dimensions(const char *fname);
        
        static void Read_noise_vtk(const char *fname, REAL_SCALAR *Noise, int Nr,const double& MSS);
        
        SolidSolutionNoiseReader(const std::string& noiseFile,const PolycrystallineMaterialBase& mat,
                                 const GridSizeType& _gridSize, const GridSpacingType& _gridSpacing_A);
        
        
    };

    class SolidSolutionNoise : public NoiseTraits<2>::NoiseContainerType
    {
        typedef typename NoiseTraitsBase::REAL_SCALAR REAL_SCALAR;
        typedef typename NoiseTraitsBase::COMPLEX COMPLEX;
        typedef typename NoiseTraitsBase::GridSizeType GridSizeType;
        typedef typename NoiseTraitsBase::GridSpacingType GridSpacingType;
        typedef typename NoiseTraits<2>::NoiseType NoiseType;
        typedef typename NoiseTraits<2>::NoiseContainerType NoiseContainerType;
        
        
        const NoiseContainerType& noiseVector() const;
        
        NoiseContainerType& noiseVector();

    public:
        
        const GridSizeType gridSize;
        const GridSpacingType gridSpacing_A;
        
        SolidSolutionNoise(const std::string& noiseFile,const PolycrystallineMaterialBase& mat,
                           const GridSizeType& _gridSize, const GridSpacingType& _gridSpacing_A, const int& solidSolutionNoiseMode);
        
    };

    struct SolidSolutionNoiseGenerator: public NoiseTraits<2>::NoiseContainerType
    {

        typedef typename NoiseTraits<2>::REAL_SCALAR REAL_SCALAR;
        typedef typename NoiseTraits<2>::COMPLEX COMPLEX;
        typedef typename NoiseTraits<2>::GridSizeType GridSizeType;
        typedef typename NoiseTraitsBase::GridSpacingType GridSpacingType;
        typedef typename NoiseTraits<2>::NoiseType NoiseType;
        typedef typename NoiseTraits<2>::NoiseContainerType NoiseContainerType;
        
        int NX, NY, NZ;
        REAL_SCALAR DX, DY, DZ;
        REAL_SCALAR a;
        REAL_SCALAR a_cai;
        int seed;
        //    int flag;
        REAL_SCALAR LX, LY, LZ;
        REAL_SCALAR DV;
        int NR;
        int NK;
        REAL_SCALAR Norm;

        SolidSolutionNoiseGenerator(const std::string& noiseFile,const PolycrystallineMaterialBase& mat,
                                    const GridSizeType& _gridSize, const GridSpacingType& _gridSpacing_A);
        
        
#ifdef _MODEL_GLIDE_PLANE_NOISE_GENERATOR_
        // Cai doubly-convoluted spreading function in Fourier space
        static REAL_SCALAR Wk_Cai(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz, REAL_SCALAR a) ;

        // Cai spreading function
        static REAL_SCALAR W_Cai(REAL_SCALAR r2, REAL_SCALAR a) ;
 
        static REAL_SCALAR W_t_Cai(REAL_SCALAR r2, REAL_SCALAR a) ;
        // normalized auto-correlation function in Fourier space for sigma_xy
        REAL_SCALAR S_xy_k(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz) const;

        // normalized auto-correlation function in Fourier space for sigma_xz
        REAL_SCALAR S_xz_k(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz) const;

        // normalized auto-correlation function in Fourier space for sigma_yz
        REAL_SCALAR S_yz_k(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz) const;
#endif
        
        void Write_field_slice(REAL_SCALAR *F, const char *fname);
        
    };




    struct StackingFaultNoise : public NoiseTraits<1>::NoiseContainerType
    {
        typedef typename NoiseTraits<1>::REAL_SCALAR REAL_SCALAR;
        typedef typename NoiseTraits<1>::GridSizeType GridSizeType;
        typedef typename NoiseTraits<1>::NoiseType NoiseType;
        typedef typename NoiseTraits<1>::NoiseContainerType NoiseContainerType;
        
        std::default_random_engine generator;


        
        StackingFaultNoise(const std::string&, // noiseFile
                           const PolycrystallineMaterialBase& mat,
                           const NoiseTraitsBase::GridSizeType& gridSize,
                           const NoiseTraitsBase::GridSpacingType& gridSpacing_SI);
        
    };


    

    struct GlidePlaneNoise : public UniformPeriodicGrid<2>
    {
        typedef typename Eigen::Matrix<double,3,1> VectorDim;
        typedef typename NoiseTraitsBase::GridSizeType GridSizeType;

        const int solidSolutionNoiseMode;
        const int stackingFaultNoiseMode;

        const std::shared_ptr<SolidSolutionNoise> solidSolution;
        const std::shared_ptr<StackingFaultNoise> stackingFault;
                
        
        GlidePlaneNoise(const std::string& noiseFile,const PolycrystallineMaterialBase& mat);
        
        GridSizeType rowAndColIndices(const int& storageIndex) const;
        
        int storageIndex(const int& i,const int& j) const;
                
        std::tuple<double,double,double> gridInterp(const Eigen::Matrix<double,2,1>& localPos) const;
        
        std::tuple<double,double,double> gridVal(const Eigen::Array<int,2,1>& idx) const;
        
        
    };


}
#endif

