/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PolyhedronInclusion_H_
#define model_PolyhedronInclusion_H_


#include <Eigen/Dense>
//#include <Material.h>
#include <StaticID.h>

// http://solidmechanics.org/text/Chapter5_4/Chapter5_4.htm
#include <SlipSystem.h>
#include <SecondPhase.h>
#include <EshelbyInclusionBase.h>
#include <PolyhedronInclusionNodeIO.h>
#include <Plane.h>
#include <VoigtTraits.h>

namespace model
{
    
    template <int dim>
    class PolyhedronInclusion : public EshelbyInclusionBase<dim>
    /*                      */, private std::map<size_t,Plane<dim>>
    {
        
    public:

        static constexpr int voigtSize=SymmetricVoigtTraits<dim>::voigtSize;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        typedef Eigen::Matrix<double,voigtSize,voigtSize> MatrixVoigt;
        
    private:
        static const SymmetricVoigtTraits<dim> voigtTraits;

        static std::map<size_t,std::vector<std::pair<size_t,VectorDim>>> getFaces(const std::map<size_t,PolyhedronInclusionNodeIO<dim>>& nodesMap,
                                                                                  const std::map<size_t,std::vector<size_t>>& faceIDs);
        
        static double PHI_ij(int i, int j, double a, double b, double lm, double lp, const VectorDim& Svnorm, const VectorDim& Vnorm, const VectorDim& Vdir);
        static double Phi_u_II_a(double a, double b, double le);
        static double Phi_u_II_b(double a, double b, double le);
        static double Phi_u_II_le(double a, double b, double le);
        
        double PHI_ij(int i, int j,const VectorDim& x) const;


    public:
        
        const std::map<size_t,PolyhedronInclusionNodeIO<dim>>& nodes;
        const std::map<size_t,std::vector<std::pair<size_t,VectorDim>>> faces;
        
        PolyhedronInclusion(const std::map<size_t,PolyhedronInclusionNodeIO<dim>>& nodesMap,
                            const std::map<size_t,std::vector<size_t>>& faceIDs,
                         const MatrixDim& _eT,
                         const double& _nu,
                         const double& _mu,
                         const double& _mobilityReduction,
                         const int& _phaseID,
                         const std::shared_ptr<SecondPhase<dim>>& sph);
  
        
        bool contains(const VectorDim& x) const override;
        MatrixDim stress(const VectorDim& x) const override;
        MatrixDim elasticStrain(const VectorDim& x) const ;
        MatrixDim strain(const VectorDim& x) const ;
        double eshelbyTensorComponent(const int&i,const int&j,const int&k,const int&l,const VectorDim& x) const;
        MatrixVoigt eshelbyTensorVoigt(const VectorDim& x) const ;
        const std::map<size_t,Plane<dim>>& planes() const;
        
    };
}
#endif
