/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PolyhedronInclusion_cpp_
#define model_PolyhedronInclusion_cpp_


#include <PolyhedronInclusion.h>
#include <numbers>

namespace model
{


//    template <int dim>
//    void PolyhedronInclusion<dim>::addSlipSystems(const std::vector<std::shared_ptr<SlipSystem>>& slipSystems)
//    {
//
//        for(const auto& slipSystem : slipSystems)
//        {
//            if(gammaSurfaceMap.find(slipSystem->gammaSurface.get())==gammaSurfaceMap.end())
//            {// current slipSystem gammaSurface not found
//
//                //                    GammaSurface temp();
//                //
//                //                    gammaSurfaceMap.emaplace(slipSystem->gammaSurface.get(),temp);
//            }
//        }
//
//    }

//    template <int dim>
//    double PolyhedronInclusion<dim>::misfitEnergy(const Eigen::Matrix<double,3,1>& b, const GammaSurface* const matrixGammaSurface) const
//    {
//        const auto iter(gammaSurfaceMap.find(matrixGammaSurface));
//        return iter==gammaSurfaceMap.end()? 0.0 : iter->second(b);
//    }

    /**********************************************************************/
    template <int dim>
    PolyhedronInclusion<dim>::PolyhedronInclusion(const std::map<size_t,PolyhedronInclusionNodeIO<dim>>& nodes_in,
                                                  const std::map<size_t,std::vector<size_t>>& faces_in,
                                            const MatrixDim& _eT,
                                            const double& _nu,
                                            const double& _mu,
                                            const double& _mobilityReduction,
                                            const int& _phaseID,
                                            const std::shared_ptr<SecondPhase<dim>>& sph) :
    /* init */ EshelbyInclusionBase<dim>(_eT,_nu,_mu,_mobilityReduction,_phaseID,sph)
    /* init */,nodes(nodes_in)
    /* init */,faces(faces_in)
    {
        std::cout<<"Creating PolyhedronInclusion "<<this->sID<<" (type "<<this->phaseID<<"):\n eT="<<this->eT<<std::endl;
        
        for(const auto& face : faces)
        {
            VectorDim C(VectorDim::Zero()); // face center
            for(const auto& nodeID : face.second)
            {
                const auto nodeIter(nodes.find(nodeID));
                if(nodeIter!=nodes.end())
                {
                    C+=nodeIter->second.P;
                }
                else
                {
                    throw std::runtime_error("PolyhedronInclusion: nodeID not found in nodes.");
                }
            }
            C/=face.second.size();
            
            VectorDim nA(VectorDim::Zero()); // face center
            const VectorDim P0(nodes.find(face.second.front())->second.P);
            for(size_t k=0;k<face.second.size();++k)
            {
                const size_t k1(k<face.second.size()-1? k+1 : 0);
                const VectorDim Pk(nodes.find(face.second[k])->second.P);
                const VectorDim Pk1(nodes.find(face.second[k1])->second.P);
                nA+= 0.5*(Pk-P0).cross(Pk1-Pk);
            }

            Plane<3> plane(C,nA);
            for(const auto& nodeID : face.second)
            {
                const auto nodeIter(nodes.find(nodeID));
                if(nodeIter!=nodes.end())
                {
                    if(!plane.contains(nodeIter->second.P))
                    {
                        throw std::runtime_error("PolyhedronInclusion: face plane does not include face vertex.");
                    }
                }
            }
            
            this->emplace(face.first,plane);
        }
        

        
    }

template <int dim>
const std::map<size_t,Plane<dim>>& PolyhedronInclusion<dim>::planes() const
{
    return *this;
}

    /**********************************************************************/
    template <int dim>
    bool PolyhedronInclusion<dim>::contains(const VectorDim& x) const
    {
        bool contained(true);
        for(const auto& plane : this->planes())
        {
            contained=(contained && plane.second.isAbove(x));
        }
        return contained;
    }

    template <int dim>
    double PolyhedronInclusion<dim>::eshelbyTensorComponent(const int&i,const int&j,const int&k,const int&l,const VectorDim& x) const
    {
//        return  0.125/(std::numbers::pi*(1.0-this->nu))*(Psi_ijkl_a[voigtIndex(i,j)][voigtIndex(k,l)]
//                                                         -2.0*this->nu*d[2][2]*Phi_ij_a[voigtIndex(i,j)]
//                                                         -(1.0-this->nu)*(Phi_ij_a[voigtIndex(i,l)]*d[j][k]
//                                                                          +Phi_ij_a[voigtIndex(j,k)]*d[i][l]
//                                                                          +Phi_ij_a[voigtIndex(j,l)]*d[i][k]
//                                                                          +Phi_ij_a[voigtIndex(i,k)]*d[j][l]));
        return 0.0;
    }

    template <int dim>
    Eigen::Matrix<double,PolyhedronInclusion<dim>::voigtSize,PolyhedronInclusion<dim>::voigtSize> PolyhedronInclusion<dim>::eshelbyTensorVoigt(const VectorDim& x) const
    {
        
        Eigen::Matrix<double,voigtSize,voigtSize> temp(Eigen::Matrix<double,voigtSize,voigtSize>::Zero());
        
        temp(voigtTraits.voigtIndex(0,0),voigtTraits.voigtIndex(0,0))=eshelbyTensorComponent(0,0,0,0,x);
        temp(voigtTraits.voigtIndex(0,0),voigtTraits.voigtIndex(1,1))=eshelbyTensorComponent(0,0,1,1,x);
        temp(voigtTraits.voigtIndex(0,0),voigtTraits.voigtIndex(2,2))=eshelbyTensorComponent(0,0,2,2,x);
        
        temp(voigtTraits.voigtIndex(1,1),voigtTraits.voigtIndex(0,0))=eshelbyTensorComponent(1,1,0,0,x);
        temp(voigtTraits.voigtIndex(1,1),voigtTraits.voigtIndex(1,1))=eshelbyTensorComponent(1,1,1,1,x);
        temp(voigtTraits.voigtIndex(1,1),voigtTraits.voigtIndex(2,2))=eshelbyTensorComponent(1,1,2,2,x);

        temp(voigtTraits.voigtIndex(2,2),voigtTraits.voigtIndex(0,0))=eshelbyTensorComponent(2,2,0,0,x);
        temp(voigtTraits.voigtIndex(2,2),voigtTraits.voigtIndex(1,1))=eshelbyTensorComponent(2,2,1,1,x);
        temp(voigtTraits.voigtIndex(2,2),voigtTraits.voigtIndex(2,2))=eshelbyTensorComponent(2,2,2,2,x);

        temp(voigtTraits.voigtIndex(0,1),voigtTraits.voigtIndex(0,1))=2.0*eshelbyTensorComponent(0,1,0,1,x);
        temp(voigtTraits.voigtIndex(0,2),voigtTraits.voigtIndex(0,2))=2.0*eshelbyTensorComponent(0,2,0,2,x);
        temp(voigtTraits.voigtIndex(1,2),voigtTraits.voigtIndex(1,2))=2.0*eshelbyTensorComponent(1,2,1,2,x);

        return temp;
    }


    template <int dim>
    typename PolyhedronInclusion<dim>::MatrixDim PolyhedronInclusion<dim>::strain(const VectorDim& x) const
    {
        return this->eTNorm>FLT_EPSILON? voigtTraits.v2m(eshelbyTensorVoigt(x)*voigtTraits.m2v(this->eT,true),true) : MatrixDim::Zero();
    }

    template <int dim>
    typename PolyhedronInclusion<dim>::MatrixDim PolyhedronInclusion<dim>::elasticStrain(const VectorDim& x) const
    {
        if(contains(x))
        {
            return strain(x)-this->eT;
        }
        else
        {
            return strain(x);
        }
    }

    template <int dim>
    typename PolyhedronInclusion<dim>::MatrixDim PolyhedronInclusion<dim>::stress(const VectorDim& x) const
    {
        const MatrixDim elStrain(elasticStrain(x));
        return this->lambda*elStrain.trace()*MatrixDim::Identity()+2.0*this->mu*elStrain;
    }

    template <int dim>
    const SymmetricVoigtTraits<dim> PolyhedronInclusion<dim>::voigtTraits=SymmetricVoigtTraits<dim>((typename SymmetricVoigtTraits<dim>::VoigtSizeMatrixType()<<0,0,1,1,2,2,1,2,0,2,0,1).finished());
    //const typename PolyhedronInclusion<dim>::VoigtSizeMatrixType PolyhedronInclusion<dim>::voigtOrder=(PolyhedronInclusion<dim>::VoigtSizeMatrixType()<<0,0,0,1,0,2,1,1,1,2,2,2).finished();

    template class PolyhedronInclusion<3>;

}
#endif
