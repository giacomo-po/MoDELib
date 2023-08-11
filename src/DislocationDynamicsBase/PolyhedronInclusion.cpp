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

    template <int dim>
    PolyhedronInclusion<dim>::PolyhedronInclusion(const std::map<size_t,PolyhedronInclusionNodeIO<dim>>& nodesMap,
                                                  const std::map<size_t,std::vector<size_t>>& faceIDs,
                                            const MatrixDim& _eT,
                                            const double& _nu,
                                            const double& _mu,
                                            const double& _mobilityReduction,
                                            const int& _phaseID,
                                            const std::shared_ptr<SecondPhase<dim>>& sph) :
    /* init */ EshelbyInclusionBase<dim>(_eT,_nu,_mu,_mobilityReduction,_phaseID,sph)
    /* init */,nodes(nodesMap)
    /* init */,faces(getFaces(nodesMap,faceIDs))
    {
        std::cout<<"Creating PolyhedronInclusion "<<this->sID<<" (type "<<this->phaseID<<"):\n eT="<<this->eT<<std::endl;
        
        for(const auto& face : faces)
        {
            VectorDim C(VectorDim::Zero()); // face center
            for(const auto& nodePair : face.second)
            {
//                const auto nodeIter(nodes.find(nodeID));
                C+=nodePair.second;
//                if(nodeIter!=nodes.end())
//                {
//                    C+=nodeIter->second.P;
//                }
//                else
//                {
//                    throw std::runtime_error("PolyhedronInclusion: nodeID not found in nodes.");
//                }
            }
            C/=face.second.size();
            
            VectorDim nA(VectorDim::Zero()); // normal
//            const VectorDim P0(nodes.find(face.second.front())->second.P);
            const VectorDim P0(face.second.front().second);
            for(size_t k=0;k<face.second.size();++k)
            {
                const size_t k1(k<face.second.size()-1? k+1 : 0);
//                const VectorDim Pk(nodes.find(face.second[k])->second.P);
//                const VectorDim Pk1(nodes.find(face.second[k1])->second.P);
                const VectorDim Pk(face.second[k].second);
                const VectorDim Pk1(face.second[k1].second);
                nA+= 0.5*(Pk-P0).cross(Pk1-Pk);
            }

            Plane<3> plane(C,nA);
            for(const auto& nodePair : face.second)
            {
//                const auto nodeIter(nodes.find(nodeID));
//                if(nodeIter!=nodes.end())
//                {
                    if(!plane.contains(nodePair.second))
                    {
                        throw std::runtime_error("PolyhedronInclusion: face plane does not include face vertex.");
                    }
//                }
            }
            this->emplace(face.first,plane);
        }
    }

    template <int dim>
    const std::map<size_t,Plane<dim>>& PolyhedronInclusion<dim>::planes() const
    {
        return *this;
    }

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
std::map<size_t,std::vector<std::pair<size_t,typename PolyhedronInclusion<dim>::VectorDim>>> PolyhedronInclusion<dim>::getFaces(const std::map<size_t,PolyhedronInclusionNodeIO<dim>>& nodesMap,
                                                                                                                                const std::map<size_t,std::vector<size_t>>& faceIDs)
{
    std::map<size_t,std::vector<std::pair<size_t,VectorDim>>> temp;
    for(const auto& fID : faceIDs)
    {
        std::vector<std::pair<size_t,VectorDim>> tempV;
        for(const auto& nID : fID.second)
        {
            const auto nodeIter(nodesMap.find(nID));
            if(nodeIter!=nodesMap.end())
            {
                tempV.emplace_back(nID,nodeIter->second.P);
            }
            else
            {
                throw std::runtime_error("PolyhedronInclusion node not found.");
            }
        }
        temp.emplace(fID.first,tempV);
    }
    return temp;
}

template <int dim>
double PolyhedronInclusion<dim>::Phi_u_II_a(double a, double b, double le)
{
//    double result;
    if (a != 0 && b != 0 && le != 0) {
        return -a / abs(a) * atan(le / b) + asin(a * le * abs(b) / (b * sqrt((a * a + b * b) * (b * b + le * le))));
    }
    else if (a == 0) {
        return 0.0;
    }
    else if (a != 0) {

        return 0.0;
    }
    else { // a, b, le equals zero
        return 0.0;
    }
}

template <int dim>
double PolyhedronInclusion<dim>::Phi_u_II_b(double a, double b, double le)
{
//    double result;
    double dis = sqrt(a * a + b * b + le * le);
    if (a != 0 && b != 0 && le != 0) {
        return (1.0 / (b * b + le * le)) * (
            -le * dis + le * abs(a) + (b * b + le * le) * atanh(le / dis));
    }
    else if (a == 0) {
        if (b != 0) {    // le may equal zero, all includes
            return -le / dis + atanh(le / dis);
        }
        else { // b ==0, le != 0
            return 0.0;
        }
    }
    else if (a != 0) {
        if (b != 0 && le == 0) {
            return 0.0;
        }
        else if (b == 0 && le == 0) {
            return 0.0;
        }
        else if (b == 0 && le != 0) {
            return (-dis + abs(a)) / le + atanh(le / dis);
        }
    }
    else { // a, b, le equals zero
        return 0.0;
    }
}

template <int dim>
double PolyhedronInclusion<dim>::Phi_u_II_le(double a, double b, double le)
{
//    double result;
    if (a != 0 && b != 0 && le != 0) {
        return b * (a * a + b * b + le * le - sqrt(a * a + b * b + le * le) * abs(a)) / ((b * b + le * le) * sqrt(a * a + b * b + le * le));
    }
    else if (a == 0) {
        if (b != 0 || le != 0) {
            return b / sqrt(b * b + le * le);
        }
        else {
            return 0.0;
        }
    }
    else if (a != 0) {
        if (b != 0 && le == 0) {
            return (sqrt(a * a + b * b) - abs(a)) / b;
        }
        else if (b == 0 && le != 0) {
            return 0;
        }
        else if (b == 0 && le == 0) {
            return 0;
        }
    }
}

template <int dim>
double PolyhedronInclusion<dim>::PHI_ij(int i, int j, double a, double b, double lm, double lp, const VectorDim& Svnorm, const VectorDim& Vnorm, const VectorDim& Vdir)
{
//    double result = 0.0;
//    result = -(Svnorm[i]) * (-(Phi_u_II_a(a, b, lp) - Phi_u_II_a(a, b, lm)) * Svnorm[j]
//        - (Phi_u_II_b(a, b, lp) - Phi_u_II_b(a, b, lm)) * Vnorm[j]
//        - (Phi_u_II_le(a, b, lp)) * Vdir[j] + (Phi_u_II_le(a, b, lm)) * Vdir[j]
//        );
    return -(Svnorm[i]) * (-(Phi_u_II_a(a, b, lp) - Phi_u_II_a(a, b, lm)) * Svnorm[j]
                           - (Phi_u_II_b(a, b, lp) - Phi_u_II_b(a, b, lm)) * Vnorm[j]
                           - (Phi_u_II_le(a, b, lp)) * Vdir[j] + (Phi_u_II_le(a, b, lm)) * Vdir[j]
                           );

}

template <int dim>
double PolyhedronInclusion<dim>::PHI_ij(int i, int j, const VectorDim& x) const
{
    double result=0.0;
    for(const auto& face : faces)
    {
        const auto& P0(face.second[0].second);
        const auto& P1(face.second[1].second);
        const auto& P2(face.second[2].second);
        const VectorDim Svnorm((P1-P0).cross(P2-P1).normalized());
            const double a(Svnorm.dot(P0-x));
            for(size_t e=0;e<face.second.size();++e)
            {
                const size_t e1(e<face.second.size()-1? e+1 : 0);
                const VectorDim& vm(face.second[e].second);
                const VectorDim& vp(face.second[e1].second);
                const VectorDim Vdir((vp-vm).normalized());
                const VectorDim Vnorm(Vdir.cross(Svnorm));
                const double b((vp-x).dot(Vnorm));
                const double lm((vm-x).dot(Vdir));
                const double lp((vp-x).dot(Vdir));
                result+=PHI_ij(i,j,a,b,lm,lp,Svnorm,Vnorm,Vdir);
            }
    }
    return result;
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
//        return PHI_ij(0,0,x);
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
