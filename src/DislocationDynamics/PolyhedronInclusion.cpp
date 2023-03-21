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

    /**********************************************************************/
    template <int dim>
    typename PolyhedronInclusion<dim>::MatrixDim PolyhedronInclusion<dim>::stress(const VectorDim& x) const
    {
        
        if(this->eTNorm>FLT_EPSILON)
        {
//            const VectorDim r(x-C);
//            const double R2(r.squaredNorm());
//            const double R(sqrt(R2));
//
//            if(R>a)
//            {
//                const double R3=std::pow(R,3);
//                const double R4=std::pow(R,4);
//                const double a2R2=std::pow(a,2)/R2;
//                const double a3R3=std::pow(a,3)/R3;
//
//                const VectorDim pTr=this->pT*r;
//                const double pTrr=pTr.dot(r);
//                const double pTt=this->pT.trace();
//                return a3R3/2.0/(1-this->nu)*( (10.0*(1.0-2.0*this->nu)+6.0*a2R2)/15.0*this->pT
//                                        +(2.0*this->nu-2.0*a2R2)/R2*(pTr*r.transpose()+r*pTr.transpose())
//                                        +((3.0*a2R2-5.0*(1.0-2.0*this->nu))/15.0*pTt + (1.0-2.0*this->nu-a2R2)/R2*pTrr)*MatrixDim::Identity()
//                                        +(-(5.0-7.0*a2R2)/R4*pTrr+(1.0-a2R2)/R2*pTt)*r*r.transpose()
//                                        );
//            }
//            else
//            {
//                return 2.0*this->mu*((L+this->nu/(1.0-2.0*this->nu)*(2.0*M+3.0*L))*this->eT.trace()*MatrixDim::Identity()+2.0*M*this->eT)-this->pT;
//            }
            return MatrixDim::Zero();
        }
        else
        {
            return MatrixDim::Zero();
        }
    }

    template class PolyhedronInclusion<3>;

}
#endif
