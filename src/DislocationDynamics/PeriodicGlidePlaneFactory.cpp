/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicGlidePlaneFactory_cpp_
#define model_PeriodicGlidePlaneFactory_cpp_

#include <Eigen/Dense>
#include <PeriodicGlidePlaneFactory.h>

namespace model
{

    template<int dim>
    PeriodicGlidePlaneFactory<dim>::PeriodicGlidePlaneFactory(const Polycrystal<dim>& poly_in,
                              GlidePlaneFactory<dim>& glidePlaneFactory_in) :
    /* init */ Lattice<dim>(get_B(poly_in),Eigen::Matrix<double,dim,dim>::Identity())
    /* init */,poly(poly_in)
    /* init */,glidePlaneFactory(glidePlaneFactory_in)
    /* init */,N(get_N(poly,this->latticeBasis))
    {
    }
    
    /**********************************************************************/
    //Excluded this during cpp implementation....Was not being used
    // template<int dim>
    // typename PeriodicGlidePlaneFactory<dim>::PeriodicGlidePlaneSharedPtrType  PeriodicGlidePlaneFactory<dim>::get(const GlidePlaneType& plane)
    // {
    //     return BaseType::getFromKey(periodicPlaneKey(plane.key));
    // }
    
    /**********************************************************************/
    template<int dim>
    typename PeriodicGlidePlaneFactory<dim>::PeriodicGlidePlaneSharedPtrType  PeriodicGlidePlaneFactory<dim>::get(const GlidePlaneKeyType& temp)
    {
//        return BaseType::get(periodicPlaneKey(temp)); // unique key for plane family (reference plane)
        return BaseType::getFromKey(temp);
    }
    
    template<int dim>
    typename PeriodicGlidePlaneFactory<dim>::BaseType&  PeriodicGlidePlaneFactory<dim>::periodicGlidePlanes()
    {
        return *this;
    }
    
    template<int dim>
    const typename PeriodicGlidePlaneFactory<dim>::BaseType&  PeriodicGlidePlaneFactory<dim>::periodicGlidePlanes() const
    {
        return *this;
    }
    
    template<int dim>
    Eigen::Matrix<double,dim,dim>  PeriodicGlidePlaneFactory<dim>::get_B(const Polycrystal<dim>& poly)
    {
        assert(poly.grains().size() == 1 && "Periodic simulations only supported in single crystals");
        const auto &grain(poly.grains().begin()->second);
        const auto &meshRegion(grain.region);
        
        Eigen::Matrix<double, dim, dim> B(Eigen::Matrix<double, dim, dim>::Zero());
        int col = 0;
        for (const auto &pair : meshRegion.parallelFaces())
        {
//            std::cout << "Checking if parallel faces " << pair.first << "<->" << pair.second << " are commensurate" << std::endl;
            const PlanarMeshFace<dim> &face1(*meshRegion.faces().at(pair.first));
            const PlanarMeshFace<dim> &face2(*meshRegion.faces().at(pair.second));
            const VectorDim cc(face1.center() - face2.center());
            const VectorDim ccc(cc.dot(face1.outNormal()) * face1.outNormal());
            if (col < dim)
            {
                B.col(col) = ccc;
                col++;
            }
            const LatticeDirection<dim> ld(grain.latticeDirection(face1.outNormal()));
            const double normRatio(ccc.norm() / ld.cartesian().norm());
            if (std::fabs(std::round(normRatio) - normRatio) > FLT_EPSILON)
            {
                //                            std::cout<<"Face outNormal="<<std::setprecision(15)<<std::scientific<<face1.outNormal().transpose()<<std::endl;
                std::cout << "Mesh in direction " << std::setprecision(15) << std::scientific << ld.cartesian().normalized().transpose() << " is not commensurate for periodicity" << std::endl;
                std::cout << "Mesh size in that direction must be a multiple of " << std::setprecision(15) << std::scientific << ld.cartesian().norm() << std::endl;
                std::cout << "Size detected=" << std::setprecision(15) << std::scientific << ccc.norm() << std::endl;
                std::cout << "Closest commensurate size=" << std::setprecision(15) << std::scientific << std::round(normRatio) * ld.cartesian().norm() << std::endl;
                assert(false && "MESH NOT COMMENSURATE");
            }
        }
        return B;
    }
    
    template<int dim>
    Eigen::Matrix<long int,dim,dim>  PeriodicGlidePlaneFactory<dim>::get_N(const Polycrystal<dim>& poly,const Eigen::Matrix<double,dim,dim>& B)
    {
        assert(poly.grains().size() == 1 && "Periodic simulations only supported in single crystals");
        const auto &grain(poly.grains().begin()->second);
        // const auto &meshRegion(grain.region);
        
        Eigen::Matrix<double,dim,dim> Nd(grain.latticeBasis.inverse()*B);
        Eigen::Matrix<double,dim,dim> Ni(Nd.array().round().matrix());
        if((Nd-Ni).norm()>FLT_EPSILON)
        {
            assert(0 && "Mesh not commensurate with the lattice");
        }
        return Ni.template cast<long int>();
    }
    template struct PeriodicGlidePlaneFactory<3>;
}
#endif
