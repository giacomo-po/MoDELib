/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicGlidePlaneFactory_H_
#define model_PeriodicGlidePlaneFactory_H_


#include <memory>
#include <string>
#include <list>
#include <WeakPtrFactories.h>

#include <GlidePlane.h>
#include <PeriodicGlidePlane.h>
#include <GlidePlaneFactory.h>
#include <TerminalColors.h>
#include <DiophantineSolver.h>
#include <PeriodicGlidePlane.h>
//#include <PeriodicGlidePlane.hpp>
namespace model
{
        
    template<int dim>
    struct PeriodicGlidePlaneFactory : public KeyConstructableWeakPtrFactory<PeriodicGlidePlaneFactory<dim>,PeriodicGlidePlane<dim>>
//    /*                              */,public Lattice<dim>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef PeriodicGlidePlaneFactory<dim> PeriodicGlidePlaneFactoryType;
        typedef PeriodicGlidePlane<dim> PeriodicGlidePlaneType;
        typedef KeyConstructableWeakPtrFactory<PeriodicGlidePlaneFactoryType,PeriodicGlidePlaneType> BaseType;
        typedef GlidePlane<dim> GlidePlaneType;
        typedef typename GlidePlaneType::KeyType GlidePlaneKeyType;
        typedef std::shared_ptr<PeriodicGlidePlaneType> PeriodicGlidePlaneSharedPtrType;
        
    public:
        
        const Polycrystal<dim>& poly;
        GlidePlaneFactory<dim>& glidePlaneFactory;
//        const Eigen::Matrix<long int,dim,dim> N; // B=AN, where A is lattice matrix
        PeriodicGlidePlaneFactory(const Polycrystal<dim>& poly_in,GlidePlaneFactory<dim>& glidePlaneFactory_in);
        // PeriodicGlidePlaneSharedPtrType get(const GlidePlaneType& plane);
        PeriodicGlidePlaneSharedPtrType get(const GlidePlaneKeyType& temp);
        BaseType& periodicGlidePlanes();
        const BaseType& periodicGlidePlanes() const;
//        static Eigen::Matrix<double,dim,dim> get_B(const Polycrystal<dim>& poly);
//        static Eigen::Matrix<long int,dim,dim> get_N(const Polycrystal<dim>& poly,const Eigen::Matrix<double,dim,dim>& B);
        
        
        //        GlidePlaneKeyType periodicPlaneKey(const GlidePlaneKeyType& key) const
        //        //        GlidePlaneKeyType periodicPlaneKey(const GlidePlaneType& plane) const
        //        {
        //
        //            const VectorDim meshCenter(0.5*(poly.mesh.xMax()+poly.mesh.xMin()));
        //            const ReciprocalLatticeVector<dim> r(key.reciprocalDirectionComponents(),poly.grains().begin()->second);
        //            const long int meshCenterPlaneIndex(r.closestPlaneIndexOfPoint(meshCenter));
        //            const auto gcd(LatticeGCD<dim>::gcd(key.reciprocalDirectionComponents().transpose()*N));
        //            const long int meshCenterPlaneIndexTemp=meshCenterPlaneIndex/gcd; // make mesh center belong to correct family by flooring
        //            const long int meshCenterPlaneIndexAct=meshCenterPlaneIndexTemp*gcd;
        //            const auto periodicPlaneIndex(key.planeIndex()>=0? key.planeIndex()%gcd : (key.planeIndex()%gcd)+gcd); // adjust sign
        //            return GlidePlaneKeyType(key.reciprocalDirectionComponents(),periodicPlaneIndex+meshCenterPlaneIndexAct,key.latticeID()); // get closest plane to center of mesh
        //        }
        
    };
    
}
#endif
