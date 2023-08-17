/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2020 by Danny Perez <danny_perez@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDconfigFields_H_
#define model_DDconfigFields_H_

#include <map>

#include <Eigen/Dense>

#include <DislocationDynamicsBase.h>
#include <DDconfigIO.h>
#include <DislocationLoopPatches.h>
#include <PeriodicGlidePlaneFactory.h>
#include <EshelbyInclusionBase.h>
#include <SphericalInclusion.h>
#include <PolyhedronInclusion.h>



namespace model
{

    template <int dim>
    class DDconfigFields : private std::map<size_t,DislocationLoopPatches<dim>>
    /*                 */, private std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>
    /*                 */, public std::map<size_t,std::shared_ptr<EshelbyInclusionBase<dim>>>
    /*                 */, public std::map<size_t,PolyhedronInclusionNodeIO<dim>>
    {
        
    public:
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef std::map<size_t,std::shared_ptr<EshelbyInclusionBase<dim>>> EshelbyInclusionContainerType;
        typedef std::map<size_t,PolyhedronInclusionNodeIO<dim>> PolyhedronInclusionNodeContainerType;
        
        DislocationDynamicsBase<dim>& ddBase;
        const std::vector<VectorDim> periodicShifts;
        const DDconfigIO<dim>& configIO;
        
        DDconfigFields(DislocationDynamicsBase<dim>& ddBase_in,const DDconfigIO<dim>& configIO_in);
        void updateConfiguration();
        const std::map<size_t,DislocationLoopPatches<dim>>& loopPatches() const;
        std::map<size_t,DislocationLoopPatches<dim>>& loopPatches();
        const std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>& segments() const;
        std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim>>& segments();
        const PolyhedronInclusionNodeContainerType& polyhedronInclusionNodes() const;
        PolyhedronInclusionNodeContainerType& polyhedronInclusionNodes();
        const EshelbyInclusionContainerType& eshelbyInclusions() const;
        EshelbyInclusionContainerType& eshelbyInclusions();
        double solidAngle(const VectorDim& x) const;
        VectorDim dislocationPlasticDisplacement(const VectorDim& x) const;
        MatrixDim dislocationStress(const VectorDim& x) const;
        MatrixDim inclusionStress(const VectorDim& x) const;

    };

}
#endif
