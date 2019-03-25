/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po       <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby     <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel         <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed  <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONSEGMENT_H
#define model_DISLOCATIONSEGMENT_H

#include <PlanarDislocationSegment.h>

#ifndef NDEBUG
#define VerboseDislocationSegment(N,x) if(verboseDislocationSegment>=N){model::cout<<x;}
#else
#define VerboseDislocationSegment(N,x)
#endif

namespace model
{
    template <int dim, short unsigned int corder, typename InterpolationType>
    struct DislocationSegment : public PlanarDislocationSegment<DislocationSegment<dim,corder,InterpolationType>>
    {
        typedef typename TypeTraits<DislocationSegment<dim,corder,InterpolationType>>::NodeType NodeType;
        typedef PlanarDislocationSegment<DislocationSegment<dim,corder,InterpolationType>> PlanarDislocationSegmentType;
        static constexpr int Ncoeff=PlanarDislocationSegmentType::Ncoeff;
        static constexpr int pOrder=PlanarDislocationSegmentType::pOrder;
        typedef typename PlanarDislocationSegmentType::VectorDim VectorDim;
        typedef typename PlanarDislocationSegmentType::MatrixDim MatrixDim;
        typedef DislocationParticle<dim> DislocationParticleType;

        /**********************************************************************/
        DislocationSegment(const std::shared_ptr<NodeType>& nI,
                           const std::shared_ptr<NodeType>& nJ) :
        /* init */ PlanarDislocationSegmentType(nI,nJ)
        {

        }

        /**********************************************************************/
        std::pair<VectorDim,VectorDim> vacancyConcentrationMatrices(const VectorDim& x) const
        {
            const double a(this->chordLengthSquared());
            const double b(-2.0*(x-this->source->get_P()).dot(this->chord()));
            const double c((x-this->source->get_P()).squaredNorm()+DislocationStress<dim>::a2);
            const double ba(b/a);
            const double ca(c/a);
            const double sqbca(sqrt(1.0+ba+ca));
            const double sqca(sqrt(ca));
            const double logTerm(log((2.0*sqbca+2.0+ba)/(2.0*sqca+ba)));
            const double M0(-1.0/(4.0*M_PI*this->network().poly.Dv)/this->chordLength()*((1.0+0.5*ba)*logTerm-sqbca+sqca));
            const double M1(-1.0/(4.0*M_PI*this->network().poly.Dv)/this->chordLength()*(-0.5*ba*logTerm+sqbca-sqca));
            
            return std::make_pair(M0*this->chord().cross(this->burgers()),M1*this->chord().cross(this->burgers()));
        }
        
        /**********************************************************************/
        double vacancyConcentration(const VectorDim& x) const
        {
//            const double a(this->chordLengthSquared());
//            const double b(-2.0*(x-this->source->get_P()).dot(this->chord()));
//            const double c((x-this->source->get_P()).squaredNorm()+DislocationStress<dim>::a2);
//            const double ba(b/a);
//            const double ca(c/a);
//            const double sqbca(sqrt(1.0+ba+ca));
//            const double sqca(sqrt(ca));
//            const double logTerm(log((2.0*sqbca+2.0+ba)/(2.0*sqca+ba)));
//            const double M0((1.0+0.5*ba)*logTerm-sqbca+sqca);
//            const double M1(-0.5*ba*logTerm+sqbca-sqca);
//            return -1.0/(4.0*M_PI*this->network().poly.Dv)/this->chordLength()*(this->chord().cross(this->burgers()).dot(M0*this->source->climbVelocity()+M1*this->sink->climbVelocity()));
            const std::pair<VectorDim,VectorDim> vCM(vacancyConcentrationMatrices(x));
            return vCM.first.dot(this->source->climbVelocity())+vCM.second.dot(this->sink->climbVelocity());
        }
        
        /**********************************************************************/
        void addToVacancyConcentrationAssembly(const VectorDim& x,
                                               std::vector<Eigen::Triplet<double>>& lhsT,
                                               size_t& globalRow)
        {
            
            const size_t gSource(dim*this->source->gID());
            const size_t gSink(dim*this->sink->gID());
            
                // LHS
                const std::pair<VectorDim,VectorDim> vCM(vacancyConcentrationMatrices(x));
            
            
                for(int k=0;k<dim;++k)
                {
                    lhsT.emplace_back(globalRow,gSource+k,vCM.first(k));
                    lhsT.emplace_back(globalRow,gSink+k,vCM.second(k));
                }
                
                // RHS
            

    
        }
        
        /**********************************************************************/
        const MatrixDim& midPointStress() const __attribute__ ((deprecated))
        {/*!\returns The stress matrix for the centre point over this segment.*/
            return this->quadraturePoints().size()? quadraturePoint(this->quadraturePoints().size()/2).stress : MatrixDim::Zero();            
        }
    };
}
#endif
