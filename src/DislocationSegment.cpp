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

#ifndef model_DislocationSegment_cpp_
#define model_DislocationSegment_cpp_

#include <DislocationSegment.h>

#ifndef NDEBUG
#define VerboseDislocationSegment(N,x) if(verboseDislocationSegment>=N){model::cout<<x;}
#else
#define VerboseDislocationSegment(N,x)
#endif

namespace model
{
    template <int dim, short unsigned int corder, typename InterpolationType>
    DislocationSegment<dim,corder,InterpolationType>::DislocationSegment(LoopNetworkType* const net,
                       const std::shared_ptr<NetworkNodeType>& nI,
                       const std::shared_ptr<NetworkNodeType>& nJ) :
    /* init */ NetworkLink<DislocationSegment>(net,nI,nJ)
    {
        //            std::cout<<"Creating DummyLoopNode "<<this->tag()<<std::endl;
    }
    
    
    /******************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationSegment<dim,corder,InterpolationType>::initFromFile(const std::string& fileName)
    {
//        LinkType::alpha=TextFileParser(fileName).readScalar<double>("parametrizationExponent",true);
//        assert((LinkType::alpha)>=0.0 && "parametrizationExponent MUST BE >= 0.0");
//        assert((LinkType::alpha)<=1.0 && "parametrizationExponent MUST BE <= 1.0");
        quadPerLength=TextFileParser(fileName).readScalar<double>("quadPerLength",true);
        //            assembleWithTangentProjection=TextFileParser(fileName).readScalar<int>("assembleWithTangentProjection",true);
        assert((NetworkLinkType::quadPerLength)>=0.0 && "quadPerLength MUST BE >= 0.0");
        verboseDislocationSegment=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseDislocationSegment",true);
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const typename DislocationSegment<dim,corder,InterpolationType>::MatrixDim DislocationSegment<dim,corder,InterpolationType>::I=typename DislocationSegment<dim,corder,InterpolationType>::MatrixDim::Identity();
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    const typename DislocationSegment<dim,corder,InterpolationType>::VectorDim DislocationSegment<dim,corder,InterpolationType>::zeroVector=typename DislocationSegment<dim,corder,InterpolationType>::VectorDim::Zero();
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    double DislocationSegment<dim,corder,InterpolationType>::quadPerLength=0.2;
    
//    template <int dim, short unsigned int corder, typename InterpolationType>
//    bool DislocationSegment<Derived>::assembleWithTangentProjection=false;
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    int DislocationSegment<dim,corder,InterpolationType>::verboseDislocationSegment=0;

}
#endif
