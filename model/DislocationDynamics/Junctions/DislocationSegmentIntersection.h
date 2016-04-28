/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Tamer Crosby <tamercrosby@gmail.com>,
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DISLOCATIONSEGMENTINTERSECTION_H_
#define model_DISLOCATIONSEGMENTINTERSECTION_H_

#include <float.h>
#include <assert.h>
#include <set>
#include <map>
#include <utility>
#include <math.h>
#include <float.h>
#include <Eigen/Dense>
#include <model/Geometry/Splines/Intersection/PlanarSplineImplicitization.h>
#include <model/Geometry/Splines/Coeff2Hermite.h>
#include <model/DislocationDynamics/DislocationLocalReference.h>


namespace model
{
    
    /******************************************************************************************************************/
    /* DislocationSegmentIntersection<3,3,1,1>: template specializatio (planar spline-spline intersection) ************/
    /******************************************************************************************************************/
    template <typename LinkType>
    class DislocationSegmentIntersection
    {
        
        enum {dim=3,polyDegree=3};
        enum {polyCoeff=polyDegree+1};
        typedef Eigen::Matrix<double,dim,polyCoeff> MatrixDimPolyCoeff;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
        
        
        /**********************************************************************/
        int planePlaneType(const VectorDim & normal1, const VectorDim & normal2, const VectorDim & pt1, const VectorDim & pt2) const
        {
            
            bool areParallelNormals(normal1.cross(normal2).norm()<FLT_EPSILON);
            bool areCoincidentPoints(std::fabs((pt1-pt2).dot(normal1))<FLT_EPSILON);
            
            int i;
            if(areParallelNormals && !areCoincidentPoints)
            {
                /*parallel planes: no intersection*/
                i= 0;
            }
            else if(areParallelNormals && areCoincidentPoints)
            {
                /*unique planes: planar intersection intersection*/
                i= 1;
            }
            else
            {
                /*angle planes: angular intersection*/
                i= 2;
            }
            
            return i;
        }
        
    public:
        //        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
        const LinkType& ds1;
        const MatrixDimPolyCoeff H1;
        const VectorDim& n1;
        const VectorDim P0;	  // end point of the spline
        const VectorDim T0;	  // tangent of the spline source
        const VectorDim P1;   // end point of the spline
        const VectorDim T1;   // tangent of the spline sink
        const VectorDim chord1;	  // end point of the spline
        //        const bool spline1Degen;
        //		const Eigen::Matrix<double,dim,dim> R1;
        
        static double compTol;
        static bool useNodeIntersection;

        /* Constructor *****************************************************/
        //        DislocationSegmentIntersection(const MatrixDimPolyCoeff& H1_in,
        //                                       /*                          */ const VectorDim& n1_in) : H1(H1_in),
        DislocationSegmentIntersection(const LinkType& ds, const VectorDim& n1_in) :
        /* init */ ds1(ds),
        /* init */ H1(ds1.hermiteCoefficients()),
        /* init */ n1(n1_in),
        /* init */ P0 (H1.col(0)),   // start point of the spline
        /* init */ T0 (H1.col(1)),   // tangent of the spline source
        /* init */ P1 (H1.col(2)),   // end point of the spline
        /* init */ T1 (H1.col(3)),   // tangent of the spline sink
        /* init */ chord1 (P1-P0)
        {   // chord vector
            const double chord1norm(chord1.norm());
            assert(chord1norm>FLT_EPSILON);
            //            assert(std::fabs(n1.dot(chord1/chord1norm))<FLT_EPSILON);
            //            assert(std::fabs(n1.dot(H1.col(1)))<FLT_EPSILON);
            //            assert(std::fabs(n1.dot(H1.col(3)))<FLT_EPSILON);
        }
        
        
        /**********************************************************************/
        std::set<std::pair<double,double> > intersectWith(const LinkType& ds2,
                                                          const VectorDim& n2,
                                                          const double& physTol,
                                                          const double& avoidNodeIntersection) const
        {
            const double OneMinusAvoidNodeIntersection=1.0-avoidNodeIntersection;
            
            const MatrixDimPolyCoeff H2=ds2.hermiteCoefficients();
            
            std::set<std::pair<double,double> > intersectionParameters;
            
            
            
            // chech that H1 is planar
            const VectorDim chord2=H2.col(2)-H2.col(0);
            const double chord2norm(chord2.norm());
            assert(chord2norm>FLT_EPSILON);
            
            const double absN2dotC2(std::fabs(n2.dot(chord2/chord2norm)));
            const double absN2dotT2source(std::fabs(n2.dot(H2.col(1))));
            const double   absN2dotT2sink(std::fabs(n2.dot(H2.col(3))));
            if (absN2dotC2>FLT_EPSILON || absN2dotT2source>FLT_EPSILON || absN2dotT2sink>FLT_EPSILON)
            {
                model::cout<<"n2="<<n2.transpose()<<std::endl;
                model::cout<<"H2="<<H2<<std::endl;
                model::cout<<"absN2dotC2="<<absN2dotC2<<"\n";
                model::cout<<"absN2dotT2source="<<absN2dotT2source<<"\n";
                model::cout<<"absN2dotT2sink="<<absN2dotT2sink<<"\n";
            }
            
            
            const VectorDim P2 = H2.col(0);	  // end point of the spline
            const VectorDim T2 = H2.col(1);	  // tangent of the spline source
            const VectorDim P3 = H2.col(2);   // end point of the spline
            const VectorDim T3 = H2.col(3);   // tangent of the spline sink
            
            if((0.5*(P0+P1)-0.5*(P2+P3)).norm()<0.75*(chord1.norm()+chord2.norm()))
            {
                
                if(useNodeIntersection)
                {
                    
                    const double nodeTol=physTol*0.25;
                    
                    if(ds1.source->sID==ds2.source->sID )
                    {// intersectionIsSourceSource
                        std::pair<double,std::pair<double,VectorDim> > map1=ds1.closestPoint(ds2.sink->get_P());
                        std::pair<double,std::pair<double,VectorDim> > map2=ds2.closestPoint(ds1.sink->get_P());
                        if(map1.first<=nodeTol && map1.second.first>avoidNodeIntersection && map2.first>nodeTol)
                        {
                            //std::cout<<"DislocationSegmentIntersection case 1 "<<map1.second.first<<std::endl;
                            intersectionParameters.emplace(map1.second.first,1.0);
                        }
                        else if(map1.first>nodeTol && map2.first<=nodeTol && map2.second.first>avoidNodeIntersection)
                        {
                            //std::cout<<"DislocationSegmentIntersection case 2 "<<map2.second.first<<std::endl;
                            intersectionParameters.emplace(1.0,map2.second.first);
                        }
                        else if(map1.first<=nodeTol && map2.first<=nodeTol)
                        {
                            //std::cout<<"DislocationSegmentIntersection case 3 "<<std::endl;
                            intersectionParameters.emplace(1.0,1.0);
                        }
                        else
                        {
                            // do nothing
                        }
                    }
                    else if(ds1.source->sID==ds2.sink->sID )
                    {// intersectionIsSourceSink
                        std::pair<double,std::pair<double,VectorDim> > map1=ds1.closestPoint(ds2.source->get_P());
                        std::pair<double,std::pair<double,VectorDim> > map2=ds2.closestPoint(ds1.sink->get_P());
                        if(map1.first<=nodeTol && map1.second.first>avoidNodeIntersection && map2.first>nodeTol)
                        {
                            //std::cout<<"DislocationSegmentIntersection case 4 "<<map1.second.first<<std::endl;
                            intersectionParameters.emplace(map1.second.first,0.0);
                        }
                        else if(map1.first>nodeTol && map2.first<=nodeTol && map2.second.first<OneMinusAvoidNodeIntersection)
                        {
                            //std::cout<<"DislocationSegmentIntersection case 5 "<<map2.second.first<<std::endl;
                            intersectionParameters.emplace(1.0,map2.second.first);
                        }
                        else if(map1.first<=nodeTol && map2.first<=nodeTol)
                        {
                            //std::cout<<"DislocationSegmentIntersection case 6 "<<std::endl;
                            intersectionParameters.emplace(1.0,0.0);
                        }
                        else
                        {
                            // do nothing
                        }
                    }
                    else if(ds1.sink->sID==ds2.source->sID )
                    {// intersectionIsSinkSource
                        std::pair<double,std::pair<double,VectorDim> > map1=ds1.closestPoint(ds2.sink->get_P());
                        std::pair<double,std::pair<double,VectorDim> > map2=ds2.closestPoint(ds1.source->get_P());
                        if(map1.first<=nodeTol && map1.second.first<OneMinusAvoidNodeIntersection && map2.first>nodeTol)
                        {
                            //std::cout<<"DislocationSegmentIntersection case 7 "<<map1.second.first<<std::endl;
                            intersectionParameters.emplace(map1.second.first,1.0);
                        }
                        else if(map1.first>nodeTol && map2.first<=nodeTol && map2.second.first>avoidNodeIntersection)
                        {
                            //std::cout<<"DislocationSegmentIntersection case 8 "<<map2.second.first<<std::endl;
                            intersectionParameters.emplace(0.0,map2.second.first);
                        }
                        else if(map1.first<=nodeTol && map2.first<=nodeTol)
                        {
                            //std::cout<<"DislocationSegmentIntersection case 9 "<<std::endl;
                            intersectionParameters.emplace(0.0,1.0);
                        }
                        else
                        {
                            // do nothing
                        }
                    }
                    else if(ds1.sink->sID==ds2.sink->sID )
                    {//intersectionIsSinkSink
                        std::pair<double,std::pair<double,VectorDim> > map1=ds1.closestPoint(ds2.source->get_P());
                        std::pair<double,std::pair<double,VectorDim> > map2=ds2.closestPoint(ds1.source->get_P());
                        if(map1.first<=nodeTol && map1.second.first<OneMinusAvoidNodeIntersection && map2.first>nodeTol)
                        {
                            //std::cout<<"DislocationSegmentIntersection case 10 "<<map1.second.first<<std::endl;
                            intersectionParameters.emplace(map1.second.first,0.0);
                        }
                        else if(map1.first>nodeTol && map2.first<=nodeTol && map2.second.first<OneMinusAvoidNodeIntersection)
                        {
                            //std::cout<<"DislocationSegmentIntersection case 11 "<<map2.second.first<<std::endl;
                            intersectionParameters.emplace(0.0,map2.second.first);
                        }
                        else if(map1.first<=nodeTol && map2.first<=nodeTol)
                        {
                            //std::cout<<"DislocationSegmentIntersection case 12 "<<map1.second.first<<std::endl;
                            intersectionParameters.emplace(0.0,0.0);
                        }
                        else
                        {
                            // do nothing
                        }
                    }
                    else
                    {
                    // no node intersection
                    }
                }
                
                

                if(intersectionParameters.empty()) // sements have no node intersections
                {
                    const int planesType=planePlaneType(n1,n2,P0,P2);
                    switch (planesType)
                    {
                        case 1: // coplanar planes
                        {
                            //                        											////std::cout<<"Coplanar case ";
                            const Eigen::Matrix<double,dim,dim> R1(DislocationLocalReference<dim>::global2local(chord1,n1));
                            
                            Eigen::Matrix<double,dim-1,polyCoeff> H1L;
                            H1L.col(0)=(R1*(P0-P0)).template segment<dim-1>(0);
                            H1L.col(1)=(R1*T0).template segment<dim-1>(0);
                            H1L.col(2)=(R1*(P1-P0)).template segment<dim-1>(0);
                            H1L.col(3)=(R1*T1).template segment<dim-1>(0);
                            
                            Eigen::Matrix<double,dim-1,polyCoeff> H2L;
                            H2L.col(0)=(R1*(P2-P0)).template segment<dim-1>(0);
                            H2L.col(1)=(R1*T2).template segment<dim-1>(0);
                            H2L.col(2)=(R1*(P3-P0)).template segment<dim-1>(0);
                            H2L.col(3)=(R1*T3).template segment<dim-1>(0);
                            
                            PlanarSplineImplicitization<polyDegree> sli(Coeff2Hermite<polyDegree>::template h2c<dim-1>(H1L));
                            PlanarSplineImplicitization<polyDegree>::physTol=physTol;
                            //						intersectionParameters=sli.template intersectWith<polyDegree>(Coeff2Hermite<polyDegree>::template h2c<dim-1>(H2L),physTol);
                            intersectionParameters=sli.template intersectWith<polyDegree>(Coeff2Hermite<polyDegree>::template h2c<dim-1>(H2L));
                            
                            break;
                        }
                            
                        case 2: // incident planes
                        {
                            //                       						////std::cout<<"Incident case ";
                            std::set<std::pair<double,double> > lineIntersectionParameters1;
                            std::set<std::pair<double,double> > lineIntersectionParameters2;
                            
                            // finding the common line
                            const double denom(1.0-std::pow(n1.dot(n2),2));
                            const double numer((P2-P0).dot(n2));
                            
                            if(std::fabs(denom)>compTol)
                            { // planes are incident
                                const double u=numer/denom;
                                const VectorDim linePoint = P0+(n2-n2.dot(n1)*n1)*u;
                                const VectorDim lineDir   = n1.cross(n2).normalized();
                                
                                Eigen::Matrix<double,4,1> projPoints;
                                projPoints<<(P0-linePoint).dot(lineDir),
                                /**********/(P1-linePoint).dot(lineDir),
                                /**********/(P2-linePoint).dot(lineDir),
                                /**********/(P3-linePoint).dot(lineDir);
                                
                                const Eigen::Matrix<double,dim,1> Pmean=linePoint+0.5*(projPoints.minCoeff()+projPoints.maxCoeff())*lineDir;
                                const Eigen::Matrix<double,dim,1> P4=Pmean+(chord1.norm()+chord2.norm())*lineDir;
                                const Eigen::Matrix<double,dim,1> P5=Pmean-(chord1.norm()+chord2.norm())*lineDir;
                                
                                const Eigen::Matrix<double,dim,dim> R1(DislocationLocalReference<dim>::global2local(P5-P4,n1));
                                const Eigen::Matrix<double,dim,dim> R2(DislocationLocalReference<dim>::global2local(P5-P4,n2));
                                
                                
                                // Intersect first spine and common line
                                Eigen::Matrix<double,dim-1,polyCoeff> H1L;
                                H1L.col(0)=(R1*(P0-P4)).template segment<dim-1>(0);
                                H1L.col(1)=(R1*T0).template segment<dim-1>(0);
                                H1L.col(2)=(R1*(P1-P4)).template segment<dim-1>(0);
                                H1L.col(3)=(R1*T1).template segment<dim-1>(0);
                                Eigen::Matrix<double,dim-1,2> H2L;
                                H2L.col(0)=(R1*(P4-P4)).template segment<dim-1>(0);
                                H2L.col(1)=(R1*(P5-P4)).template segment<dim-1>(0);
                                PlanarSplineImplicitization<polyDegree> sli1(Coeff2Hermite<polyDegree>::template h2c<dim-1>(H1L));
                                PlanarSplineImplicitization<polyDegree>::physTol=physTol;
                                lineIntersectionParameters1=sli1.template intersectWith<1>(Coeff2Hermite<1>::template h2c<dim-1>(H2L));
                                
                                
                                // Intersect second spline and common line
                                Eigen::Matrix<double,dim-1,polyCoeff> H3L;
                                H3L.col(0)=(R2*(P2-P4)).template segment<dim-1>(0);
                                H3L.col(1)=(R2*T2).template segment<dim-1>(0);
                                H3L.col(2)=(R2*(P3-P4)).template segment<dim-1>(0);
                                H3L.col(3)=(R2*T3).template segment<dim-1>(0);
                                H2L.col(0)=(R2*(P4-P4)).template segment<dim-1>(0);
                                H2L.col(1)=(R2*(P5-P4)).template segment<dim-1>(0);
                                PlanarSplineImplicitization<polyDegree> sli2(Coeff2Hermite<polyDegree>::template h2c<dim-1>(H3L));
                                PlanarSplineImplicitization<polyDegree>::physTol=physTol;
                                lineIntersectionParameters2=sli2.template intersectWith<1>(Coeff2Hermite<1>::template h2c<dim-1>(H2L));
                                
                                
                                
                                // du = dl / j = dl/L for a line
                                double tolU=physTol/(P5-P4).norm();
                                
                                
                                for (std::set<std::pair<double,double> >::const_iterator iter1 =lineIntersectionParameters1.begin();
                                     /*                                               */ iter1!=lineIntersectionParameters1.end();)
                                {
                                    
                                    // find the closest to iter1
                                    
                                    std::map<double,std::set<std::pair<double,double> >::const_iterator> compareTo1;
                                    
                                    for (std::set<std::pair<double,double> >::const_iterator iter2 =lineIntersectionParameters2.begin();
                                         /*                                               */ iter2!=lineIntersectionParameters2.end();
                                         /*                                               */ iter2++)
                                    {
                                        // sort the difference in the intersections on the line
                                        compareTo1.insert(std::make_pair(std::fabs(iter1->second-iter2->second),iter2));
                                    }
                                    
                                    if(compareTo1.size())
                                    {
                                        if(compareTo1.begin()->first<tolU)
                                        {
                                            double u1=iter1->first;
                                            double u2=compareTo1.begin()->second->first;
                                            intersectionParameters.insert(std::make_pair(u1,u2));
                                            std::set<std::pair<double,double> >::const_iterator toBeErased(iter1);
                                            ++iter1;
                                            assert(lineIntersectionParameters1.erase(*toBeErased) && "NONE ERASED");
                                            assert(lineIntersectionParameters2.erase(*(compareTo1.begin()->second)) && "NONE DELETED");
                                        }
                                        else
                                        {
                                            ++iter1;
                                        }
                                    }
                                    else{
                                        ++iter1;
                                    }
                                }
                            }
                            else
                            { // planes are parallel
                                assert(0 && "SOMETHING WENT REALLY, REALLY WRONG HERE, YOU SHOULD HAVE FOUND THIS ABOVE");
                            }
                            
                            
                            break;
                        }
                            
                        default:
                            
                            break;
                    }
                }
                
                
                
                
            }
            return intersectionParameters;
        } // close constructor
        
    };
    
    
    
    // Declare statica data members
    template <typename LinkType>
    double DislocationSegmentIntersection<LinkType>::compTol=FLT_EPSILON;
 
    template <typename LinkType>
    bool DislocationSegmentIntersection<LinkType>::useNodeIntersection=true;
    
} // namespace model
#endif

