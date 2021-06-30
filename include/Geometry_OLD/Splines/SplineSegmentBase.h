/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SplineSegmentBase_H_
#define model_SplineSegmentBase_H_

#include <math.h>
#include <Eigen/Dense>


namespace model
{
    
    /************************************************************************/
    /* SplineSegmentBase, general case **************************************/
    /************************************************************************/
    template <short unsigned int dim, short unsigned int corder>
    class SplineSegmentBase {};
    
    
    /************************************************************************/
    /* SplineSegmentBase, template specialization corder=0 ******************/
    /************************************************************************/
    template <short unsigned int dim>
    struct SplineSegmentBase<dim,0>
    {
        
        
        static constexpr int corder = 0;
        static constexpr int Ncoeff= 2*(corder+1);
        typedef Eigen::Matrix<double,Ncoeff,1>     VectorNcoeff;
        typedef Eigen::Matrix<double, 1, Ncoeff> RowNcoeff;
        typedef Eigen::Matrix<double, 1, Ncoeff-1> RowNcoeffu;
        typedef Eigen::Matrix<double, 1, Ncoeff-2> RowNcoeffuu;
        typedef Eigen::Matrix<double, Ncoeff, Ncoeff> MatrixNcoeff;
        typedef Eigen::Matrix<double, Ncoeff, dim> MatrixNcoeffDim;
        typedef Eigen::Matrix<double,dim,1>     VectorDim;
        typedef Eigen::Matrix<double,dim,dim>   MatrixDim;

        typedef std::map<size_t,
        /*            */ std::pair<VectorNcoeff,VectorDim>,
        /*            */ std::less<size_t>
        /*            */ > H2PmapType;
        
        /**********************************************************************/
        static RowNcoeff powers(const double& uin)
        {
            return (RowNcoeff()<<1.0, uin).finished();
        }
        
        /**********************************************************************/
        static RowNcoeffu powersDiff1(const double& )
        {
            return (RowNcoeffu()<<1.0).finished();
        }
        
        /**********************************************************************/
        static RowNcoeffuu powersDiff2(const double& )
        {
            return RowNcoeffuu::Zero();
        }
        
        /**********************************************************************/
        static MatrixNcoeff sfCoeffs(const double&)
        {/*! The matrix of shape function coefficients in Hermite form of this
          *  spline segment.
          */
            /*                         P0   P1 */
            return (MatrixNcoeff()<<  1.0, 0.0,             // u^0
                    /*            */ -1.0, 1.0).finished(); // u^1;
        }
        
        /**********************************************************************/
        static RowNcoeff sf(const double& uin,const double& g)
        {/*\param[in] uin parameter value in [0,1]
          *\returns The row vector of shape function at uin
          */
            return powers(uin)*sfCoeffs(g);
        }
        
        /******************************************************************************/
        static RowNcoeff sfDiff1(const double& uin,const double& g)
        {
            return powersDiff1(uin)*sfCoeffs(g).template block<Ncoeff-1,Ncoeff>(1,0);
        }
        
        /******************************************************************************/
        static RowNcoeff sfDiff2(const double& uin,const double& g)
        {
            return  powersDiff2(uin)*sfCoeffs(g).template block<Ncoeff-2,Ncoeff>(2,0);
        }

        
        /**********************************************************************/
        template<typename LinkType>
        static MatrixNcoeffDim hermiteDofs(const LinkType& link)
        {/*!\returns The matrix of Hermite dofs of this spline segment.
          *  [P0x P0y P0z;P1x P1y P1z]
          */
            return (MatrixNcoeffDim()<< link.sourceP.transpose(),
                    /*            */	link.  sinkP.transpose()).finished();
        }
        
        /**********************************************************************/
        template<typename LinkType>
        static VectorDim sourceT(const LinkType&)
        {
            return VectorDim::Zero();
        }
        
        /**********************************************************************/
        template<typename LinkType>
        static VectorDim sinkT(const LinkType&)
        {
            return VectorDim::Zero();
        }
        
//        /**********************************************************************/
//        template<typename LinkType>
//        static H2PmapType hermite2posMap(const LinkType& link) // move to CatmullRom?
//        {
//
//            //            std::cout<<link.source->sID<<"->"<<link.sink->sID<<std::endl;
//
//            //            std::map<size_t,VectorDim> posMap;
//            //            std::map<size_t,std::pair<double,double>> temp;
//
//            H2PmapType temp;
//            temp.emplace(link.source->snID(),std::make_pair((VectorNcoeff()<<1.0,0.0).finished(),link.sourceP));
//            temp.emplace(link.  sink->snID(),std::make_pair((VectorNcoeff()<<0.0,1.0).finished(),link.  sinkP));
//
//
//
//            return temp;
//        }
        
           /**********************************************************************/
            template<typename LinkType>
           static Eigen::Matrix<double,dim,Ncoeff> hermiteCoefficients(const LinkType& link)
           {
               Eigen::Matrix<double,dim,Ncoeff> HrCf;
               HrCf.col(0)= link.source->get_P();
               HrCf.col(1)= link.sink->get_P();
               return HrCf;
           }
    
        //        MatrixNcoeffDim get_qH() const {
        //            return (MatrixNcoeffDim()<< link.source->get_P().transpose(),
        //                    /*               */ link.  sink->get_P().transpose()).finished();
        //        }
        
        
        
        //        //////////////////////////////////////////////////////////////
        //        //BezierCoefficients: Bezier coefficients (uniform parametrization)
        //        Eigen::Matrix<double,dim+1,Ncoeff> BezierCoefficients() const {
        //            Eigen::Matrix<double,dim+1,Ncoeff> BzCf;
        //            double w =1.0;
        //            BzCf.col(0)<< link.source->get_P(), w;
        //            BzCf.col(2)<< link.sink  ->get_P(), w;
        //            return BzCf;
        //        }
        
        
    };
    
    
    
    /**************************************************************************/
    /* SplineSegmentBase, template specialization corder=1 ******************/
    /**************************************************************************/
    template <short unsigned int dim>
    struct SplineSegmentBase<dim,1>
    {
        
        static constexpr int corder = 1;
        static constexpr int Ncoeff= 2*(corder+1);
        typedef Eigen::Matrix<double, 1, Ncoeff> RowNcoeff;
        typedef Eigen::Matrix<double, 1, Ncoeff-1> RowNcoeffu;
        typedef Eigen::Matrix<double, 1, Ncoeff-2> RowNcoeffuu;
        typedef Eigen::Matrix<double, Ncoeff, Ncoeff> MatrixNcoeff;
        typedef Eigen::Matrix<double, Ncoeff, dim> MatrixNcoeffDim;
        typedef Eigen::Matrix<double,Ncoeff,1>     VectorNcoeff;
        typedef Eigen::Matrix<double,dim,1>     VectorDim;
        typedef Eigen::Matrix<double,dim,dim>   MatrixDim;
        
        typedef std::map<size_t,
        /*            */ std::pair<VectorNcoeff,VectorDim>,
        /*            */ std::less<size_t>> H2PmapType;
        
        /**********************************************************************/
        static RowNcoeff powers(const double & uin)
        {
            return (RowNcoeff()<<1.0, uin, std::pow(uin,2), std::pow(uin,3)).finished();
        }
        
        static RowNcoeffu powersDiff1(const double & uin)
        {
            return (RowNcoeffu()<<1.0, 2.0*uin, 3.0*std::pow(uin,2)).finished();
        }
        
        static RowNcoeffuu powersDiff2(const double & uin)
        {
            return (RowNcoeffuu()<<2.0, 6.0*uin).finished();
        }
        
        /**********************************************************************/
        static MatrixNcoeff sfCoeffs(const double& g)
        {/*!\returns The matrix of shape function coefficients in Hermite form of this
          *  spline segment.
          */
            /*                         P0      T0    P1   T1  */
            return (MatrixNcoeff()<<  1.0,    0.0,  0.0, 0.0,             // u^0
                    /*            */  0.0,      g,  0.0, 0.0,             // u^1
                    /*            */ -3.0, -2.0*g,  3.0,  -g,             // u^2
                    /*            */  2.0,      g, -2.0,   g).finished(); // u^3;
        }
        
        /******************************************************************************/
        static RowNcoeff sf(const double& uin,const double& g)
        {
            return powers(uin)*sfCoeffs(g);
        }
        
        /******************************************************************************/
        static RowNcoeff sfDiff1(const double& uin,const double& g)
        {
            return powersDiff1(uin)*sfCoeffs(g).template block<Ncoeff-1,Ncoeff>(1,0);
        }
        
        /******************************************************************************/
        static RowNcoeff sfDiff2(const double& uin,const double& g)
        {
            return  powersDiff2(uin)*sfCoeffs(g).template block<Ncoeff-2,Ncoeff>(2,0);
        }
        
        /**********************************************************************/
        template<typename LinkType>
        static MatrixNcoeffDim hermiteDofs(const LinkType& link)
        {/*!\returns The 4xdim matrix of Hermite dof of this spline segment.
          * [P0x P0y P0z;
          * T0x T0y T0z;
          * P1x P1y P1z;
          * T1x T1y T1z]
          */
            return (MatrixNcoeffDim()<< link.source->get_P().transpose(),
                    /*            */	sourceT(link).transpose(),
                    /*            */	link.  sink->get_P().transpose(),
                    /*            */	sinkT(link).transpose()).finished();
        }
        
        
        
        /**********************************************************************/
        template<typename LinkType>
        static VectorDim sourceT(const LinkType& link)
        {
            VectorDim temp=VectorDim::Zero();
            for(const auto& loopLink : link.loopLinks())
            {
                if(link.source->sID==loopLink->source()->sID && link.sink->sID==loopLink->sink()->sID)
                {
                    temp+=link.source->tangents().at(loopLink->loop()->sID);
                }
                else if(link.source->sID==loopLink->sink()->sID && link.sink->sID==loopLink->source()->sID)
                {
                    temp-=link.source->tangents().at(loopLink->loop()->sID);
                }
                else
                {
                    assert(0);
                }
            }
            
            return link.loopLinks().size()==0? temp : temp/link.loopLinks().size();
        }
        
        /**********************************************************************/
        template<typename LinkType>
        static VectorDim sinkT(const LinkType& link)
        {
            VectorDim temp=VectorDim::Zero();
            for(const auto& loopLink : link.loopLinks())
            {
                if(link.source->sID==loopLink->source()->sID && link.sink->sID==loopLink->sink()->sID)
                {
                    temp+=link.sink->tangents().at(loopLink->loop()->sID);
                }
                else if(link.source->sID==loopLink->sink()->sID && link.sink->sID==loopLink->source()->sID)
                {
                    temp-=link.sink->tangents().at(loopLink->loop()->sID);
                }
                else
                {
                    assert(0);
                }
            }
            
            return link.loopLinks().size()==0? temp : temp/link.loopLinks().size();
        }
        
        
        
        
        /**********************************************************************/
        template<typename LinkType>
        static std::pair<Eigen::Matrix<double,Ncoeff,Eigen::Dynamic>,Eigen::Matrix<double,Eigen::Dynamic,dim>> hermite2posMatrix(const LinkType& link)
        {
            H2PmapType temp=hermite2posMap(link);
            
            Eigen::Matrix<double,Ncoeff,Eigen::Dynamic> m(Ncoeff,temp.size());
            Eigen::Matrix<double,Eigen::Dynamic,dim>    n(temp.size(),dim);
            
            size_t k=0;
            for (const auto& ele : temp)
            {
                m.col(k)=ele.second.first;
                n.row(k)=ele.second.second.transpose();
                k++;
            }
            
            return std::make_pair(m,n);
        }
        
        /**********************************************************************/
        template<typename LinkType>
        static H2PmapType hermite2posMap(const LinkType& link) // move to CatmullRom?
        {
            
            //            std::cout<<link.source->sID<<"->"<<link.sink->sID<<std::endl;
            
            //            std::map<size_t,VectorDim> posMap;
            //            std::map<size_t,std::pair<double,double>> temp;
            
            H2PmapType temp;
            
            const std::map<size_t,std::map<size_t,std::pair<double,VectorDim>>> sourceMapMap=link.source->loopTangentCoeffs();
            const std::map<size_t,std::map<size_t,std::pair<double,VectorDim>>>   sinkMapMap=link.  sink->loopTangentCoeffs();
            
            
            for(const auto& loopLink : link.loopLinks())
            {
                
                temp.emplace(link.source->snID(),std::make_pair((VectorNcoeff()<<1.0,0.0,0.0,0.0).finished(),link.source->get_P()));
                temp.emplace(link.  sink->snID(),std::make_pair((VectorNcoeff()<<0.0,0.0,1.0,0.0).finished(),link.  sink->get_P()));
                
                const auto& sourceMap=sourceMapMap.at(loopLink->loop()->sID);
                const auto&   sinkMap=  sinkMapMap.at(loopLink->loop()->sID);
                
                for(const auto& pair : sourceMap)
                {
                    // initialize if not existent
                    const auto mapIter=temp.find(pair.first);
                    if(mapIter==temp.end())
                    {
                        temp.emplace(pair.first,std::make_pair(VectorNcoeff::Zero(),pair.second.second));
                    }
                    
                    if(link.source->sID==loopLink->source()->sID && link.sink->sID==loopLink->sink()->sID)
                    {
                        temp[pair.first].first(1)+=pair.second.first/link.loopLinks().size();
                    }
                    else if(link.source->sID==loopLink->sink()->sID && link.sink->sID==loopLink->source()->sID)
                    {
                        temp[pair.first].first(1)-=pair.second.first/link.loopLinks().size();
                    }
                    else
                    {
                        assert(0);
                    }
                }
                
                for(const auto& pair : sinkMap)
                {
                    
                    // initialize if not existent
                    const auto mapIter=temp.find(pair.first);
                    if(mapIter==temp.end())
                    {
                        temp.emplace(pair.first,std::make_pair(VectorNcoeff::Zero(),pair.second.second));
                    }
                    
                    if(link.source->sID==loopLink->source()->sID && link.sink->sID==loopLink->sink()->sID)
                    {
                        temp[pair.first].first(3)+=pair.second.first/link.loopLinks().size();
                    }
                    else if(link.source->sID==loopLink->sink()->sID && link.sink->sID==loopLink->source()->sID)
                    {
                        temp[pair.first].first(3)-=pair.second.first/link.loopLinks().size();
                    }
                    else
                    {
                        assert(0);
                    }
                }
                
            }
            
            return temp;
        }
        
        //
        //        /**********************************************************************/
        //        Eigen::Matrix<double,dim+1,Ncoeff> BezierCoefficients() const {
        //            /*!\returns The [dim x 4] matrix of Bezier coefficients of this spline segment
        //             * [P0 P0+T0/3 P1-T1/3 P1]
        //             */
        //            Eigen::Matrix<double,dim+1,Ncoeff> BzCf;
        //            double w =1.0;
        //            BzCf.col(0)<< link.source->get_P(), w;
        //            BzCf.col(1)<< link.source->get_P()+sourceT()/3.0*chordParametricLength(), w;
        //            BzCf.col(2)<< link.sink  ->get_P()-  sinkT()/3.0*chordParametricLength(), w;
        //            BzCf.col(3)<< link.sink  ->get_P(), w;
        //            return BzCf;
        //        }
        //
        //        /**********************************************************************/
        //        Eigen::Matrix<double,dim,Ncoeff> hermiteCoefficients() const
        //        {
        //            Eigen::Matrix<double,dim,Ncoeff> HrCf;
        //            HrCf.col(0)= link.source->get_P();
        //            HrCf.col(1)= sourceT()*chordParametricLength();
        //            HrCf.col(2)= link.sink->get_P();
        //            HrCf.col(3)= sinkT()*chordParametricLength();
        //            return HrCf;
        //        }
        //
        //        /************************************************************************/
        //        Eigen::Matrix<double,dim,Ncoeff> polynomialCoeff() const
        //        {/*!\returns The matrix of coefficients of the polynomial associated to this
        //          *  SplineSegmentBase. If C=polynomialCoeff() then the polynomial is:
        //          *  P(u)=C.col(0)+u*C.col(1)+u^2*C.col(2)+...
        //          */
        //            return Coeff2Hermite<pOrder>::template h2c<dim>(hermiteCoefficients());
        //        }
        
    };
    
    
    
    
    
    
    //	/************************************************************************/
    //	/* SplineSegmentBase, template specialization corder=2 ******************/
    //	/************************************************************************/
    //	template <typename Derived, short unsigned int dim>
    //	class SplineSegmentBase<Derived,dim,2> : public NetworkLink<Derived>,
    //	/*	                                        */ public ParametricCurve<Derived,dim> {
    //
    //
    //		enum {corder=2};
    //#include<model/Geometry/Splines/SplineEnums.h>
    //
    //
    //
    //        ////////////////////////////////////////////////////////////////
    //        //! Returns the length of the chord vector to the power alpha
    //        double chordParametricLength() const
    //        {
    //        	return std::pow(chordLength(),alpha);;
    //        }
    //
    //	public:
    //
    //        static double alpha;
    //
    //
    //		RowNcoeff powers(const double & uin){
    //			return (RowNcoeff()<<1.0, uin, std::pow(uin,2), std::pow(uin,3), std::pow(uin,4), std::pow(uin,5)).finished();
    //		}
    //
    //
    //		RowNcoeffu powersu(const double & uin){
    //			return (RowNcoeffu()<<1.0, 2.0*uin, 3.0*std::pow(uin,2), 4.0*std::pow(uin,3), 5.0*std::pow(uin,4)).finished();
    //		}
    //
    //
    //		RowNcoeffuu powersuu(const double & uin){
    //			return (RowNcoeffuu()<<2.0, 6.0*uin, 12.0*std::pow(uin,2), 20.0*std::pow(uin,3)).finished();
    //		}
    //
    //
    //		//////////////////////////////////////////////////////////////
    //		//get_SFCH
    //		MatrixNcoeff get_SFCH() const {
    //			//! !!!!!!!! FINISH HERE !!!!!!!! //
    //			assert(0);
    //			return MatrixNcoeff::Zero();
    //		}
    //
    //		MatrixNcoeffDim get_qH() const {
    //			return (MatrixNcoeffDim()<< link.source->get_P().transpose(),
    //					/*               */ sourceT().transpose(),
    //					/*               */ link.source->get_K().transpose(),
    //					/*               */ link.  sink->get_P().transpose(),
    //					/*               */ sinkT().transpose(),
    //					/*               */ link.  sink->get_K().transpose()).finished();
    //		}
    //
    //#include "SplineSegmentBase_common.h"
    //
    //	};
    //
    //    //static data
    //    template <typename Derived, short unsigned int dim>
    //	double SplineSegmentBase<Derived,dim,2>::alpha=0.5;
    ////	double SplineSegmentBase<Derived,dim,2>::alpha=1.0;
    //
    //
    //	//////////////////////////////////////////////////////////////s
    //} // namespace model
}
#endif

