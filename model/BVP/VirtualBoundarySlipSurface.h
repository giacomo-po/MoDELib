/* This file is part of model, the Mechanics of Defects Evolution Library.
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_VIRTUALBOUNDARYSLIPSURFACE_H_
#define  model_VIRTUALBOUNDARYSLIPSURFACE_H_

#include <boost/ptr_container/ptr_vector.hpp>
#include <Eigen/Dense>
#include <stdio.h>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/DislocationDynamics/Materials/CrystalOrientation.h>


namespace model {
    
    template <short unsigned int dim>
    struct RadialVirtualSegment {
        
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
        const VectorDim sourceP;
        const VectorDim sinkP;
        const MatrixDim Qt;
        const VectorDim Burg;
        
        RadialVirtualSegment (const VectorDim sourceP_in, const VectorDim sinkP_in, const MatrixDim Qt_in, const VectorDim Burg_in) :
        sourceP(sourceP_in), sinkP(sinkP_in), Qt(Qt_in), Burg(Burg_in) {}
        
    };
    
    //   template <>
    //    template <typename DislocationSegmentType, short unsigned int dim>
    template <typename DislocationSegmentType>
    class VirtualBoundarySlipSurface {
        
        enum{dim=TypeTraits<DislocationSegmentType>::dim};
        
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
        
        
        typedef typename GlidePlaneObserver<DislocationSegmentType>::GlidePlaneSharedPtrType GlidePlaneSharedPtrType;
    	typedef typename CrystalOrientation<dim>::PlaneNormalContainerType PlaneNormalContainerType;
        
    public:
        const GlidePlaneSharedPtrType pGlidePlane;
        const static double L;
        const VectorDim sourceP;
        const VectorDim sinkP;
        const VectorDim sourceN;
        const VectorDim sinkN;
        const VectorDim sourcePatL;
        const VectorDim   sinkPatL;
        const VectorDim Burgers;
        const double coreL;
        
        static double mu, nu ;
        
        const double gp_H;
        const VectorDim  gp_N;         // point on the glide plane, and the glide plane normal
        
        MatrixDim Qt1 , Qt2 , Qt3;           // orthogonal transformation matrices for the 3 straight segments
        static double C1, C2, C3, C4;
        
        boost::ptr_vector< RadialVirtualSegment<dim> > radialSegmentsVector;
        
        
        
        /* Constructor ********************************************************/
        VirtualBoundarySlipSurface(const DislocationSegmentType& ds) :
        /* init list                                     */ pGlidePlane(ds.pGlidePlane),
        /* init list                                     */ sourceP(ds.source->get_P()),
        /* init list                                     */   sinkP(ds.sink->get_P()),
        /* init list                                     */ sourceN(getOutDirection(ds,*(ds.source))),
        /* init list                                     */   sinkN(getOutDirection(ds,*(ds.sink))),
        /* init list                                     */ sourcePatL(sourceP+L*sourceN),
        /* init list                                     */   sinkPatL(  sinkP+L*sinkN),
        /* init list                                     */ Burgers(ds.Burgers),
        /* init list                                     */ coreL(std::pow(DislocationSegmentType::coreLsquared,0.5)),
        /* init list                                     */ gp_H(ds.pGlidePlane->height),
        /* init list                                     */ gp_N(ds.pGlidePlane->planeNormal){
            
            
            setElasticConstants();
            pGlidePlane->addToGLidePlane(this);
            //std::cout << sourceN.transpose() << "    :    " << sinkN.transpose() << std::endl;
            //----------calculate the displacement field induced by those segments on the boundary nodes ----
            for (unsigned int dN=0; dN<ds.shared.domain.nodeContainer.size();++dN){
                if(ds.shared.domain.nodeContainer[dN].isBoundaryNode) {
                    if(ds.shared.domain.nodeContainer[dN].triIDs.size() == 0) assert(0&&"Error: Boundary node without triangle element index array");
                    ds.shared.domain.nodeContainer[dN].uVir+= displacement(ds.shared.domain.nodeContainer[dN].P , ds.shared.domain.triContainer[ds.shared.domain.nodeContainer[dN].triIDs[0]]->outNormal);
                }
            }
        }
        
        
        /* Destructor *************************************************/
        ~VirtualBoundarySlipSurface(){
            pGlidePlane->removeFromGlidePlane(this);
        }
        
        
        //=========================================================================
        // function to set the elastic constants
        //========================================================================
        // template <typename DislocationSegmentType>
        void setElasticConstants () {
            //            mu=ds.shared.material.mu;
            //            nu=ds.shared.material.nu;
            mu =Material<Isotropic>::mu;
			nu =Material<Isotropic>::nu;
            
            C1=1.0-nu;
            C2=mu/(4.0*M_PI*C1);
            C3=1.0-2.0*nu;
            C4=1.0/(8.0*M_PI*C1);
        }
        
        //        //=========================================================================
        //        // function to set the elastic constants
        //        //========================================================================
        //        template <typename SharedType>
        //        void setElasticConstants (const SharedType* sharedPtr) {
        //            mu=sharedPtr->material.mu;
        //            nu=sharedPtr->material.nu;
        //
        //            C1=1.0-nu;
        //            C2=mu/(4.0*M_PI*C1);
        //            C3=1.0-2.0*nu;
        //            C4=1.0/(8.0*M_PI*C1);
        //        }
        
        //==========================================================================
        // function to add radial segments to the radialSegmentsVector
        // index =2: add at both source and sink
        // index =0: add only at source
        // index =1: add only at sink
        //=========================================================================
        void addRadialSegments (const unsigned int& index = 2) {
            
            if (index==0 || index==2){
                std::auto_ptr<RadialVirtualSegment<dim> > pRVS1 (new RadialVirtualSegment<dim> (sourceP,sourcePatL, getOrthogonalMatrix(sourceP,sourcePatL),Burgers) );
                radialSegmentsVector.push_back(pRVS1);
            }
            
            if (index==1 || index==2){
                std::auto_ptr<RadialVirtualSegment<dim> > pRVS2 (new RadialVirtualSegment<dim> (sinkP,sinkPatL, getOrthogonalMatrix(sinkP,sinkPatL), -Burgers) );
                radialSegmentsVector.push_back(pRVS2);
            }
        }
        
        //==========================================================================
        // function to find the direction in which virtual dislocation segments will be displaced outside the domain
        //==========================================================================
        // template <typename DislocationSegmentType , typename DislocationNodeType>
        template <typename DislocationNodeType>
        VectorDim getOutDirection(const DislocationSegmentType& ds, const DislocationNodeType& dn){
            
            //----------- displace the node in the direction of boundary normal, then project this to the segments glide plane -----
            VectorDim P_BN = dn.get_P() + 1000.0*dn.boundaryNormal;
            
            //-------------- project this point to the glide plane ------------
            VectorDim P_Prj = P_BN + ( ( (ds.pGlidePlane->height*ds.pGlidePlane->planeNormal)-P_BN ).dot(ds.pGlidePlane->planeNormal) )*ds.pGlidePlane->planeNormal;
            
            VectorDim temp = VectorDim::Zero();
            bool calculated = false;
	    
	    //VectorDim conjugatePlaneNormal = ds.conjugatePlaneNormal(ds.chord());
            //VectorDim conjugatePlaneNormal = ds.get_sessileNormal(ds.chord(),ds.Burgers);
	    PlaneNormalContainerType allowedSlipSystems;
	    CrystalOrientation<dim>::find_slipSystem(ds.chord(),ds.Burgers,allowedSlipSystems);
	    
	    VectorDim conjugatePlaneNormal;
	    
            if((P_Prj - dn.get_P()).norm() > 1.0e-7) {
                temp =  (P_Prj - dn.get_P()).normalized() ;
                temp = dn.boundaryNormal.dot(temp) / std::abs(dn.boundaryNormal.dot(temp)) * temp;
                calculated = true;
            }
            
            //else if (conjugatePlaneNormal.norm() > 0.0) {
            else if (allowedSlipSystems.size() > 1) {
	      
	      for (unsigned int k=0;k<allowedSlipSystems.size();k++) {
		if ((ds.pGlidePlane->planeNormal.cross(allowedSlipSystems[k])).norm() < FLT_EPSILON) continue;
		conjugatePlaneNormal = allowedSlipSystems[k].normalized();
	      }
	      
	      
	      double h = dn.get_P().dot(conjugatePlaneNormal);
                P_Prj = P_BN + ( ( (h*conjugatePlaneNormal)-P_BN ).dot(conjugatePlaneNormal) )*conjugatePlaneNormal;
                
                if((P_Prj - dn.get_P()).norm() > 1.0e-7) {
                    temp =  (P_Prj - dn.get_P()).normalized() ;
                    temp = dn.boundaryNormal.dot(temp) / std::abs(dn.boundaryNormal.dot(temp)) * temp;
                    calculated = true;
                }
            }
            
            if (!calculated) {
                std::cout << dn.get_P().transpose() << std::endl;
                std::cout << ds.pGlidePlane->planeNormal.transpose() << std::endl;
                std::cout << conjugatePlaneNormal.transpose() << std::endl;
                assert(0&& "Radial projection direction for virtual dislocations could not be calculated");
            }
            
            
            return temp;
        }
        
        //==================================================================================================
        // function to calculate the coordinate orthogonal transformation matrix for a given straight dislocation segment
        //==================================================================================================
        
        MatrixDim getOrthogonalMatrix (const VectorDim& source, const VectorDim& sink ) const {
            
            VectorDim linUV = (sink-source).normalized();
            
            MatrixDim TM;
            
            double cosTheta, theta;
            
            //-- the two angles, phi and theta , defining the dislocation line direction
            
            double phi = acos(linUV(2));
            
            if   (sin(phi)!=0.0) 	cosTheta = linUV(0) / sin(phi);           // = cos (theta)
            else 	                cosTheta = 1.0e0;
            
            if (cosTheta >  1.0) cosTheta = 1.0e0;
            if (cosTheta < -1.0) cosTheta = -1.0e0;
            
            if ( linUV(1) >= 0.0 ) theta = acos (cosTheta);
            else                   theta = (2.0e0*M_PI)-acos(cosTheta);
            
            //--- the orthogonal transformation matrix elements (this is actually transpose(Qy)*transpose(Qz))---
            //----------------, which is the transpose of the total rotation matrix. ----------------------------
            //--- So to rotate any vector to the new coordinate multiply it directly by this matrix ------------
            
            
            TM(0,0) =  cos(phi) * cos(theta);
            TM(0,1) =  cos(phi) * sin(theta);
            TM(0,2) = -sin(phi);
            TM(1,0) = -sin(theta);
            TM(1,1) =  cos(theta);
            TM(1,2) =  0.0e0;
            TM(2,0) =  sin(phi) * cos(theta);
            TM(2,1) =  sin(phi) * sin(theta);
            TM(2,2) =  cos(phi);
            
            return TM;
        }
        
        
        
        /* stress ****************************************************/
        MatrixDim stress(const VectorDim& Rfield) const{
            //return ( stress_CoorDpndt(Rfield,sourceP,sourcePatL,Qt1) + stress_CoorDpndt(Rfield,sourcePatL,sinkPatL,Qt2) + stress_CoorDpndt(Rfield,sinkPatL,sinkP,Qt3) );
            //return ( stress_CoorDpndt(Rfield,sourceP,sourcePatL,Qt1) + stress_CoorDpndt(Rfield,sinkPatL,sinkP,Qt3) );
            MatrixDim temp = MatrixDim::Zero();
            
            for (unsigned int i = 0; i<radialSegmentsVector.size(); i++) temp+= stress_CoorDpndt(Rfield,radialSegmentsVector[i].sourceP,radialSegmentsVector[i].sinkP,
                                                                                                 radialSegmentsVector[i].Qt, radialSegmentsVector[i].Burg);
            
            return temp;
        }
        
        
        //==================================================================================================
        // stress field at point R induced by straight segment x1 -> x2 using the coordinate dependent forms
        //==================================================================================================
        
        MatrixDim stress_CoorDpndt(const VectorDim & R, const VectorDim & x1 , const VectorDim & x2 , const MatrixDim & Qt) const {
            
            VectorDim B = Qt * Burgers;
            VectorDim Rm = Qt*(R-x1);
            
            double l1 = -Rm(2);
            double l2 = (x2-x1).norm()-Rm(2);
            
            double Ra1 = sqrt ( Rm(0)*Rm(0) + Rm(1)*Rm(1) + l1*l1 + coreL*coreL );
            double Ra2 = sqrt ( Rm(0)*Rm(0) + Rm(1)*Rm(1) + l2*l2 + coreL*coreL );
            
            return (2.0e00*C4* Qt.transpose() * ( stress_straight_func (l2,Ra2,Rm,B) - stress_straight_func(l1,Ra1,Rm,B) ) * Qt );
        }
        
        //==================================================================================================
        // stress field at point R induced by straight segment x1 -> x2 using the coordinate dependent forms
        //==================================================================================================
        
        MatrixDim stress_CoorDpndt(const VectorDim & R, const VectorDim & x1 , const VectorDim & x2 , const MatrixDim & Qt, const VectorDim & segBurgers) const {
            
            VectorDim B = Qt * segBurgers;
            VectorDim Rm = Qt*(R-x1);
            
            double l1 = -Rm(2);
            double l2 = (x2-x1).norm()-Rm(2);
            
            double Ra1 = sqrt ( Rm(0)*Rm(0) + Rm(1)*Rm(1) + l1*l1 + coreL*coreL );
            double Ra2 = sqrt ( Rm(0)*Rm(0) + Rm(1)*Rm(1) + l2*l2 + coreL*coreL );
            
            return (2.0e00*C4* Qt.transpose() * ( stress_straight_func (l2,Ra2,Rm,B) - stress_straight_func(l1,Ra1,Rm,B) ) * Qt );
        }
        
        
        //===============================================================================================
        // function to calculate the stress field at one end of a straight dislocation segment using the coordinate-dependent forms
        //===============================================================================================
        
        MatrixDim stress_straight_func (const double& l1, const double& R1 , const VectorDim& Rm, const VectorDim& B) const {
            
            double p1 = Rm(0);    double p2 = Rm(1);          //double p3 = Rm(2);
            double b1 = B(0);     double b2 = B(1);           double b3 = B(2);
            
            double a = coreL;
            
            double RaxRaPlusa1 = R1*(R1+l1);
            double Ra2 = R1*R1;
            double Ra3 = R1*R1*R1;
            double a2 = a*a;
            double x2 = p1*p1;
            double y2 = p2*p2;
            
            MatrixDim temp;
            
            temp(0,0) = b1*p2/RaxRaPlusa1*(1.0e00 +(x2+a2)/Ra2+(x2+a2)/RaxRaPlusa1) +b2*p1/RaxRaPlusa1* (1.0 -(x2+a2)/Ra2-(x2+a2)/RaxRaPlusa1);
            
            temp(0,1) = -b1*p1/RaxRaPlusa1*(1.0 -y2/Ra2- y2/RaxRaPlusa1) + b2*p2/RaxRaPlusa1*  (1.0 -x2/Ra2-x2/RaxRaPlusa1);
            
            temp(0,2) = -b1*p1*p2/(Ra3) + b2*(x2/(Ra3)-nu/R1+(1.0 -nu)*0.5 *a2/(Ra3)) + b3*p2*(1.0 -nu)/RaxRaPlusa1*(1.0 +0.5 *a2/Ra2+0.5 * a2/RaxRaPlusa1);
            
            temp(1,1) = -b1*p2/RaxRaPlusa1*(1.0 -(y2+a2)/Ra2-(y2+a2)/ RaxRaPlusa1) -b2*p1/RaxRaPlusa1*  (1.0 +(y2+a2)/Ra2+(y2+a2)/RaxRaPlusa1) ;
            
            temp(1,2) = b1*(nu/R1-y2/(Ra3)-(1.0 -nu)*0.5 *a2/(Ra3))  + b2*p1*p2/(Ra3) - b3*p1*(1.0 -nu)/RaxRaPlusa1* (1.0 +0.5 *a2/Ra2+0.5 *a2/RaxRaPlusa1);
            
            temp(2,2) = b1*(2.0 *nu*p2/RaxRaPlusa1*(1.0 +0.5 *a2/Ra2+0.5 *a2/RaxRaPlusa1)+p2*l1/(Ra3))-b2*(2.0 *nu*p1/RaxRaPlusa1*(1.0 +0.5 *a2/Ra2+0.5 *a2 /RaxRaPlusa1)+p1*l1/(Ra3));
            
            temp(1,0)=temp(0,1);
            temp(2,0)=temp(0,2);
            temp(2,1)=temp(1,2);
            
            return temp;
        }
        
        //====================================================================================================
        // function to calculate the infinite medium displacement field of the 3 straight segments
        //====================================================================================================
        VectorDim displacement (const VectorDim& Rfield, const VectorDim& S) const {
            
            //return ( displacementStraight(Rfield,S,sourceP,sourcePatL) + displacementStraight(Rfield,S,sourcePatL,sinkPatL) + displacementStraight(Rfield,S,sinkPatL,sinkP) /*+
            //         (-piercingCompensation(Rfield,S)/(4.0*M_PI)*Burgers ) */);
            //return ( displacementStraight(Rfield,S,sourcePatL,sinkPatL) + (-piercingCompensation(Rfield,S)/(4.0*M_PI)*Burgers ) );
            
            return ( displacementStraight(Rfield,S,sourceP,sinkP) + displacement_Triangular(Rfield,sinkP,sourceP,sourcePatL) + displacement_Triangular(Rfield,sinkP,sourcePatL,sinkPatL) );
        }
        
        //====================================================================================================
        // function to calculate the infinite medium displacement field induced by a straight segment
        //====================================================================================================
        VectorDim displacementStraight (const VectorDim& Rfield, const VectorDim& S , const VectorDim& x1, const VectorDim& x2 ) const {
            
            VectorDim temp(VectorDim::Zero());
            
            VectorDim Ra = x1 - Rfield;    		  VectorDim Rb = x2 - Rfield;
            
            //VectorDim la = Ra.normalized();                VectorDim lb = Rb.normalized();
            
            VectorDim la  = VectorDim::Zero();             VectorDim lb  = VectorDim::Zero();
            if (Ra.norm()>1.0e-7)  la = Ra.normalized();
            if (Rb.norm()>1.0e-7)  lb = Rb.normalized();
            
            VectorDim tab = (x2-x1).normalized();
            
            //------------------ Solid angle ---------------
            double la_dot_lb = la.dot(lb);
            VectorDim la_cross_lb = la.cross(lb);
            
            
            double omega = 0.0e00;
            
            if (la_cross_lb.norm()>1.0e-7) omega = 2.0 * atan ( (S.dot(la_cross_lb)*sqrt((1-la_dot_lb)/(1+la_dot_lb)) / la_cross_lb.norm()) /
                                                               (1.0-(S.dot(la))- ( S.dot(la.cross(la_cross_lb))*sqrt((1-la_dot_lb)/(1+la_dot_lb)) / la_cross_lb.norm() ) ) );
            
            if (omega!= 0.0e00 && omega/omega != 1.0 ) {
                std::cout  << "Error in Solid angle calculation. OMEGA = " << omega << std::endl;
                assert(0&& "Error in Solid angle calculations ");
            }
            
            temp =  (-(omega/(4.0*M_PI))*Burgers) - (C3*C4*( f_function(tab,Ra,Rb))) +  (C4*(g_function(la,lb)  ));
            
            return temp;
        }
        
        //=========================================================
        // auxiliary function for displacement calculations
        //=========================================================
        VectorDim f_function (VectorDim t, VectorDim Ra , VectorDim Rb) const {
            
            double ra = sqrt(Ra.squaredNorm()+std::pow(coreL,2.0));
            double rb = sqrt(Rb.squaredNorm()+std::pow(coreL,2.0));
            
            VectorDim temp = log( (rb+Rb.dot(t)) / (ra+Ra.dot(t))  ) * Burgers.cross(t);
            
            double mag = temp.norm();
            
            if (mag!=0.0e00 && (mag/mag)!=1.0e00 ) {
                std::cout  << "Error in f_function calculation. f_function = "  << std::endl;
                std::cout  << temp.transpose() << std::endl;
                std::cout  << t.transpose() << std::endl;
                std::cout  << Ra.transpose() << std::endl;
                std::cout  << Rb.transpose() << std::endl;
                assert(0&&"Error in f_function calculation @ VirtualBoundarySlipSurface");
            }
            
            return temp;
        }
        
        //=========================================================
        // auxiliary function for displacement calculations
        //=========================================================
        VectorDim g_function ( VectorDim la , VectorDim lb) const {
            
            VectorDim temp = VectorDim::Zero();
            
            if ((la.cross(lb)).norm() > 1.0e-7) temp = (Burgers.dot(la.cross(lb))) /  (1.0 + la.dot(lb)) * (la+lb) ;
            
            double mag = temp.norm();
            
            if (mag!=0.0e00 && (mag/mag)!=1.0e00 ) {
                std::cout  << "Error in g_function calculation. g_function = "  << std::endl;
                std::cout  << temp.transpose() << std::endl;
                std::cout  << la.transpose() << std::endl;
                std::cout  << lb.transpose() << std::endl;
                assert(0&&"Error in g_function calculation @ VirtualBoundarySlipSurface");
            }
            
            return temp;
        }
        
        
        //==================================================================
        // function to calculate the displacement for a triangular loop
        //==================================================================
        
        VectorDim displacement_Triangular (const VectorDim& Rfield, const VectorDim& x0 , const VectorDim& x1, const VectorDim& x2) const {
            
            VectorDim temp(VectorDim::Zero());
            
            VectorDim Ra = x0 - Rfield;    		  VectorDim Rb = x1 - Rfield;    		  VectorDim Rc = x2 - Rfield;
            
            //VectorDim la = Ra.normalized();    		  VectorDim lb = Rb.normalized();                VectorDim lc = Rc.normalized();
            
            VectorDim tab = (x1-x0).normalized();         VectorDim tbc = (x2-x1).normalized();    VectorDim tca = (x0-x2).normalized();
            
            VectorDim N = (tab.cross(tbc)).normalized();
            
            VectorDim la  = VectorDim::Zero();             VectorDim lb  = VectorDim::Zero();        VectorDim lc  = VectorDim::Zero();
            
            if (Ra.norm()>1.0e-7)  la = Ra.normalized();
            if (Rb.norm()>1.0e-7)  lb = Rb.normalized();
            if (Rc.norm()>1.0e-7)  lc = Rc.normalized();
            
            double sign = 1.0e00;
            
            if      ( std::abs(la.dot(N)) > 1.0e-7 )  sign =  la.dot(N)/std::abs(la.dot(N));
            else if ( std::abs(lb.dot(N)) > 1.0e-7 )  sign =  lb.dot(N)/std::abs(lb.dot(N));
            else if ( std::abs(lc.dot(N)) > 1.0e-7 )  sign =  lc.dot(N)/std::abs(lc.dot(N));
            
            
            //----------- the solid angle calculations -----------
            //double a = acos(lb.dot(lc));            double b = acos(la.dot(lc));            double c = acos(la.dot(lb));
            double temp0;
            
            temp0 = lb.dot(lc);   if (temp0>1.0e00) temp0=1.0e00;       if (temp0<-1.0e00) temp0=-1.0e00;
            double a = acos(temp0);
            
            temp0 = la.dot(lc);   if (temp0>1.0e00) temp0=1.0e00;       if (temp0<-1.0e00) temp0=-1.0e00;
            double b = acos(temp0);
            
            temp0 = la.dot(lb);   if (temp0>1.0e00) temp0=1.0e00;       if (temp0<-1.0e00) temp0=-1.0e00;
            double c = acos(temp0);
            
            
            double s = (a+b+c)/2.0;
            
            double E4 = sqrt(std::abs (  tan(0.5*s) * tan(0.5*(s-a)) *  tan(0.5*(s-b)) * tan(0.5*(s-c)) ) );   // the absolute value is because some times one of the terms tan(0.5*(...-..)) give -1.0e-15 for example instead of zero
            
            double omega = -4.0* sign * atan(E4);
            
            if (omega!= 0.0e00 && omega/omega != 1.0 ){
                std::cout  << "Error in Solid angle calculation. OMEGA = " << omega << " " << E4 << " " << a << " " << lb.dot(lc) << " " <<
                b << " " << la.dot(lc) << " " <<
                c << " " << la.dot(lb) << " " << sign << std::endl;
                std::cout << la.transpose() << "  :  " << lb.transpose() << "  :  " << lc.transpose() <<std::endl;
                
                assert(0&& "Error in Solid angle calculations ");
            }
            
            //------------ the displacement expression ----------------------------
            
            temp =  (-(omega/(4.0*M_PI))*Burgers)
            -(C3*C4*( f_function(tab,Ra,Rb) +  f_function(tbc,Rb,Rc) + f_function(tca,Rc,Ra) ))
            +  (C4*( g_function(la,lb) + g_function(lb,lc) + g_function(lc,la) ));
            
            return temp;
        }
        
    };
    
    // init static data
    //    template <typename DislocationSegmentType, short unsigned int dim>
    template <typename DislocationSegmentType>
    const double VirtualBoundarySlipSurface<DislocationSegmentType>::L=10000000.0;
    
    template <typename DislocationSegmentType>
    double VirtualBoundarySlipSurface<DislocationSegmentType>::mu ;
    
    template <typename DislocationSegmentType>
    double VirtualBoundarySlipSurface<DislocationSegmentType>::nu ;
    
    template <typename DislocationSegmentType>
    double VirtualBoundarySlipSurface<DislocationSegmentType>::C1 ;
    
    template <typename DislocationSegmentType>
    double VirtualBoundarySlipSurface<DislocationSegmentType>::C2 ;
    
    template <typename DislocationSegmentType>
    double VirtualBoundarySlipSurface<DislocationSegmentType>::C3 ;
    
    template <typename DislocationSegmentType>
    double VirtualBoundarySlipSurface<DislocationSegmentType>::C4 ;
    
    /**************************************************************************/
    /**************************************************************************/
} // namespace model
#endif


