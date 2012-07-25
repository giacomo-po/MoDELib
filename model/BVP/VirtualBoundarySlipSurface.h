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
#include <model/Dislocations/DislocationNetworkTraits.h>
#include <model/Dislocations/GlidePlanes/GlidePlaneObserver.h>

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
    	const GlidePlaneSharedPtrType pGlidePlane;
        
    public:
        
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
        /* init list                                     */ coreL(ds.coreL),
        /* init list                                     */ gp_H(ds.pGlidePlane->height),
        /* init list                                     */ gp_N(ds.pGlidePlane->planeNormal){
            
            
            setElasticConstants(ds);
            
            //----------calculate the displacement field induced by those segments on the boundary nodes ----
            for (unsigned int dN=0; dN<ds.shared.domain.nodeContainer.size();++dN){
                if(ds.shared.domain.nodeContainer[dN].isBoundaryNode) {
                    if(ds.shared.domain.nodeContainer[dN].triIDs.size() == 0) assert(0&&"Error: Boundary node without triangle element index array");                    
                    ds.shared.domain.nodeContainer[dN].uVir+= displacement(ds.shared.domain.nodeContainer[dN].P , ds.shared.domain.triContainer[ds.shared.domain.nodeContainer[dN].triIDs[0]]->outNormal);  
                }
            }
        }
        
        
        //=========================================================================
        // function to set the elastic constants
        //========================================================================
        // template <typename DislocationSegmentType>
        void setElasticConstants (const DislocationSegmentType& ds) {
            mu=ds.shared.material.mu;
            nu=ds.shared.material.nu;
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
            
            VectorDim temp =  (P_Prj - dn.get_P()).normalized() ;
            
            temp = dn.boundaryNormal.dot(temp) / std::abs(dn.boundaryNormal.dot(temp)) * temp;
            
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
            
            
            temp =  (-(omega/(4.0*M_PI))*Burgers) - (C3*C4*( f_function(tab,Ra,Rb))) +  (C4*(g_function(la,lb)  ));
            
            return temp;
        }
        
        //=========================================================
        // auxiliary function for displacement calculations 
        //=========================================================
        VectorDim f_function (VectorDim t, VectorDim Ra , VectorDim Rb) const {
            
            double ra = sqrt(Ra.squaredNorm()+std::pow(coreL,2.0));
            double rb = sqrt(Rb.squaredNorm()+std::pow(coreL,2.0));
            
            return ( log( (rb+Rb.dot(t)) / (ra+Ra.dot(t))  ) * Burgers.cross(t) );
        }
        
        //=========================================================
        // auxiliary function for displacement calculations 
        //=========================================================
        VectorDim g_function ( VectorDim la , VectorDim lb) const {
            
            VectorDim temp = VectorDim::Zero();
            
            if ((la.cross(lb)).norm() > 1.0e-7) temp = (Burgers.dot(la.cross(lb))) /  (1.0 + la.dot(lb)) * (la+lb) ;
            
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
            double a = acos(lb.dot(lc));            double b = acos(la.dot(lc));            double c = acos(la.dot(lb));
            
            double s = (a+b+c)/2.0;
            
            double E4 = sqrt(std::abs (  tan(0.5*s) * tan(0.5*(s-a)) *  tan(0.5*(s-b)) * tan(0.5*(s-c)) ) );   // the absolute value is because some times one of the terms tan(0.5*(...-..)) give -1.0e-15 for example instead of zero
            
            double omega = -4.0* sign * atan(E4);
            
            if (omega!= 0.0e00 && omega/omega != 1.0 ) std::cout  << "Error in Solid angle calculation. OMEGA = " << omega << " " << E4 << " " << a << " " << b << " " << c << " " << sign << std::endl;
            
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






/*
 //==================================================================================
 // add a compensation for the solid angle term if a vector from the field point in the direction of S is piercing the slipped area
 //==================================================================================
 double piercingCompensation (const VectorDim& Rfield, const VectorDim& S) const {
 
 double dOmega = 0.0;
 
 if (isPiercing(Rfield,S)) {std::cout<<Rfield.transpose()<< std::endl;  dOmega = 4.0*M_PI;}
 
 return dOmega; 
 }
 
 //==================================================================================
 // function to check if a vector from the field point in the direction of S is piercing the slipped area
 //==================================================================================
 bool isPiercing(const VectorDim& Rfield, const VectorDim& S) const {
 
 bool itPierces = false;
 
 //-------- the point at which the vector from Rfield in S direction intersects the slipped area --
 double l=0.0;
 double temp = S.dot(gp_N);
 if(std::abs(temp) > 1.0e-8) l = (gp_P-Rfield).dot(gp_N)/temp;
 
 if (l>0.0e00) {      //  if -ve this means that no piercing
 VectorDim P = Rfield + l*S;
 itPierces = isInsdieTriangle(P,sourceP,sourcePatL,sinkP) || isInsdieTriangle(P,sinkP,sourcePatL,sinkPatL) ;
 }
 
 return itPierces;
 }
 
 //==================================================================================
 // function to check if a point is included inside a triangle, given its 3 vertexes positions 
 //==================================================================================
 
 bool isInsdieTriangle (const VectorDim& P, const VectorDim& x0, const VectorDim& x1, const VectorDim& x2) const {
 
 Eigen::Matrix<unsigned int,2,1> pP_index = find2DProjectionPlane(x0,x1,x2);
 
 Eigen::Matrix<double,3,1>::Index ii;
 Eigen::Matrix<double,3,1> Bary = getBarycentric(P , pP_index(0) , pP_index(1) , x0,x1,x2);
 double baryMin = Bary.minCoeff(&ii);
 
 return (baryMin>=0.0e00);
 }
 
 //======================================================================================================
 // function to find the 2D projection plane for the triangle. On this plane, the Barycentric coordinates
 // in 2D will be calculated. If the returned values ix = 2 , iy = 3 so the projection plane is y-z plane
 //=======================================================================================================
 Eigen::Matrix<unsigned int,2,1> find2DProjectionPlane(const VectorDim& x0, const VectorDim& x1, const VectorDim& x2) const{
 
 VectorDim uV,outNormal;
 std::vector<int> normal , not_normal;
 Eigen::Matrix<unsigned int,2,1> ppindx;
 
 outNormal = ((x1-x0).cross(x2-x0)).normalized();
 
 for(int i=0; i<3; i++){
 uV = VectorDim::Zero();     uV(i) = 1.0;
 if(std::abs(uV.dot(outNormal)) < 0.15e00){     // the plane normal makes angle of > 80 degrees with direction x_i
 normal.push_back(i);
 }
 else{not_normal.push_back(i);}
 }
 
 switch (normal.size()){
 case 2:
 ppindx(0) = normal[0];
 ppindx(1) = normal[1];
 break;
 case 1:
 ppindx(0) = normal[0];
 ppindx(1) = not_normal[0];
 break;
 case 0:
 ppindx(0) = not_normal[0];
 ppindx(1) = not_normal[1];
 break;
 }
 return ppindx;
 }
 
 //=====================================================================================
 // function to return the barycentric coordinates of point P w.r.t. triangle (x0,x1,x2)
 //=====================================================================================
 Eigen::Matrix<double,3,1>  getBarycentric (VectorDim P , unsigned int ix  , unsigned int iy , VectorDim x0 , VectorDim x1, VectorDim x2  ) const
 {
 Eigen::Matrix<double,3,1> bary;
 double V0;
 
 V0 = getVol(x0,x1,x2 ,ix, iy);
 
 bary(0) = getVol(P , x1, x2, ix, iy)/V0;
 bary(1) = getVol(x0, P , x2, ix, iy)/V0;
 bary(2) = getVol(x0, x1, P , ix, iy)/V0;
 
 return bary;
 }
 
 //======================================================================
 // return the volume of the tetradedron, given its points coordinates in order
 //======================================================================
 double getVol(VectorDim a, VectorDim b, VectorDim c , int ix, int iy) const {
 
 MatrixDim temp;
 temp << a(ix), b(ix), c(ix),
 a(iy), b(iy), c(iy),
 1.0   , 1.0   , 1.0;  
 return temp.determinant()/2.0;
 }
 */


//        //=========================================================================================
//        // Constructor to be used in the restart case
//        //========================================================================================
//        template <typename SharedType>
//        VirtualBoundarySlipSurface(const double gp_H_in, const VectorDim gp_N_in, const VectorDim sourceP_in, const VectorDim sourceN_in, const VectorDim sinkP_in,
//                                   const VectorDim sinkN_in, const VectorDim Burgers_in, const double coreL_in ,const SharedType* sharedPtr) : 
//        /*                                               */  pGlidePlane(GlidePlaneObserver<dim,DislocationSegmentType>::findExistingGlidePlane(gp_N_in,gp_H_in)),
//        /*                                               */ sourceP(sourceP_in),
//        /*                                               */   sinkP(sinkP_in),
//        /*                                               */ sourceN(sourceN_in),
//        /*                                               */   sinkN(sinkN_in),
//        /*                                               */   sourcePatL(sourceP+L*sourceN),
//        /*                                               */     sinkPatL(  sinkP+L*sinkN),
//        /*                                               */   Burgers(Burgers_in),
//        /*                                               */   coreL(coreL_in),
//        /*                                               */  // mu(sharedPtr->material.mu),
//        /*                                               */ //  nu(sharedPtr->material.nu),
//        /*                                               */   gp_H(gp_H_in),
//        /*                                               */   gp_N(gp_N_in){
//            
//            
//            
//            //std::cout<<"Creating VirtualBoundarySlipSurface "<<std::endl;
//            
//            setElasticConstants(sharedPtr);
//            
//            /*
//             mu=sharedPtr->material.mu;
//             nu=sharedPtr->material.nu;
//             
//             C1=1.0-nu;
//             C2=mu/(4.0*M_PI*C1);
//             C3=1.0-2.0*nu;
//             C4=1.0/(8.0*M_PI*C1);
//             */
//            //----------calculate the displacement field induced by those segments on the boundary nodes ----
//            for (unsigned int dN=0; dN<sharedPtr->domain.nodeContainer.size();++dN){
//                if(sharedPtr->domain.nodeContainer[dN].isBoundaryNode) {
//                    if(sharedPtr->domain.nodeContainer[dN].triIDs.size() == 0) assert(0&&"Error: Boundary node without triangle element index array");
//                    
//                    sharedPtr->domain.nodeContainer[dN].uVir+= displacement(sharedPtr->domain.nodeContainer[dN].P , 
//                                                                            sharedPtr->domain.triContainer[sharedPtr->domain.nodeContainer[dN].triIDs[0]]->outNormal);  
//                }
//            }
//            
//        }
