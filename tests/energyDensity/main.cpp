#include <iostream>
#include <DislocatedMaterial.h>
#include <StraightDislocationSegment.h>
#include <DislocationStress.h>
#include <TextFileParser.h>
#include <vector>


using namespace model;
const double Wcai(const Eigen::Vector3d R, const Eigen::Vector3d b1, const Eigen::Vector3d b2, const Eigen::Vector3d t1, const Eigen::Vector3d t2);

int main(int argc, char** argv)
{
    std::cout<<"a2 is "<<DislocationStress<3>::a2<<std::endl;
	std::cout<<"nu is "<<DislocatedMaterial<3,Isotropic>::Nu<<std::endl;
			
    const bool computeSelfEnergiesOnly= argc > 1 ? atoi(argv[1]) : false;
    
    if(computeSelfEnergiesOnly)
    {
        const int nGP(TextFileParser("input.txt").readScalar<int>("nGP",false));
        
        const Eigen::MatrixXd segments(TextFileParser("input.txt").readMatrixCols<double>("S",9,false));
        std::ofstream pointsFile("points.txt");
        for(int r=0;r<segments.rows();++r)
        {
            const Eigen::Vector3d P0(segments.template block<1,3>(r,0));
            const Eigen::Vector3d P1(segments.template block<1,3>(r,3));
            const Eigen::Vector3d b(segments.template block<1,3>(r,6));
            const double length((P1-P0).norm());
            const Eigen::Vector3d t((P1-P0)/length);
            
            StraightDislocationSegment sds(P0,P1,b,length,t);
            
            const double Wns=Wcai(P0-P0,b,b,t,t)-Wcai(P1-P0,b,b,t,t);
            
            const double dL(length/nGP);
            for(int k=0;k<nGP;++k)
            {
                const Eigen::Vector3d x(P0+(0.5+k)*dL*t);
                const double e(sds.elasticInteractionEnergy(x,t,b));
                pointsFile<<r<<" "<<k<<" "<<x.transpose()<<" "<<e<<" "<<dL<<" "<<Wns<<"\n";
            }
            
        }
    }
    else
    {
        const int nGP(TextFileParser("input.txt").readScalar<int>("nGP",false));
        
        const Eigen::MatrixXd segments(TextFileParser("input.txt").readMatrixCols<double>("S",9,false));
		const int paires=segments.rows()/2;
        std::ofstream pointsFile("points.txt");
        for(int r=0;r<paires;++r)
        {
            const Eigen::Vector3d P0(segments.template block<1,3>(2*r,0));
            const Eigen::Vector3d P1(segments.template block<1,3>(2*r,3));
            const Eigen::Vector3d b1(segments.template block<1,3>(2*r,6));
            const double length1((P1-P0).norm());
            const Eigen::Vector3d t1((P1-P0)/length1);
            
            StraightDislocationSegment sds1(P0,P1,b1,length1,t1);
			
			const Eigen::Vector3d P2(segments.template block<1,3>(2*r+1,0));
            const Eigen::Vector3d P3(segments.template block<1,3>(2*r+1,3));
            const Eigen::Vector3d b2(segments.template block<1,3>(2*r+1,6));
            const double length2((P3-P2).norm());
            const Eigen::Vector3d t2((P3-P2)/length2);
            
            StraightDislocationSegment sds2(P2,P3,b2,length2,t2);
            
            const double Wself1=Wcai(P0-P0,b1,b1,t1,t1)-Wcai(P1-P0,b1,b1,t1,t1);
			const double Wself2=Wcai(P2-P2,b2,b2,t2,t2)-Wcai(P3-P2,b2,b2,t2,t2);
			const double W12=Wcai(P3-P1,b1,b2,t1,t2)+Wcai(P2-P0,b1,b2,t1,t2)-Wcai(P3-P0,b1,b2,t1,t2)-Wcai(P2-P1,b1,b2,t1,t2);
			const double Wns=Wself1+Wself2+W12;
            
            const double dL1(length1/nGP);
			const double dL2(length2/nGP);
            for(int k=0;k<nGP;++k)
            {
                const Eigen::Vector3d x1(P0+(0.5+k)*dL1*t1);
				const Eigen::Vector3d x2(P2+(0.5+k)*dL2*t2);
                const double e11(sds1.elasticInteractionEnergy(x1,t1,b1));
				const double e12(sds2.elasticInteractionEnergy(x1,t1,b1));
				
				const double e21(sds1.elasticInteractionEnergy(x2,t2,b2));
				const double e22(sds2.elasticInteractionEnergy(x2,t2,b2));
                pointsFile<<2*r<<" "<<k<<" "<<x1.transpose()<<" "<<e11+e12<<" "<<dL1<<" "<<Wself1+0.5*W12<<"\n";
				pointsFile<<2*r+1<<" "<<k<<" "<<x2.transpose()<<" "<<e21+e22<<" "<<dL2<<" "<<Wself2+0.5*W12<<"\n";
            }
		}
    }
    
    
    return 0;
}

const double Wcai(const Eigen::Vector3d R, const Eigen::Vector3d b1, const Eigen::Vector3d b2, const Eigen::Vector3d t1, const Eigen::Vector3d t2)
{
	if( (t1.cross(t2)).norm() < 1.e-3 )
	{
		// Parallel case
		const double Ra=sqrt(R.squaredNorm()+DislocationStress<3>::a2);
		const double W0=DislocatedMaterial<3,Isotropic>::C2;
		const double b1t1=b1.dot(t1);
		const double b1R=b1.dot(R);
		const double b1b2=b1.dot(b2);
		const double b2t1=b2.dot(t1);
		const double b2R=b2.dot(R);
		const double Rt1=R.dot(t1);
		const double logRaRt1=log(Ra+Rt1);
		const double RaRt=pow(Ra,2.0)-pow(Rt1,2.0);
		
		return W0*((b1t1*b2R+b1R*b2t1-((2.0-DislocatedMaterial<3,Isotropic>::Nu)*b1t1*b2t1+b1b2)*Rt1)*logRaRt1
		/*  */ +((1.0-DislocatedMaterial<3,Isotropic>::Nu)*b1t1*b2t1+b1b2)*Ra
		/*  */ -(b1R-Rt1*b1t1)*(b2R-Rt1*b2t1)/RaRt*Ra
		/*  */ +DislocationStress<3>::a2*( (1.0+DislocatedMaterial<3,Isotropic>::Nu)*b1t1*b2t1-2.0*b1b2 )/2.0/RaRt*Ra);
		
	}
	else{
		// Non-parallel case
		const Eigen::Vector3d u(t1.cross(t2));
		const Eigen::Vector3d v(u.cross(t1)); 
		const Eigen::Vector3d vp(t2.cross(u));
		const double Ra=sqrt(R.squaredNorm()+DislocationStress<3>::a2);
		const double W0=DislocatedMaterial<3,Isotropic>::C2/u.squaredNorm();
		const double b1t1=b1.dot(t1);
		const double b1t2=b1.dot(t2);
		const double b1b2=b1.dot(b2);
		const double b2t1=b2.dot(t1);
		const double b2t2=b2.dot(t2);
		const double t1t2=t1.dot(t2);
		const double b1u=b1.dot(u);
		const double b2u=b2.dot(u);
		const double b1v=b1.dot(v);
		const double b2v=b2.dot(v);
		const double b1vp=b1.dot(vp);
		const double b2vp=b2.dot(vp);
		const double A1=DislocatedMaterial<3,Isotropic>::C1*b1t1*b2t2+2.0*DislocatedMaterial<3,Isotropic>::Nu*b2t1*b1t2;
		const double A2=(b1b2+b1t1*b2t1)*t1t2;
		const double A2p=(b1b2+b1t2*b2t2)*t1t2;
		const double A3p=(b1u*b2vp+b2u*b1vp)*t1t2/u.squaredNorm();
		const double A3=(b1u*b2v+b2u*b1v)*t1t2/u.squaredNorm();
		const double A4=(b1t1*b2v+b1t2*b2vp)*t1t2;
		const double A5=2.0* (b1.cross(u)).dot( (b2.cross(u)) )*t1t2/u.squaredNorm();
		const double logRaRt2=log(Ra+R.dot(t2));
		const double logRaRt1=log(Ra+R.dot(t1));
		const double uua2Ru2=u.squaredNorm()*DislocationStress<3>::a2 + pow(R.dot(u),2.0);
		const double act=atan( ((1.0+t1t2)*Ra+R.dot(t1+t2))/sqrt(uua2Ru2) );
		
		return W0*((A1-A2p)*R.dot(vp)*logRaRt2+A3p*R.dot(u)*logRaRt2
		/*  */  +(A1-A2)*R.dot(v)*logRaRt1+A3*R.dot(u)*logRaRt1 + A4*Ra
		/*  */	+(A1-A5)*(u.squaredNorm()*DislocationStress<3>::a2 + 2.0*pow(R.dot(u),2.0))/sqrt(uua2Ru2)*act);	
	}
	
	
}


