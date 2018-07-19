

//double axialStrain0;
//double axialStrainRate;
//
//void initUserLoad()
//{
//
//}

MatrixDim userStress(const VectorDim& x) const
{

    // We want to add the torsion stress in a clinder along z
    
    // Compute radius of cylinder
    const double R(0.5*(this->network().mesh.xMax()(0)-this->network().mesh.xMin()(0)));
    //const double J(0.5*M_PI*pow(R,4)); // polar moment of inertia of cross section
    
    MatrixDim temp(-this->network().extStressController.externalStress()); // remove applied stress
    
    
    const double s_33=this->network().extStressController.externalStress()(2,2); // re-apply only zz component
    
    
    const double theta=atan2 (x(1),x(0));
    const double r=sqrt(x(1)*x(1)+x(0)*x(0));
    
    const double a=1.0;
    
    const double tau_31=-a*s_33*sin(theta);
    const double tau_32=+a*s_33*cos(theta);
    
    
    temp(2,2)+=s_33;
    temp(2,0)+=tau_31*r/R;
    temp(2,1)+=tau_32*r/R;
    
    temp(0,2)=temp(2,0);
    temp(1,2)=temp(2,1);
    

    
    return temp;
    
}
