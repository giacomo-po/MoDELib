/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// The following line is needed for shape functions or order 5 or more
//#define EIGEN_STACK_ALLOCATION_LIMIT 1000000


#include <iostream>
#include <model/IO/SequentialOutputFile.h>
#include <model/FEM/FiniteElement.h>
#include <model/FEM/Boundaries/AtXmin.h>
#include <model/FEM/Boundaries/AtXmax.h>
#include <model/FEM/BoundaryConditions/Fix.h>
#include <model/FEM/Domains/SubDomain.h>
#include <model/FEM/SymTensorToVoigt.h>
#include <model/FEM/TensorToVoigt.h>
#include <HookeanMaterial.h>


using namespace model;


template<int dim>
struct VoigtIdentity
{
    
    //    static Eigen::Map<Eigen::Matrix<double,dim*dim,1> getI()
    //    {
    //
    //        Eigen::Matrix<double,dim,dim>::Identity()
    //
    //    }
    
    
    static const Eigen::Matrix<double,dim,dim> I;
    static const Eigen::Map<const Eigen::Matrix<double,dim*dim,1>> voigtI;
    
};

template<int dim>
const Eigen::Matrix<double,dim,dim> VoigtIdentity<dim>::I=Eigen::Matrix<double,dim,dim>::Identity();
template<int dim>
const Eigen::Map<const Eigen::Matrix<double,dim*dim,1>> VoigtIdentity<dim>::voigtI=Eigen::Map<const Eigen::Matrix<double,dim*dim,1>>(VoigtIdentity<dim>::I.data());


template<typename FiniteElementType>
struct DeformationGradient : public EvalFunction<DeformationGradient<FiniteElementType> >
{
    constexpr static int dim=FiniteElementType::dim;
    constexpr static int rows=dim*dim;
    constexpr static int cols=1;
    typedef TrialFunction<'u',dim,FiniteElementType> TrialDisplacementType;
    const TrialDisplacementType& U;                    // displacement
    const TrialGrad<TrialDisplacementType> gradU;       // distortion
    
    /**********************************************************************/
    DeformationGradient(const TrialDisplacementType& u_in) :
    /* init */ U(u_in),
    /* init */ gradU(grad(U))
    {
        
    }
    
    /**********************************************************************/
    template<typename ElementType, typename BaryType>
    const Eigen::Matrix<double,rows,cols> operator() (const ElementType& ele, const BaryType& bary) const
    {
        return eval(gradU)(ele,bary)+VoigtIdentity<dim>::voigtI;
        //        return VoigtIdentity<dim>::voigtI;
    }

    /**********************************************************************/
    template<typename ElementType, typename BaryType>
    const Eigen::Matrix<double,dim,dim> asTensor(const ElementType& ele, const BaryType& bary) const
    {
        return Eigen::Map<const Eigen::Matrix<double,dim,dim,Eigen::RowMajor>>(this->operator()(ele,bary).data());
    }

    
};


template<typename FiniteElementType>
struct InversePlasticDeformation : public EvalFunction<InversePlasticDeformation<FiniteElementType> >
{
    constexpr static int dim=FiniteElementType::dim;
    constexpr static int rows=dim*dim;
    constexpr static int cols=1;
    const TrialFunction<'A',dim*dim,FiniteElementType>& Fp;       // Plastic distortion
    
    /**********************************************************************/
    InversePlasticDeformation(const TrialFunction<'A',dim*dim,FiniteElementType>& Fp_in) :
    /* init */ Fp(Fp_in)
    {
        
    }
    
    /**********************************************************************/
    template<typename ElementType, typename BaryType>
    const Eigen::Matrix<double,dim,dim> asTensor(const ElementType& ele, const BaryType& bary) const
    {// inv(Fp)
        const Eigen::Matrix<double,rows,cols> fp(eval(Fp)(ele,bary));
        const Eigen::Map<const Eigen::Matrix<double,dim,dim,Eigen::RowMajor>> Mfp(fp.data());
        return Mfp.inverse();
        
    }
    
    /**********************************************************************/
    template<typename ElementType, typename BaryType>
    const Eigen::Matrix<double,rows,cols> operator() (const ElementType& ele, const BaryType& bary) const
    {// inv(Fp)
        return Eigen::Map<const Eigen::Matrix<double,rows,cols>>(asTensor(ele,bary).data());
    }
    
    
};

template<typename FiniteElementType>
struct ElasticDeformation : public EvalFunction<ElasticDeformation<FiniteElementType> >
{
    constexpr static int dim=FiniteElementType::dim;
    constexpr static int rows=dim*dim;
    constexpr static int cols=1;
    const DeformationGradient<FiniteElementType>& F;
    const InversePlasticDeformation<FiniteElementType>& Gp;       // Plastic distortion
    
    /**********************************************************************/
    ElasticDeformation(const DeformationGradient<FiniteElementType>& F_in,
                               const InversePlasticDeformation<FiniteElementType>& Gp_in) :
    /* init */ F(F_in),
    /* init */ Gp(Gp_in)
    {
        
    }
    
    /**********************************************************************/
    template<typename ElementType, typename BaryType>
    const Eigen::Matrix<double,dim,dim> asTensor(const ElementType& ele, const BaryType& bary) const
    {// F=Fe*Fp, so Fe=F*inv(Fp)
        return F.asTensor(ele,bary)*Gp.asTensor(ele,bary);
    }
    
    /**********************************************************************/
    template<typename ElementType, typename BaryType>
    const Eigen::Matrix<double,rows,cols> operator() (const ElementType& ele, const BaryType& bary) const
    {// F=Fe*Fp, so Fe=F*inv(Fp)
        return Eigen::Map<const Eigen::Matrix<double,rows,cols>>(asTensor(ele,bary).data());
    }
};


/**************************************************************************/
/**************************************************************************/
template <typename FiniteElementType, template<int> class MaterialType>
struct JacobianMatrix : public EvalFunction<JacobianMatrix<FiniteElementType,HookeanMaterial> >
{
    constexpr static int dim=FiniteElementType::dim;
    constexpr static int rows=dim*dim;
    constexpr static int cols=dim*dim;
    static constexpr int StrainVoigtSize=dim*(dim+1)/2;
    typedef Eigen::Matrix<double,StrainVoigtSize,1> StrainVoigtVector;
    typedef Eigen::Matrix<double,StrainVoigtSize,StrainVoigtSize> StrainVoigtMatrix;
    const DeformationGradient<FiniteElementType>& F;
    const InversePlasticDeformation<FiniteElementType>& Gp;
    const MaterialType<dim> material;
    
    /**********************************************************************/
    template<typename...ArgTypes>
    JacobianMatrix(const DeformationGradient<FiniteElementType>& F_in,
                    const InversePlasticDeformation<FiniteElementType>& Gp_in,
                    const ArgTypes&...materialArgs) :
    /* init */ F(F_in),
    /* init */ Gp(Gp_in),
    /* init */ material(materialArgs...)
    {
        //            std::cout<<"Constant Constructor 1"<<std::endl;
        //            std::cout<<"c="<<c<<std::endl;
    }
    
    
    /**********************************************************************/
    template<typename ElementType, typename BaryType>
    const Eigen::Matrix<double,rows,cols> operator() (const ElementType& ele, const BaryType& bary) const
    {/*!@param[in] elem the element
      * @param[in] bary the barycentric cooridinate
      *\returns the current stiffness C, which in general is a funciton of C0 and grad(u)
      */
        
        // see this for reshaping (https://eigen.tuxfamily.org/dox/group__TutorialReshapeSlicing.html)
        //        Eigen::Matrix<double,dim,dim> Fp(fp(ele,bary));
        
        //        StrainVoigtVector Ee; //FINISH HERE
        //        const StrainVoigtMatrix C(material.C(Ee));
        //
        
        const Eigen::Matrix<double,dim,dim> f(F.asTensor(ele,bary));
        const Eigen::Matrix<double,dim,dim> gp(Gp.asTensor(ele,bary));
        const Eigen::Matrix<double,dim,dim> fe(f*gp);
        
        const Eigen::Matrix<double,dim,dim> ee=0.5*(fe.transpose()*fe-VoigtIdentity<dim>::I);
        
        Eigen::Matrix<double,rows,cols> temp(Eigen::Matrix<double,rows,cols>::Zero());


        for(int i=0;i<dim;++i)
        {
            for(int J=0;J<dim;++J)
            {
                const int voigtRow=TensorToVoigt<dim,dim>::t2v(i,J);
                for(int k=0;k<dim;++k)
                {
                    for(int L=0;L<dim;++L)
                    {
                        const int voigtCol=TensorToVoigt<dim,dim>::t2v(k,L);

                        for(int a=0;a<dim;++a)
                        {
                            for(int b=0;b<dim;++b)
                            {
                                const int ab=SymTensorToVoigt<dim>::t2v(a,b);

                                for(int g=0;g<dim;++g)
                                {
                                    for(int d=0;d<dim;++d)
                                    {
                                        const int gd=symTensorToVoigt<dim>::t2v(g,d);

                                        temp(voigtRow,voigtCol)+= C(ab,gd)*Fe(i,a)*Gp(J,b)*Fe(k,g)*Gp(L,d);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        //
        //
        //        // second term S_ab*\delta_ik*Gp_Lb*Gp_Ja
        //        for(int i=0;i<dim;++i)
        //        {
        //            for(int J=0;J<dim;++J)
        //            {
        //                const voigtRow=tensorToVoigt<dim>(i,J);
        //                const int& k=i; // avoid computing delta_ik
        //
        //                for(int L=0;L<dim;++L)
        //                {
        //                    const voigtCol=tensorToVoigt<dim>(k,L);
        //
        //                    for(int a=0;a<dim;++a)
        //                    {
        //                        for(int b=0;b<dim;++b)
        //                        {
        //                            temp(voigtRow,voigtCol)+= S(a,b)*Gp(L,b)*Gp(J,a);
        //                        }
        //                    }
        //                }
        //            }
        //        }
        
        return temp;
    }
    
};




int main(int argc, char** argv)
{
    const int dim=2;
    
    // Take meshID as a user input
    int meshID(0);
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }
    
    // Create a 3d-SimplicialMesh object
    SimplicialMesh<dim> mesh;
    // Read the mesh files ./T/T_meshID.txt and ./N/N_meshID.txt
    mesh.readMesh(meshID);
    
    
    // Create a FiniteElement object on the mesh
    typedef LagrangeElement<dim,2> ElementType;
    typedef FiniteElement<ElementType> FiniteElementType;
    FiniteElementType fe(mesh);
    
    /**************************************************************************/
    // Define some constants for the 3d elastic problem
    const double mu =75.6;  // GPa (for Cu)
    const double lam=119.9; // GPa (for Cu)
    const double C11(lam+2.0*mu);
    const double C12(lam);
    const double C44(mu);
    
    Eigen::Matrix<double,dim*(dim+1)/2,dim*(dim+1)/2> C;
    C  << C11, C12, 0.0,
    /***/ C12, C11, 0.0,
    /***/ 0.0, 0.0, C44;
    C/=mu; // make dimensionless
    
    //    C  << C11, C12, C12, 0.0, 0.0, 0.0,
    //    /***/ C12, C11, C12, 0.0, 0.0, 0.0,
    //    /***/ C12, C12, C11, 0.0, 0.0, 0.0,
    //    /***/ 0.0, 0.0, 0.0, C44, 0.0, 0.0,
    //    /***/ 0.0, 0.0, 0.0, 0.0, C44, 0.0,
    //    /***/ 0.0, 0.0, 0.0, 0.0, 0.0, C44;
    //    C/=mu; // make dimensionless
    
    /**************************************************************************/
    // Define trial function (displacement field) u and related expressions
    auto u=fe.trial<'u',dim>();   // displacement field u=[u1; u2; u3]
    auto Fp=fe.trial<'A',dim*dim>();   // displacement field u=[u1; u2; u3]
    DeformationGradient<FiniteElementType> F(u);
    InversePlasticDeformation<FiniteElementType> Gp(Fp);
    ElasticDeformation<FiniteElementType> Fe(F,Gp);
    JacobianMatrix<FiniteElementType,HookeanMaterial> J(Fe,Gp,C);

    
    u=Eigen::VectorXd::Random(20402);
    Fp.setNodeDof(VoigtIdentity<dim>::voigtI);

    
    //    auto b=grad(u);             // displacement gradient b=[u1,1; u1,2; u2,1; u2,2]
    //auto e=def(u);              // engineering strain e=[u1,1; u2,2; u3,3; u1,2+u2,1; u2,3+u3,2; u1,3+u3,1]
    //auto s=C*e;                 // stress field s=[s11; s22; s33; s12; s23; s13]
    
    /**************************************************************************/
    // Create the BilinearWeakForm bWF_u=int(test(e)^T*s)dV
    //    auto dV=fe.domain<EntireDomain,4,GaussLegendre>();
    //    auto bWF_u=(test(b),B*b)*dV;
    
    //    for(int i=0;i<dim;++i)
    //    {
    //        for(int j=0;j<dim;++j)
    //        {
    //            std::cout<<SymTensorToVoigt<dim>::t2v(i,j)<<" ";
    //        }
    //        std::cout<<std::endl;
    //    }
    //
    //    std::cout<<std::endl;
    //    for(int i=0;i<dim;++i)
    //    {
    //        for(int j=0;j<dim;++j)
    //        {
    //            std::cout<<TensorToVoigt<dim,dim>::t2v(i,j)<<" ";
    //        }
    //        std::cout<<std::endl;
    //    }
    
    //    /**************************************************************************/
    //    // Create the LinearWeakForm lWF_1=int(test(u)^T*f)ndA
    //    // f is a constant traction vector.
    //    //    Eigen::Matrix<double,3,1> f((Eigen::Matrix<double,3,1>()<<0.0,0.000,0.00).finished());
    //    auto f=make_constant((Eigen::Matrix<double,3,1>()<<0.0,0.000,0.00).finished());
    //    auto dA_1=fe.boundary<AtXmax<2>,3,GaussLegendre>();
    //    auto lWF_1=(test(u),f)*dA_1;
    //
    //    /**************************************************************************/
    //    // Create the LinearWeakForm lWF_2=int(test(u)^T*p)ndA
    //    // p is a constant hydrostatic tensor. The traction vector will be t=p*n;
    //    Eigen::Matrix<double,3,3> p(-0.01*Eigen::Matrix<double,3,3>::Identity());
    //    //    auto p=make_constant(-0.01*Eigen::Matrix<double,3,3>::Identity());
    //    auto ndA_2=fe.boundary<ExternalBoundary,3,GaussLegendre>();
    //    auto lWF_2=(test(u),make_constant(p))*ndA_2;
    //
    //    /**************************************************************************/
    //    // Create the LinearWeakForm lWF_1=int(test(u)^T*f)ndA
    //    // a is a constant boby force vector.
    //    //    Eigen::Matrix<double,3,1> a((Eigen::Matrix<double,3,1>()<<0.0,0.000,-0.00001).finished());
    //    auto a=make_constant((Eigen::Matrix<double,3,1>()<<0.0,0.000,-0.00001).finished());
    //    auto lWF_3=(test(u),a)*dV;
    //
    //    /**************************************************************************/
    //    // Create the WeakProblem
    //    auto weakProblem(bWF_u==lWF_1+lWF_2+lWF_3); //  weak problem
    //
    //    /**************************************************************************/
    //    // Set up Dirichlet boundary conditions
    //    // Create a list of nodes having x(0)=x0_min, where x0_min is the minimum value among the fe nodes
    //    const size_t nodeList_0=fe.createNodeList<AtXmin<2>>();
    //    // Fix those those nodes
    //    Fix fix;
    //    u.addDirichletCondition(nodeList_0,fix,{1, 1, 1}); // fix u1, u2, u3
    //
    //    // Create a list of nodes having x(2)=x2_max, where x2_max is the maximum value among the fe nodes
    //    //const size_t nodeList_1=fe.createNodeList<AtXmax<2>>();
    //    //u.addDirichletCondition(nodeList_1,fix,{0, 0, 1}); // fix only u3
    //
    //
    //    /**************************************************************************/
    //    // Assemble
    //    //    weakProblem.assembleWithLagrangeConstraints();
    //    //    weakProblem.assembleWithPenaltyConstraints(1000000.0);
    //    weakProblem.assembleWithMasterSlaveConstraints();
    //
    //    /**************************************************************************/
    //    // Solve for u using current value as guess solution
    //    const double solverTolerance=1.0e-6;
    //    u=weakProblem.solveWithGuess(solverTolerance,u);
    //    
    /**************************************************************************/
    // Output displacement and stress on external mesh faces
    SequentialOutputFile<'U',1> uFile;
    uFile<<u.onDomain();
    
    SequentialOutputFile<'F',1> fFile;
    
    const Eigen::Matrix<double,dim+1,dim+1> vertexBary(Eigen::Matrix<double,dim+1,dim+1>::Identity());
    for(const auto& ele : fe.elements())
    {
        for (int v=0;v<dim+1;++v)
        {
            fFile<<ele.second.position(vertexBary.col(v)).transpose()<<" "
            <<Fe(ele.second,vertexBary.col(v)).transpose()<<"\n";
        }
    }
    
    //
    //    SequentialOutputFile<'F',1> sFile;
    //    sFile<<s.onBoundary();
    
    return 0;
}




