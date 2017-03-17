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
#include <stdlib.h>     /* atoi,atof */
#include <Eigen/Dense>
#include <model/FEM/FiniteElement.h>
#include <model/FEM/Elements/LagrangeElement.h>

using namespace model;

/**************************************************************************/
/**************************************************************************/
struct MyFunction : public EvalFunction<MyFunction>
{
    
    constexpr static int rows=1;
    constexpr static int cols=1;
    
    
    const double dt;
    
    MyFunction(const double& dt_in) : dt(dt_in){}
    
    /**********************************************************************/
    template<typename ElementType, typename BaryType>
    double operator() (const ElementType& ele, const BaryType& bary) const
    {/*!@param[in] elem the element
      * @param[in] bary the barycentric cooridinate
      *\returns the constant c.
      */
        
        const Eigen::Matrix<double,ElementType::dim,1> P=ele.position(bary);
        
        return dt*(P.norm()<0.1);
    }
    
};


/**************************************************************************/
template <int dim>
void generalTheta(const double& theta,
                  const double& D,
                  const int& nSteps)
{/*!Conside the PDE of the function \f$f=f(x,t)\f$:
  * \f[
  * \dot{f}=\nabla\left(D\nabla f\right)+p
  * \f]
  * Using the general-\f$\theta\f$ approximation the time discretization:
  * \f[
  * \frac{f_{n+1}-f_n}{dt}=\theta\left[\nabla(D\nabla f_{n+1})+p_{n+1}\right]+(1-\theta)
  * \left[\nabla(D\nabla f_{n})+p_{n}\right]
  * \f]
  * Let \f$\delta\f$ be the test operator, and \f$(X,Y)=\int_\Omega XY\, d\Omega\f$
  * indicate integration over the domain \f$\Omega\f$. The Galerkin
  * form of the equation above is obtained by multiplying it by \f$\delta f\f$,
  * and integration over the domain:
  * \f[
  * (\delta f,f_{n+1})=(\delta f,f_{n})
  * +(\delta f,dt\,\theta\,\nabla(D\nabla f_{n+1}))+(\delta f,dt\,\theta\,p_{n+1})
  * +(\delta f,dt\,(1-\theta)\,\nabla(D\nabla f_n))+(\delta f,dt\,(1-\theta)\,p_n)
  * \f]
  * We can reduce the order of differentiation in the second and ourth terms on the rhs using
  * integration by parts. This leads to:
  * \f[
  * (\delta f,f_{n+1})=(\delta f,f_{n})
  * +\langle\delta f,dt\,\theta\,D\nabla f_{n+1}\cdot N\rangle-(\nabla\delta f,dt\,\theta\,D\nabla f_{n+1})+(\delta f,dt\,\theta\,p_{n+1})
  * +\langle\delta f,dt\,(1-\theta)\,D\nabla f_n\cdot N\rangle-(\nabla\delta f,dt\,(1-\theta)\,D\nabla f_n)+(\delta f,dt\,(1-\theta)\,p_n)
  * \f]
  * Rearranging:
  * \f[
  * \underbrace{
  * (\delta f,f_{n+1})
  * -\langle\delta f,dt\,\theta\,D\nabla f_{n+1}\cdot N\rangle
  * +(\nabla\delta f,dt\,\theta\,D\nabla f_{n+1})
  * }_{\mathbf{A}\mathbf{f}_{n+1}}
  * =\underbrace{
  * (\delta f,f_{n})
  * +\langle\delta f,dt\,(1-\theta)\,D\nabla f_n\cdot N\rangle-(\nabla\delta f,dt\,(1-\theta)\,D\nabla f_n)+(\delta f,dt\,\theta\,p_{n+1}+dt\,(1-\theta)\,p_n)
  * }_{\mathbf{b}_n}
  * \f]
  * Here we used th divergence theorem and introduced the boundary integral
  * \f$<X,Y>=\int_{\partial\Omega} XY\, d\partial\Omega\f$.
  * Here the trial function is \f$f_{n+1}\f$, while \f$f_n\f$ is known. Therefore
  * - the lhs is a weak form bilinear in the test funtion and trial function.
  * - the rhs is a weak forms linear in the test function.
  *
  * Minimal code:
  * \code
  *
  * \endcode
  * Upon FEM discretization, the bilinear weak forms gives rise to the "stiffness" matrix \f$\mathbf{A}\f$,
  * while the linear weak forms contribute to the "nodal force vector" \f$\mathbf{b}\f$. The solution is found
  * solving the system of equations
  * \f[
  * \mathbf{A}\mathbf{f}_{n+1}=\mathbf{b}_n
  * \f]
  */
    
    const double dt=1.0/nSteps;
    const int SForder=1;
    const int meshID=dim;
    
    SimplicialMesh<dim> mesh;
    // Read the mesh files ./T/T_meshID.txt and ./N/N_meshID.txt
    mesh.readMesh(meshID);
    
    // Create a FiniteElement object on the mesh
    typedef LagrangeElement<dim,SForder> ElementType;
    typedef FiniteElement<ElementType> FiniteElementType;
    FiniteElementType fe(mesh);
    
    auto rho=fe.template trial<'r',1>();
    rho=0.0;
    SequentialOutputFile<'F',1> uFile;
    uFile<<rho.onDomain();
    
    // Define integration domain
    auto dV=fe.template domain<EntireDomain,3,GaussLegendre>();
    auto lhs1=(test(rho),rho)*dV;
    auto lhs2=(test(grad(rho)),theta*dt*D*grad(rho))*dV;

    auto rhs1=(test(rho),eval(rho))*dV;
    auto rhs2=(test(grad(rho)),eval((1.0-theta)*dt*D*grad(rho)))*dV;
    auto rhs3=(test(rho),MyFunction(dt))*dV;
    auto weakProblem(lhs1+lhs2==rhs1-rhs2+rhs3); //  weak problem
    
    const double solverTolerance=1.0e-6;
    for (int k=0;k<nSteps;++k)
    {
        std::cout<<"general-theta STEP "<<k<<std::endl;
        weakProblem.assemble();
        rho=weakProblem.solve(solverTolerance);
        
        SequentialOutputFile<'F',1> uFile;
        uFile<<rho.onDomain(); // Output function on domain
    }
    
}


int main(int argc, char** argv)
{
    
    int dim=1;          // dimensionality
    double theta=0.5;   // 0=forward Euler, 1=backward Euler, 0.5=Crank-Nicolson
    int nSteps=100;     // number of time steps
    double D=0.001;     // diffusion coefficient
    
    if (argc>1)
    {
        dim=atoi(argv[1]);
    }
    if (argc>2)
    {
        theta=atof(argv[2]);
    }
    if (argc>3)
    {
        nSteps=atoi(argv[3]);
    }
    if (argc>4)
    {
        D=atof(argv[4]);
    }
    
    assert(theta>=0.0 && theta<=1.0 && "time integration parameter theta must be 0<=theta<=1");
    
    if(dim==1)
    {
        generalTheta<1>(theta,D,nSteps);
    }
    else if(dim==2)
    {
        generalTheta<2>(theta,D,nSteps);
    }
    else
    {
        assert(0 && "ONLY CASES dim=1 and dim=2 are implemented.");
    }
    
    return 0;
}



///**************************************************************************/
//template <int dim>
//void forwardEuler(const double& D,
//                  const int& nSteps)
//{/*!Conside the PDE of the function \f$f=f(x,t)\f$:
//  * \f[
//  * \dot{f}=\nabla\left(D\nabla f\right)+p
//  * \f]
//  * Using the forward-Euler approximation of the time derivative, the PDE reads:
//  * \f[
//  * \frac{f_{n+1}-f_n}{dt}=\nabla(D\nabla f_n)+p_n
//  * \f]
//  * Let \f$\delta\f$ be the test operator, and \f$(X,Y)=\int_\Omega XY\, d\Omega\f$
//  * indicate integration over the domain \f$\Omega\f$. The Galerkin
//  * form of the equation above is obtained by multiplying it by \f$\delta f\f$,
//  * and integration over the domain:
//  * \f[
//  * (\delta f,f_{n+1})=(\delta f,f_{n})+(\delta f,dt\, \nabla(D\nabla f_n))+(\delta f,dt\,p_n)
//  * \f]
//  * We can reduce the order of differentiation in the second term on the rhs using
//  * integration by parts. This leads to:
//  * \f[
//  * \underbrace{(\delta f,f_{n+1})}_{\mathbf{A}\mathbf{f}_{n+1}}=\underbrace{(\delta f,f_{n})+<\delta f,dt\, D\nabla f_n\cdot N>-(\nabla\delta f,dt\, D\nabla f_n)+(\delta f,dt\,p_n)}_{\mathbf{b}_n}
//  * \f]
//  * Here we used th divergence theorem and introduced the boundary integral
//  * \f$<X,Y>=\int_{\partial\Omega} XY\, d\partial\Omega\f$.
//  * Here the trial function is \f$f_{n+1}\f$, while \f$f_n\f$ is known. Therefore
//  * - \f$(\delta f,f_{n+1})\f$ is a weak form bilinear in the test funtion and trial function.
//  * - \f$(\delta f,f_{n})\f$, \f$<\delta f,dt\, D\nabla f_n\cdot N>\f$,\f$(\nabla\delta f,dt\, D\nabla f_n)\f$, and \f$(\delta f,dt\,p_n)\f$
//  * are weak forms linear in the test function.
//  *
//  * Minimal code:
//  * \code
//  * const int dim=1;                                        // set problem dimensionality
//  * SimplicialMesh<dim> mesh;                               // create mesh
//  * mesh.readMesh(meshID);                                  // read mesh from file
//  *
//  * typedef LagrangeElement<dim,SForder> ElementType;       // define type of elements
//  * typedef FiniteElement<ElementType> FiniteElementType;   // defined FiniteElement
//  * FiniteElementType fe(mesh);                             // create FiniteElement from mesh
//  *
//  * auto f=fe.trial<'f',1>();                               // create trial function from FiniteElement
//  * auto gradf=grad(f);                                     // create gradient of f
//  * auto dV=fe.domain<EntireDomain,3,GaussLegendre>();      // create integration domain and integration method
//  * auto lhs=(test(f),f)*dV;                               // left-hand-side of the weak problem
//  * auto rhs=(test(f),eval(f))*dV+(test(gradf),eval(dt*D*gradf))*dV+(test(f),dp*p)    // right-hand-side of weakProblem
//  * auto weakProblem(lhs==rhs);                              // create weak problem from lhs and rhs
//  * f=weakProblem.solve();                                  // solve weak problem
//  * \endcode
//  * Upon FEM discretization, the bilinear weak forms gives rise to the "stiffness" matrix \f$\mathbf{A}\f$,
//  * while the linear weak forms contribute to the "nodal force vector" \f$\mathbf{b}\f$. The solution is found
//  * solving the system of equations
//  * \f[
//  * \mathbf{A}\mathbf{f}_{n+1}=\mathbf{b}_n
//  * \f]
//  */
//
//    const double dt=1.0/nSteps;
//    const int SForder=1;
//    const int meshID=dim;
//
//    SimplicialMesh<dim> mesh;
//    // Read the mesh files ./T/T_meshID.txt and ./N/N_meshID.txt
//    mesh.readMesh(meshID);
//
//    // Create a FiniteElement object on the mesh
//    typedef LagrangeElement<dim,SForder> ElementType;
//    typedef FiniteElement<ElementType> FiniteElementType;
//    FiniteElementType fe(mesh);
//
//    auto rho=fe.template trial<'r',1>();
//    rho=0.0;
//    SequentialOutputFile<'F',1> uFile;
//    uFile<<rho.onDomain();
//
//    // Define integration domain
//    auto dV=fe.template domain<EntireDomain,3,GaussLegendre>();
//    auto lhs=(test(rho),rho)*dV;
//    auto rhs1=(test(rho),eval(rho))*dV;
//    auto rhs2=(test(grad(rho)),eval(dt*D*grad(rho)))*dV;
//    auto rhs3=(test(rho),MyFunction(dt))*dV;
//    auto weakProblem(lhs==rhs1-rhs2+rhs3); //  weak problem
//
//    const double solverTolerance=1.0e-6;
//    for (int k=0;k<nSteps;++k)
//    {
//        std::cout<<"forwardEuler STEP "<<k<<std::endl;
//        weakProblem.assemble();
//        rho=weakProblem.solve(solverTolerance);
//
//        SequentialOutputFile<'F',1> uFile;
//        uFile<<rho.onDomain(); // Output function on domain
//    }
//
//}

///**************************************************************************/
//template <int dim>
//void backwardEuler(const double& D,
//                   const int& nSteps)
//{/*!Conside the PDE of the function \f$f=f(x,t)\f$:
//  * \f[
//  * \dot{f}=D\nabla^2f+p
//  * \f]
//  * Using the forward-Euler approximation of the time derivative, the PDE reads:
//  * \f[
//  * \frac{f_{n+1}-f_n}{dt}=D\nabla^2f_{n+1}+p_{n+1}
//  * \f]
//  * Let \f$\delta\f$ be the test operator, and \f$(X,Y)=\int_\Omega XY\, d\Omega\f$
//  * indicate integration over the domain \f$\Omega\f$. The Galerkin
//  * form of the equation above is obtained by multiplying it by \f$\delta f\f$,
//  * and integration over the domain:
//  * \f[
//  * (\delta f,f_{n+1})-(\delta f,dt\, D\nabla^2f_{n+1})=(\delta f,f_{n})+(\delta f,dt\,p_{n+1})
//  * \f]
//  * We can reduce the order of differentiation in the second term on the rhs using
//  * integration by parts. This leads to:
//  * \f[
//  * \underbrace{(\delta f,f_{n+1})-<\delta f,dt\, D\nabla f_{n+1}\cdot N>+(\nabla\delta f,dt\, D\nabla f_{n+1})}_{\mathbf{A}\mathbf{f}_{n+1}}=\underbrace{(\delta f,f_{n})+(\delta f,dt\,p_n)+(\delta f,dt^2\,\left.\frac{\partial p}{\partial t}\right|_n)}_{\mathbf{b}_n}
//  * \f]
//  * Here we used th divergence theorem and introduced the boundary integral
//  * \f$<X,Y>=\int_{\partial\Omega} XY\, d\partial\Omega\f$.
//  * Here the trial function is \f$f_{n+1}\f$, while \f$f_n\f$ is known. Therefore
//  * - \f$(\delta f,f_{n+1})\f$, \f$<\delta f,dt\, D\nabla f_{n+1}\cdot N>\f$, and \f$(\nabla\delta f,dt\, D\nabla f_{n+1})\f$ are weak form bilinear in the test funtion and trial function.
//  * - \f$(\delta f,f_{n})\f$ and \f$(\delta f,dt\,p_n)\f$
//  * are weak forms linear in the test function.
//  *
//  * Upon FEM discretization, the bilinear weak forms gives rise to the "stiffness" matrix \f$\mathbf{A}\f$,
//  * while the linear weak forms contribute to the "nodal force vector" \f$\mathbf{b}\f$. The solution is found
//  * solving the system of equations
//  * \f[
//  * \mathbf{A}\mathbf{f}_{n+1}=\mathbf{b}
//  * \f]
//  */
//
//    const double dt=1.0/nSteps;
//    const int SForder=1;
//    const int meshID=dim;
//
//    SimplicialMesh<dim> mesh;
//    // Read the mesh files ./T/T_meshID.txt and ./N/N_meshID.txt
//    mesh.readMesh(meshID);
//
//    // Create a FiniteElement object on the mesh
//    typedef LagrangeElement<dim,SForder> ElementType;
//    typedef FiniteElement<ElementType> FiniteElementType;
//    FiniteElementType fe(mesh);
//
//    auto rho=fe.template trial<'r',1>();
//    rho=0.0;
//    SequentialOutputFile<'B',1> uFile;
//    uFile<<rho.onDomain();
//
//    // Define integration domain
//    auto dV=fe.template domain<EntireDomain,3,GaussLegendre>();
//    auto lhs1=(test(rho),rho)*dV;
//    auto lhs2=(test(grad(rho)),dt*D*grad(rho))*dV;
//
//    auto rhs1=(test(rho),eval(rho))*dV;
//    auto rhs2=(test(rho),MyFunction(dt))*dV;
//    auto weakProblem(lhs1+lhs2==rhs1+rhs2); //  weak problem
//
//    const double solverTolerance=1.0e-6;
//    for (int k=0;k<nSteps;++k)
//    {
//        std::cout<<"backwardEuler STEP "<<k<<std::endl;
//        weakProblem.assemble();
//        rho=weakProblem.solve(solverTolerance);
//
//
//        SequentialOutputFile<'B',1> uFile;
//        uFile<<rho.onDomain(); // Output function on domain
//    }
//
//}

