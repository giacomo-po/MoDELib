/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2014 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DiscontinuousLagrangeElement_H_
#define model_DiscontinuousLagrangeElement_H_

#include <vector>
#include <map>
#include <deque>

#include <Eigen/Dense>

#include <TypeTraits.h>
#include <CompareVectorsByComponent.h>
#include <CombinationWithRepetition.h>
#include <StarsAndBars.h>
#include <BarycentricTraits.h>
#include <IsoparametricMapping.h>
#include <LagrangeNode.h>
#include <MulticomponentExpander.h>
#include <Simplex.h>


namespace model
{
    
    /**************************************************************************/
	/**************************************************************************/
    /*!\brief Class template for a Lagrange finite element of degree p on a
     * d-dimensional Simplex.
     *
     * The degrees of freedom for a LagrangeElement are the functinal values 
     * computed at the nodes (LagrangeNode). In a d-dimensional Simplex,
     * each node is identified by a set of (d+1) positive integers \f$I_k\f$ 
     * which satisfy the Diophantine equation:
     *	\f[
     *		\sum_{k=0}^{d+1}I_k=p
     *	\f]
     * These integeres are found using the StarsAndBars algorithm.
     * The baricentric (area) coordinate of node \f$\{I_k\}\f$ are obtained as
     * \f$A_k=I_k/p\f$.
     *
     * The shape function of node \f$\{I_k\}\f$ are:
     *	\f[
     *		N_{\{I_k\}}(A_0,\dots A_{d})=\prod_{k=0}^{d+1}\frac{1}{I_k!}\prod_{i=0}^{I_k-1}\left(pA_k-i\right)
     *	\f]
     * This can also be written as:
     *	\f[
     *		N_{\{I_k\}}(A_0,\dots A_{d})=\prod_{i=0}^{d+1}\left([c_{i0},\dots,c_{ip}]\cdot[A_i^0,\dots,A_i^p]\right)
     *	\f]
     * The matrix \f$c_{ij}\f$ is stored in the field sfCoeffs;
     */
    template<int _dim,int degree, template<typename T> class MappingType=IsoparametricMapping>
	class DiscontinuousLagrangeElement : private std::vector<LagrangeNode<DiscontinuousLagrangeElement<_dim,degree,MappingType>>*>
    {
        
        typedef DiscontinuousLagrangeElement<_dim,degree,MappingType> ElementType;
        
    public:
        
        constexpr static int dim=_dim;
        constexpr static int order=degree;
        constexpr static int nodesPerElement=CombinationWithRepetition<_dim+1,degree>::value;
        constexpr static int dofPerNode(int c) { return c; }
        
        typedef LagrangeNode<ElementType> NodeType;
        
    private:
        
        /**********************************************************************/
        static Eigen::Matrix<double,1,Eigen::Dynamic> shapeFunctionCoeff(const int& N)
        {/*!@param[in] N the degree of the polynomial
          *\returns the coefficients of the polynomial
          * \f$\frac{1}{N!}\prod_{n=0}^{N-1}(px-n)\f$, where p is the degree of
          * the LagrangeElement.
          */
            Eigen::Matrix<int,1,Eigen::Dynamic> prod(1,1);
            prod<<1; // initialize the polynomial
            int f(1); // initialize the factorial
            
            for (int n=0;n<N;++n)
            {
                Eigen::Matrix<int,1,2> temp2(-n,degree);
                Eigen::Matrix<int,1,Eigen::Dynamic> temp1(1,prod.cols()+1);
                temp1.setZero();
                for (int i = 0;i < prod.cols();++i)
                {
                    for (int j = 0;j < temp2.cols();++j)
                    {
                        temp1(i + j) += prod(i) * temp2(j); // convolution
                    }
                }
                prod=temp1; // overwrite product
                if (n>0)
                {
                    f*=n; // compute factorial
                }
            }
            if (N>0)
            {
                f*=N; // add last term to factorial
            }
            return prod.cast<double>()/f;
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim+1,degree+1>  getSFcoeffs(const Eigen::Matrix<int,1,dim+1> bary)
        {/*!@param[in] bary baricenric coordinates of a LagrangeElement
          *\returns A matrix whose k-th row is the coefficients of polynomial
          * of the the k-th baricentric coordinate.
          */
            Eigen::Matrix<double,dim+1,degree+1>  temp(Eigen::Matrix<double,dim+1,degree+1>::Zero());
            for (int k=0;k<dim+1;++k)
            {
                temp.row(k).segment(0,bary(k)+1)=shapeFunctionCoeff(bary(k));
            }
            return temp;
        }
        
        /**********************************************************************/
        static std::vector<Eigen::Matrix<double,dim+1,degree+1>>  get_sfCoeffsVector()
        {/*!\returns a vector of Matrices of shape function coefficients. Each
          * entry in the vector corresponds to a node.
          */
            std::vector<Eigen::Matrix<double,dim+1,degree+1>> temp;
            for (int n=0;n<nodesPerElement;++n)
            {
                temp.push_back(getSFcoeffs(baryStarsAndBars.row(n)));
            }
            return temp;
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,1,dim+1> sfcTimesBary(const int& n, const Eigen::Matrix<double,1,dim+1>& bary)
        {
            Eigen::Matrix<double,1,dim+1> temp(Eigen::Matrix<double,1,dim+1>::Zero());
            for(int d=0;d<dim+1;++d) // this loop can be parallelized
            {
                Eigen::Matrix<double,degree+1,1> bPow(Eigen::Matrix<double,degree+1,1>::Ones());
                for (int p=1;p<degree+1;++p) // start at 1 since bPow(0)=1 from initialization
                {
                    //                    bPow(p)=std::pow(bary(d),p);
                    bPow(p)=bPow(p-1)*bary(d);
                }
                temp(d)=sfCoeffsVector[n].row(d).dot(bPow);
            }
            return temp;
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,1,dim+1> sfcTimesBaryDiff(const int& n, const Eigen::Matrix<double,1,dim+1>& bary)
        {
            Eigen::Matrix<double,1,dim+1> temp(Eigen::Matrix<double,1,dim+1>::Zero());
            for(int d=0;d<dim+1;++d) // this loop can be parallelized
            {
                
                Eigen::Matrix<double,degree+1,1> bPow(Eigen::Matrix<double,degree+1,1>::Ones());
                for (int p=2;p<degree+1;++p) // start at 1 since bPow(0)=0 (unused) and bPow(1)=1 from initialization
                {
                    bPow(p)=p*std::pow(bary(d),p-1);
                }
                //                temp(d)=sfCoeffs.row(d).dot(bPow);
                temp(d)=sfCoeffsVector[n].template block<1,degree>(d,1).dot(bPow.template segment<degree>(1));
            }
            return temp;
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim,dim> get_FL(const Simplex<dim,dim>& s)
        {
            Eigen::Matrix<double,dim,dim> temp;
            
            typename SimplexTraits<dim,0>::SimplexIDType last;
            last[0]=s.xID(dim); // a 1x1 matrix
            
            //! Compute linear part of the deformation gradient
            for (int k=0;k<dim;++k)
            {
                typename SimplexTraits<dim,0>::SimplexIDType current;
                current[0]=s.xID(k); // a 1x1 matrix
                temp.col(k) << s.mesh->template observer<0>().simplex(current).P0 - s.mesh->template observer<0>().simplex(last).P0;
            }
            return temp;
        }
        
        
    public:
        
        
        static const Eigen::Matrix<int   ,CombinationWithRepetition<_dim+1,degree>::value,_dim+1> baryStarsAndBars;
        static const Eigen::Matrix<double,CombinationWithRepetition<_dim+1,degree>::value,_dim+1> baryNodalCoordinates;
        static const std::vector<Eigen::Matrix<double,_dim+1,degree+1>>  sfCoeffsVector;
        
        //! A const reference to the Simplex that this element refers to
        const Simplex<dim,dim>& simplex;
        
        //! The matrix of nodal coordinates
        Eigen::Matrix<double,_dim,nodesPerElement> Xe; // this should be private
        
        
        /**********************************************************************/
        DiscontinuousLagrangeElement(const Simplex<dim,dim>& s,
                        std::deque<NodeType>& nodeContainer,
                        std::map<Eigen::Matrix<double,dim,1>, NodeType* const,
                        CompareVectorsByComponent<double,dim,float> >& nodeFinder,
                        std::map<size_t,const NodeType* const>& mesh2femIDmap) :
        /* init list */ simplex(s)
        {/*!@param[in] s A const reference to a Simplex<dim,dim>
          */
            
            typename SimplexTraits<dim,0>::SimplexIDType last;
            last[0]=s.xID(dim); // a 1x1 matrix

            const Eigen::Matrix<double,dim,dim> FL(get_FL(simplex));
            
            
            for (int n=0;n<nodesPerElement;++n)
            {
                // Place nodes linearly
                const Eigen::Matrix<double,dim,1> P(FL*BarycentricTraits<dim>::l2x(baryNodalCoordinates.row(n))+s.mesh->template observer<0>().simplex(last).P0);
                
//                typename std::map<Eigen::Matrix<double,dim,1>, NodeType* const, CompareVectorsByComponent<double,dim,float> >::const_iterator nIter(nodeFinder.find(P));
                
//                NodeType* pN;
                
//                if (nIter!=nodeFinder.end())
//                {// A node already exists at P.
//                    pN=nIter->second; // grab its pointer.
//                }
//                else
//                {// No nodes exist at P.
                    nodeContainer.emplace_back(P,nodeContainer.size()); // create a new one in nodeContainer
                    NodeType* pN=&*nodeContainer.rbegin(); // grab node pointer

//                    const bool success(nodeFinder.insert(std::make_pair(P,pN)).second); // insert pointer in nodeFinder
//                    assert(success && "NODE NOT INSERTED");
                
//                    size_t maxID=0;
//                    const double maxVal(baryNodalCoordinates.row(n).maxCoeff(&maxID));
//                    if(maxVal==1.0)
//                    {
//                        mesh2femIDmap.emplace(simplex.vertices()[maxID]->xID(0),pN);
//                    }
                
//                }
                
                this->emplace_back(pN); // add node to this
                
                Xe.col(n)=P; // insert node poision in matrix Xe
                
                

                
            }
            
        }
        
//        /**********************************************************************/
//        ~LagrangeElement()
//        {
//            for(int n=0;n<this->size();++n)
//            {
//                int nErased=this->operator[](n)->erase(this);
//                assert(nErased==1);
//            }
//        }
        
        /**********************************************************************/
        //template<typename TrialFunctionType>
        static Eigen::Matrix<double,1,nodesPerElement> sf(const Eigen::Matrix<double,1,dim+1>& bary)
        {/*!\param[in] bary the vector of barycentric coordinates
          *\returns a row vector of the the shape-functions, evaluated at bary.
          * This is:
          *	\f[
          *		N_{\{I_k\}}(A_0,\dots A_{d})=\prod_{i=0}^{d+1}\left([c_{i0},\dots,c_{ip}]\cdot[A_i^0,\dots,A_i^p]\right)
          *	\f]
          */
            Eigen::Matrix<double,1,nodesPerElement> temp;
            for(int n=0;n<nodesPerElement;++n)
            {
                temp(n)=sfcTimesBary(n,bary).prod();
            }
            return temp;
        }
        
        /**********************************************************************/
        template<typename TrialFunctionType>
        static typename TypeTraits<TrialFunctionType>::ShapeFunctionMatrixType sfm(const Eigen::Matrix<double,1,dim+1>& bary)
        {/*!\param[in] bary the vector of barycentric coordinates
          *\returns a row vector of the the shape-functions, evaluated at bary.
          * This is:
          *	\f[
          *		N_{\{I_k\}}(A_0,\dots A_{d})=\prod_{i=0}^{d+1}\left([c_{i0},\dots,c_{ip}]\cdot[A_i^0,\dots,A_i^p]\right)
          *	\f]
          */
            return MulticomponentExpander<TypeTraits<TrialFunctionType>::nComponents>::expandSF(sf(bary));
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim+1,nodesPerElement> diff(const Eigen::Matrix<double,1,dim+1>& bary)
        {/*!\param[in] bary the barycentric coordinate
          *\returns the value of the derivatives of this shape-function
          * with respect to the baricentric coordinates, computed at bary. This is:
          *	\f[
          *		dN_{\{I_k\}}/dA_j(A_0,\dots A_{d})
          *	\f]
          */
            Eigen::Matrix<double,dim+1,nodesPerElement> temp(Eigen::Matrix<double,dim+1,nodesPerElement>::Zero());
            for (int n=0;n<nodesPerElement;++n)
            {
                const Eigen::Matrix<double,1,dim+1> temp1(sfcTimesBary(n,bary));
                const Eigen::Matrix<double,1,dim+1> temp2(sfcTimesBaryDiff(n,bary));
                
                for(int d=0;d<dim+1;++d)
                {
                    Eigen::Matrix<double,1,dim+1> temp3(temp1);
                    temp3(d)=temp2(d);
                    temp(d,n)=temp3.prod();
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        static Eigen::Matrix<double,dim,nodesPerElement> gradS(const Eigen::Matrix<double,1,dim+1>& bary)
        {/*!\param[in] bary the barycentric coordinate
          *\returns the gradient of this shape-function
          * with respect to the standard coordinates, computed at bary. This is:
          *	\f[
          *		dN_{\{I_k\}}/dS_j(A_0,\dots A_{d})
          *	\f]
          */
            return BarycentricTraits<dim>::dLdX.transpose()*diff(bary);
        }
        
        /**********************************************************************/
        template<typename TrialFunctionType>
        typename TypeTraits<TrialFunctionType>::ShapeFunctionGradMatrixType sfmGrad(const Eigen::Matrix<double,1,dim+1>& bary) const
        {/*!\param[in] bary the vector of barycentric coordinates
          *\returns the gradient of the shape-functions, evaluated at bary.
          * This is the matrix :
          *	\f[
          *		M_{ij}=\frac{\partial N_j}{\partial x_i}
          *	\f]
          */
            Eigen::Matrix<double,dim,nodesPerElement> temp(Gs(bary).transpose()*gradS(bary));
            return MulticomponentExpander<TypeTraits<TrialFunctionType>::nComponents>::expandSFgrad(temp);
        }
        
        
        /**********************************************************************/
        template<typename TrialFunctionType>
        typename TypeTraits<TrialFunctionType>::ShapeFunctionDefMatrixType sfmDef(const Eigen::Matrix<double,1,dim+1>& bary) const
        {/*!\param[in] bary the vector of barycentric coordinates
          * \returns A matrix of the symmetric gradient of the shape-functions,
          * evaluated at bary. Symmetric gradients are interpreted in egineering sense.
          */
            static_assert(TypeTraits<TrialFunctionType>::nComponents==dim,"SYMMETRIC GRADIENT (DEF) CAN ONLY BE COMPUTED IF nComponents==dim.");
            Eigen::Matrix<double,dim,nodesPerElement> temp(Gs(bary).transpose()*gradS(bary));
            return MulticomponentExpander<TypeTraits<TrialFunctionType>::nComponents>::expandSFdef(temp);
        }
        
        /**********************************************************************/
        const NodeType& node(const size_t& n) const
        {/*!@param[in] the local node ID in the element
          *\returns a reference to the n-th LagrangeNode
          */
            assert(n<nodesPerElement && "NODE ID EXCEEDS nodesPerElement");
            return *this->operator[](n);
        }
        
        /**********************************************************************/
        NodeType& node(const size_t& n)
        {/*!@param[in] the local node ID in the element
          *\returns a reference to the n-th LagrangeNode
          */
            assert(n<nodesPerElement && "NODE ID EXCEEDS nodesPerElement");
            return *this->operator[](n);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,dim,1> position(const Eigen::Matrix<double,1,dim+1>& bary) const
        {/*!@param[in] bary the vector of barycentric coordinate
          * \returns the position vector corresponding to bary
          */
            return MappingType<ElementType>(*this).position(bary);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,_dim,_dim> Fs(const Eigen::Matrix<double,1,dim+1>& bary) const
        {/*!@param bary the vector of baricentric coordinates
          * \returns the Jacobian matrix dx_i / ds_j, evaulated at bary.
          */
            return MappingType<ElementType>(*this).Fs(bary);
        }
        
        /**********************************************************************/
        double absJ(const Eigen::Matrix<double,1,dim+1>& bary) const
        {/*!@param bary the vector of baricentric coordinates
          * \returns the absolute value of the determinant of the Jacobian
          * matrix dx_i / ds_j, evaulated at bary.
          */
            return MappingType<ElementType>(*this).absJ(bary);
        }
        
        /**********************************************************************/
        Eigen::Matrix<double,_dim,_dim> Gs(const Eigen::Matrix<double,1,dim+1>& bary) const
        {/*!@param bary the vector of baricentric coordinates
          * \returns the inverse Jacobian matrix ds_i / dx_j, evaulated at bary.
          */
            return MappingType<ElementType>(*this).Gs(bary);
        }

        /**********************************************************************/
        Eigen::Matrix<double,dim,1> jGN(const Eigen::Matrix<double,1,dim+1>& bary, const int& n) const
        {
            return MappingType<ElementType>(*this).jGN(bary,n);
        }

        /**********************************************************************/
        bool isBoundaryElement() const
        {/*!\returns true if *this is a boundary element
          */
            return simplex.isBoundarySimplex();
        }
        
        /**********************************************************************/
        std::vector<int> boundaryFaces() const
        {/*!\returns a vector of IDs of the faces of *this that are on the boundary
          */
            return simplex.boundaryFaces();
        }
        
    };
    
    
    // Declare static data members
    template<int dim,int degree, template<typename T> class MappingType>
    const Eigen::Matrix<int,CombinationWithRepetition<dim+1,degree>::value,dim+1> DiscontinuousLagrangeElement<dim,degree,MappingType>::baryStarsAndBars=StarsAndBars<dim+1,degree>::sAb();
    
    template<int dim,int degree, template<typename T> class MappingType>
    const Eigen::Matrix<double,CombinationWithRepetition<dim+1,degree>::value,dim+1> DiscontinuousLagrangeElement<dim,degree,MappingType>::baryNodalCoordinates=DiscontinuousLagrangeElement<dim,degree>::baryStarsAndBars.template cast<double>()/degree;
    
    template<int dim,int degree, template<typename T> class MappingType>
    const std::vector<Eigen::Matrix<double,dim+1,degree+1>>  DiscontinuousLagrangeElement<dim,degree,MappingType>::sfCoeffsVector=DiscontinuousLagrangeElement<dim,degree,MappingType>::get_sfCoeffsVector();
    
}	// close namespace
#endif



