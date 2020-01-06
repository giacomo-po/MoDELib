/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <FiniteElement.h>
#include <AtXmin.h>
#include <AtXmax.h>
#include <IntegrationDomain.h>

using namespace model;


struct TestIntegrator
{
    static constexpr int dim=3;
    static constexpr int qOrder=37;
    typedef Eigen::Matrix<double,dim,1> VectorDim;
    typedef LagrangeElement<dim,2> ElementType;
    typedef FiniteElement<ElementType> FiniteElementType;
    
    
    /**************************************/
    TestIntegrator(const int& meshID)
    {
        
        SimplicialMesh<dim> mesh;
        mesh.readMesh(meshID);

        FiniteElementType fe(mesh);
        const IntegrationDomain<FiniteElementType,1,qOrder,GaussLegendre> top(topBoundary(fe));
        const IntegrationDomain<FiniteElementType,1,qOrder,GaussLegendre> bottom(bottomBoundary(fe));

        std::cout<<"top.size()="<<top.size()<<std::endl;
        VectorDim forceTop(Eigen::Matrix<double,3,1>::Zero());
        top.integrate(this,forceTop,&TestIntegrator::testFunction);    // integrate the bvp correction
        std::cout<<forceTop<<std::endl;
        
        std::cout<<"bottom.size()="<<bottom.size()<<std::endl;
        VectorDim forceBottom(Eigen::Matrix<double,3,1>::Zero());
        bottom.integrate(this,forceBottom,&TestIntegrator::testFunction);    // integrate the bvp correction
        std::cout<<forceBottom<<std::endl;
    }
    
    /**************************************/
    template <typename FiniteElementType>
    IntegrationDomain<FiniteElementType,1,qOrder,GaussLegendre> topBoundary(const FiniteElementType& fe)
    {
        IntegrationDomain<FiniteElementType,1,qOrder,GaussLegendre> temp;
        
        for(const auto& eIter : fe.elements())
        {
            if(eIter.second.isBoundaryElement())
            {
                const std::vector<int> boundaryFaces=eIter.second.boundaryFaces();
                for (unsigned int f=0;f<boundaryFaces.size();++f) // loop ever bonudary faces of the current Elements
                {
                    bool isExternalBoundaryFace(true);
                    std::array<const Simplex<FiniteElementType::dim,0>*, SimplexTraits<FiniteElementType::dim,FiniteElementType::dim-1>::nVertices> vertices=eIter.second.simplex.child(boundaryFaces[f]).vertices();
                    for(unsigned int v=0;v<vertices.size();++v) // loop over vertices of the current face
                    {
                        isExternalBoundaryFace *= (std::fabs(vertices[v]->P0(2)-fe.xMax()(2))<FLT_EPSILON); // check if the current vertices satisfies operator()
                    }
                    if(isExternalBoundaryFace)
                    {
                        temp.emplace_back(&eIter.second,boundaryFaces[f]);
                    }
                }
            }
        }
        
        return temp;
    }
    
    /**************************************/
    template <typename FiniteElementType>
    IntegrationDomain<FiniteElementType,1,qOrder,GaussLegendre> bottomBoundary(const FiniteElementType& fe)
    {
        IntegrationDomain<FiniteElementType,1,qOrder,GaussLegendre> temp;
        
        for(const auto& eIter : fe.elements())
        {
            if(eIter.second.isBoundaryElement())
            {
                const std::vector<int> boundaryFaces=eIter.second.boundaryFaces();
                for (unsigned int f=0;f<boundaryFaces.size();++f) // loop ever bonudary faces of the current Elements
                {
                    bool isExternalBoundaryFace(true);
                    std::array<const Simplex<FiniteElementType::dim,0>*, SimplexTraits<FiniteElementType::dim,FiniteElementType::dim-1>::nVertices> vertices=eIter.second.simplex.child(boundaryFaces[f]).vertices();
                    for(unsigned int v=0;v<vertices.size();++v) // loop over vertices of the current face
                    {
                        isExternalBoundaryFace *= (std::fabs(vertices[v]->P0(2)-fe.xMin()(2))<FLT_EPSILON); // check if the current vertices satisfies operator()
                    }
                    if(isExternalBoundaryFace)
                    {
                        temp.emplace_back(&eIter.second,boundaryFaces[f]);
                    }
                }
            }
        }
        
        return temp;
    }
    
    /**********************************************************************/
    Eigen::Matrix<double,dim+1,1> face2domainBary(const Eigen::Matrix<double,dim,1>& b1,
                                                  /*                                         */ const int& boundaryFace) const
    {
        // Transform to barycentric coordinate on the volume, adding a zero on the boundaryFace-face
        Eigen::Matrix<double,dim+1,1> bary;
        for (int k=0;k<dim;++k)
        {
            bary((k<boundaryFace)? k : k+1)=b1(k);
        }
        bary(boundaryFace)=0.0;
        return bary;
    }
    
    /**********************************************************************/
    VectorDim testFunction(const Eigen::Matrix<double,dim-1,1>& a1, const ElementType& ele, const int& boundaryFace) const
    {
        const Eigen::Matrix<double,dim,1> b1(BarycentricTraits<dim-1>::x2l(a1));
        const Eigen::Matrix<double,dim+1,1> bary(face2domainBary(b1,boundaryFace));
        const VectorDim pos(ele.position(bary));
//        std::cout<<pos.transpose()<<std::endl;
        
        Eigen::Matrix<double,dim,dim> stress(Eigen::Matrix<double,dim,dim>::Zero());
        stress(2,2)=std::pow(pos(0),3);
        return stress*JGNselector<dim>::jGN(ele.jGN(bary,boundaryFace));
    }

    
};


int main (int argc, char * const argv[])
{
    
    int meshID(0);
    if (argc>1)
    {
        meshID=atoi(argv[1]);
    }

    
    // Create a FiniteElement object on the mesh
    
    
    TestIntegrator ti(meshID);

    
    return 0;
}


