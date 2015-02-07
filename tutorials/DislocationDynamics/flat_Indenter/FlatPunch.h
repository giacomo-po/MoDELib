
#ifndef _FlatPunch_h_
#define _FlatPunch_h_

#include <model/FEM/Domains/IntegrationDomain.h>

using namespace model;

//template <typename DislocationNetworkType>
struct FlatPunch
{
    //! const reference to the DislocationNetwork
    //    const DislocationNetworkType& DN;
    
    //! z-coordinate of face at which FlatPunch is applied
    const double zFace;
    
    //! Height of sample in z direction
    const double Lz;
    
    static double centerX;
    
    static double centerY;
    
    //! Height of sample in z direction
    const double geometricTol;
    
    //! side length of the FlatPunch
    static double punchSize;
    
    //! Applied strain rate
    static double strainRate;
    
    //! Initial value of the punch displacement
    static double initialDisplacement;
    
    // Number of DD steps before applying strain.
    static int relaxSteps;
    
    static double last_update_time;
    
    static double deltaDisplacement;
    
    /**************************************/
    template <typename FiniteElementType>
    FlatPunch(const FiniteElementType& fe) :
    //    /* init list */ DN(DN_in),
    /* init list */ zFace(fe.xMax()(2)), // zFace is the one at x(2)= max among finite element node coordinates
    /* init list */ Lz(fe.xMax()(2)-fe.xMin()(2)),
    /* init list */ geometricTol(0.01)
    {
    }
    
    /**************************************/
    template <typename NodeType>
    bool operator()(const NodeType& node) const
    {/*!@param[in] node a finite element node
      *\returns true if the node coordinate is under the FlatPunch
      */

        return (   std::fabs(node.P0(0)-centerX)<(0.5*punchSize+geometricTol)
                /*   */ && std::fabs(node.P0(1)-centerY)<(0.5*punchSize+geometricTol)
                /*   */ && std::fabs(node.P0(2)-zFace)<geometricTol);
        
    }
    
    /**************************************/
    template <typename NodeType,typename TrialFunctionType>
    Eigen::Matrix<double,3,1> operator()(const NodeType&,const TrialFunctionType&) const
    {
        return (Eigen::Matrix<double,3,1>()<<0.0,0.0,initialDisplacement+deltaDisplacement).finished();
    }
    
    /**************************************/
    template <typename DislocationNetworkType>
    void updateDisplacement(const DislocationNetworkType& DN)
    {
        if(DN.runningID()>=relaxSteps)
        {
            const double deltaT = DN.get_totalTime() - last_update_time;
            deltaDisplacement += strainRate*deltaT*Lz;
            last_update_time += deltaT;
        }
    }
    
    /**************************************/
    template <typename FiniteElementType>
    static IntegrationDomain<FiniteElementType,1,3,GaussLegendre> boundary(const FiniteElementType& fe)
    {
        IntegrationDomain<FiniteElementType,1,3,GaussLegendre> temp;
        
        // loop ever Elements
        for (typename FiniteElementType::ElementContainerType::const_iterator eIter =fe.elementBegin();
             /*                                                            */ eIter!=fe.elementEnd();
             /*                                                            */ eIter++)
        {
            if(eIter->second.isBoundaryElement())
            {
                const std::vector<int> boundaryFaces=eIter->second.boundaryFaces();
                for (int f=0;f<boundaryFaces.size();++f) // loop ever bonudary faces of the current Elements
                {
                    bool isExternalBoundaryFace(true);
                    std::array<const Simplex<FiniteElementType::dim,0>*, SimplexTraits<FiniteElementType::dim,FiniteElementType::dim-1>::nVertices> vertices=eIter->second.simplex.child(boundaryFaces[f]).vertices();
                    for(int v=0;v<vertices.size();++v) // loop over vertices of the current face
                    {
                        isExternalBoundaryFace *= FlatPunch(fe)(*vertices[v]); // check if the current vertices satisfies operator()
                    }
                    if(isExternalBoundaryFace)
                    {
                        temp.emplace_back(&eIter->second,boundaryFaces[f]);
                    }
                }
            }
        }
        
        return temp;
    }
    
};

// Static data
double FlatPunch::punchSize=300.0;
double FlatPunch::centerX=0.0;
double FlatPunch::centerY=0.0;
double FlatPunch::strainRate=-1.0e-9;
double FlatPunch::initialDisplacement=0; // 0.0;
int FlatPunch::relaxSteps=0;
double FlatPunch::last_update_time=0.0;
double FlatPunch::deltaDisplacement=0.0;

#endif

