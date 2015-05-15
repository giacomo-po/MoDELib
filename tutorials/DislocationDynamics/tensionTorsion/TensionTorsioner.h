
#ifndef _TensionTorsioner_h_
#define _TensionTorsioner_h_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <model/FEM/Domains/IntegrationDomain.h>

using namespace model;

//template <typename DislocationNetworkType>
struct TensionTorsioner
{
    
    //! z-coordinate of face at which TensionTorsioner is applied
    const double zFace;
    
    //! Height of sample in z direction
    const double Lz;

    const double R;
    
    const Eigen::Matrix<double,3,1> pivot;
    
    const double thetaDot;
    
    //! Height of sample in z direction
    const double geometricTol;
    
    //! Applied axial strain rate
    static double epsilonDot;

    static double gammaDot;

    
    //! Initial value of the punch displacement
    static double initialDisplacement;
    static double initialTwist_Rad;
    
    // Number of DD steps before applying strain.
    static int relaxSteps;
    
    static double last_update_time;
    
    static double deltaDisplacement;
    static double deltaTheta;

    static bool apply_tension;
    static bool apply_torsion;
//    static double target_VonMisesStrain;

    
    /**************************************/
    template <typename FiniteElementType>
    TensionTorsioner(const FiniteElementType& fe) :
    //    /* init list */ DN(DN_in),
    /* init list */ zFace(fe.xMax()(2)), // zFace is the one at x(2)= max among finite element node coordinates
    /* init list */ Lz(fe.xMax()(2)-fe.xMin()(2)),
    /* init list */ R((fe.xMax()(0)-fe.xMin()(0))*0.5),
    /* init list */ pivot((fe.xMax()(0)+fe.xMin()(0))*0.5,(fe.xMax()(1)+fe.xMin()(1))*0.5,fe.xMax()(2)),
    /* init list */ thetaDot(Lz/R*gammaDot),
    /* init list */ geometricTol(0.01)
    {
    }
    
    /**************************************/
    template <typename NodeType>
    bool operator()(const NodeType& node) const
    {/*!@param[in] node a finite element node
      *\returns true if the node coordinate is under the TensionTorsioner
      */

        return (std::fabs(node.P0(2)-zFace)<geometricTol);
        
    }
    
    /**************************************/
    template <typename NodeType,typename TrialFunctionType>
    Eigen::Matrix<double,3,1> operator()(const NodeType& node,const TrialFunctionType&) const
    {
        Eigen::Matrix<double,3,3> rot(Eigen::AngleAxisd(initialTwist_Rad+deltaTheta, Eigen::Vector3d::UnitZ()));
        return (Eigen::Matrix<double,3,1>()<<(rot*(node.P0-pivot)-(node.P0-pivot)).template segment<2>(0),initialDisplacement+deltaDisplacement).finished();
    }
    
    /**************************************/
    template <typename DislocationNetworkType>
    void updateDisplacement(const DislocationNetworkType& DN)
    {
        if(DN.runningID()>=relaxSteps)
        {
            const double deltaT = DN.get_totalTime() - last_update_time;
            deltaDisplacement += epsilonDot*deltaT*Lz;
            deltaTheta += thetaDot*deltaT;
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
                        isExternalBoundaryFace *= TensionTorsioner(fe)(*vertices[v]); // check if the current vertices satisfies operator()
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
double TensionTorsioner::epsilonDot=1.0e-9;
double TensionTorsioner::gammaDot=1.0e-9;
double TensionTorsioner::initialDisplacement=0.0; // 0.0;
double TensionTorsioner::initialTwist_Rad=0.0;
int TensionTorsioner::relaxSteps=0;
double TensionTorsioner::last_update_time=0.0;
double TensionTorsioner::deltaDisplacement=0.0;
double TensionTorsioner::deltaTheta=0.0;
bool TensionTorsioner::apply_tension=true;
bool TensionTorsioner::apply_torsion=true;
//double TensionTorsioner::target_VonMisesStrain=0.01;

#endif

