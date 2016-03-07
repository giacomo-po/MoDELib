
#ifndef _Torsioner_h_
#define _Torsioner_h_

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <model/FEM/Domains/IntegrationDomain.h>
#include <model/Utilities/EigenDataReader.h>
#include <model/FEM/Boundaries/AtXmin.h>
#include <model/FEM/Boundaries/AtXmax.h>
#include <model/FEM/Boundaries/AtPoint.h>
#include <model/FEM/BoundaryConditions/Fix.h>


using namespace model;

struct Torsioner
{
    
    //! z-coordinate of face at which Torsioner is applied
    const double zFace;
    
    //! Height of sample in z direction
    const double Lz;
    
    const Eigen::Matrix<double,3,1> pivot;

    //! Rate of rotation
    static double thetaDot;

    //! Tolerance used to detect nodes on zFace
    const double geometricTol;
    
    //! Initial value of rotation
    static double initialTwist_Rad;
    
    // Number of DD steps before applying strain.
    static int relaxSteps;
    
    static double last_update_time;
    
    static double deltaTheta;

    static bool apply_torsion;

    
    /**************************************/
    template <typename FiniteElementType>
    Torsioner(const FiniteElementType& fe) :
    /* init list */ zFace(fe.xMax()(2)), // zFace is the one at x(2)= max among finite element node coordinates
    /* init list */ Lz(fe.xMax()(2)-fe.xMin()(2)),
    /* init list */ pivot((fe.xMax()(0)+fe.xMin()(0))*0.5,(fe.xMax()(1)+fe.xMin()(1))*0.5,fe.xMax()(2)),
    /* init list */ geometricTol(0.01)
    {
    }
    
    /**************************************/
    template <typename NodeType>
    bool operator()(const NodeType& node) const
    {/*!@param[in] node a finite element node
      *\returns true if the node coordinate is under the Torsioner
      */

        return (std::fabs(node.P0(2)-zFace)<geometricTol);
        
    }
    
    /**************************************/
    template <typename NodeType,typename TrialFunctionType>
    Eigen::Matrix<double,3,1> operator()(const NodeType& node,const TrialFunctionType&) const
    {/*!@param[in] node the FiniteElement node to be displaced
      *\returns the displacement of node
      */
        Eigen::Matrix<double,3,3> rot(Eigen::AngleAxisd(initialTwist_Rad+deltaTheta, Eigen::Vector3d::UnitZ()));
        return (Eigen::Matrix<double,3,1>()<<(rot*(node.P0-pivot)-(node.P0-pivot)).template segment<2>(0),0.0).finished();
    }
    
    /**************************************/
    template <typename DislocationNetworkType>
    void updateDisplacement(const DislocationNetworkType& DN)
    {/*!@param[in] DN the DislocationNetwork
      * 
      * Computes the change in rotation angle based on the time allapsed
      * since last update.
      */
        if(DN.runningID()>=relaxSteps)
        {
            const double deltaT = DN.get_totalTime() - last_update_time;
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
                for (unsigned int f=0;f<boundaryFaces.size();++f) // loop ever bonudary faces of the current Elements
                {
                    bool isExternalBoundaryFace(true);
                    std::array<const Simplex<FiniteElementType::dim,0>*, SimplexTraits<FiniteElementType::dim,FiniteElementType::dim-1>::nVertices> vertices=eIter->second.simplex.child(boundaryFaces[f]).vertices();
                    for(unsigned int v=0;v<vertices.size();++v) // loop over vertices of the current face
                    {
                        isExternalBoundaryFace *= Torsioner(fe)(*vertices[v]); // check if the current vertices satisfies operator()
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
    
    /**************************************/
    static void init(const long int& runID,
                     const unsigned int& userOutputColumn)
    {
        
        EigenDataReader EDR;
        
        EDR.readScalarInFile("./loadInput.txt","thetaDot",thetaDot);
        EDR.readScalarInFile("./loadInput.txt","initialTwist_Rad",initialTwist_Rad);
        EDR.readScalarInFile("./loadInput.txt","relaxSteps",relaxSteps);
        EDR.readScalarInFile("./loadInput.txt","apply_torsion",apply_torsion);

        
        
        VertexReader<'F',201,double> vReader;
        if (vReader.isGood(0,true))
        {
            vReader.read(0,true);
            const auto iter=vReader.find(runID);
            if (iter!=vReader.end())
            {
                model::cout<<"Initializing Torsioner at runID="<<runID<<std::endl;

                initialTwist_Rad=iter->second(userOutputColumn-1);
                model::cout<<"initialTwist_Rad="<<initialTwist_Rad<<std::endl;
                
                last_update_time=iter->second(0);
                model::cout<<"last_update_time="<<last_update_time<<std::endl;
                
            }
            else
            {
                //                assert(0 && "Torsioner::init runID not found inf F file");
            }
        }
        else
        {
            model::cout<<"Torsioner F/F_0.txt cannot be opened."<<std::endl;
        }
        
    }
    
};

// Static data
double Torsioner::thetaDot=1.0e-12;
double Torsioner::initialTwist_Rad=0.0;
int Torsioner::relaxSteps=0;
double Torsioner::last_update_time=0.0;
double Torsioner::deltaTheta=0.0;
bool Torsioner::apply_torsion=true;

#endif

