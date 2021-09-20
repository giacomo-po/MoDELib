/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureGenerator_cpp_
#define model_MicrostructureGenerator_cpp_

#include <MicrostructureGenerator.h>


namespace model
{

    
      
        PeriodicPlanePatch<3>* PolyPoint::periodicPlanePatch() const
        {
            return nullptr;
        }
        
    
        
        /**********************************************************************/
        double MicrostructureGenerator::min(const double& a,const double& b)
        {
            return a<b? a : b;
        }

        /**********************************************************************/
        double MicrostructureGenerator::max(const double& a,const double& b)
        {
            return a>b? a : b;
        }
        
        /**********************************************************************/
        typename MicrostructureGenerator::VectorDimD MicrostructureGenerator::randomOrthogonalUnitVector(VectorDimD v)
        {
            const double vNorm(v.norm());
            assert(vNorm>FLT_EPSILON);
            v/=vNorm;
            
            VectorDimD temp(v.cross(VectorDimD::Random()));
            double tempNorm(temp.norm());
            while(tempNorm<FLT_EPSILON)
            {
                temp=v.cross(VectorDimD::Random());
                tempNorm=temp.norm();
            }
            return temp/tempNorm;
        }
        
        /**********************************************************************/
        bool MicrostructureGenerator::addSingleLoop(const bool randomizeBurgersSense,
                           const std::vector<VectorDimD>& nodePos,
                           const std::vector<VectorDimD>& loopNodePos,
                           VectorDimD b,
                           const VectorDimD& unitNormal,
                           const VectorDimD& P0,
                           const int& grainID,
                           const int& loopType,
                           const std::vector<VectorDimD>& loopNodeShifts,
                           const std::vector<short int>& periodicEdgeIDs // default -1 for non-periodic nodes
                           )
//                           const long int& periodicLoopID,
//                           const VectorDimD& periodicShift)
        {
            
            assert(nodePos.size()==loopNodePos.size());
            
            if(allPointsInGrain(nodePos,grainID))
            {
                if(randomizeBurgersSense)
                {
                    b*=randomSign();
                }
                
                for(size_t k=0;k<nodePos.size();++k)
                {
//                    const size_t nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
                    const size_t nextLoopNodeID=(k+1)<nodePos.size()? loopNodeID+k+1 : loopNodeID;
                    const size_t k1=(k+1)<nodePos.size()? k+1 : 0;
                    configIO.nodes().emplace_back(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                    configIO.loopNodes().emplace_back(loopNodeID+k,loopID,loopNodePos[k],nodeID+k,loopNodeShifts[k],periodicEdgeIDs[k]);
                    configIO.loopLinks().emplace_back(loopID,loopNodeID+k,nextLoopNodeID,(loopNodePos[k]-loopNodePos[k1]).norm()>FLT_EPSILON,0);
                }
                configIO.loops().emplace_back(loopID, b,unitNormal,P0,grainID,loopType);
                nodeID+=nodePos.size();
                loopNodeID+=loopNodePos.size();
                loopID+=1;
//                snID+=1;
                return true;
            }
            else
            {
                return false;
            }
        }

        /**********************************************************************/
        bool MicrostructureGenerator::addSingleLoopwithJunction(const bool randomizeBurgersSense,
                           const std::vector<VectorDimD>& nodePos,
                           const std::vector<VectorDimD>& loopNodePos,
                           VectorDimD b,
                           const VectorDimD& unitNormal,
                           const VectorDimD& P0,
                           const int& grainID,
                           const int& loopType,
                           const std::vector<VectorDimD>& loopNodeShifts,
                           const std::vector<short int>& periodicEdgeIDs // default -1 for non-periodic nodes
                           )
//                           const long int& periodicLoopID,
//                           const VectorDimD& periodicShift)
        {
            
            assert(nodePos.size()==loopNodePos.size());
            
            if(allPointsInGrain(nodePos,grainID))
            {
                if(randomizeBurgersSense)
                {
                    b*=randomSign();
                }
                size_t newNodeAdded(0);
                size_t oldNodeUsed(0);

                size_t newNodeSize(nodePos.size());

                for(size_t k=0;k<nodePos.size();++k)
                {
//                    const size_t nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
                    const size_t nextLoopNodeID=(k+1)<nodePos.size()? loopNodeID+k+1 : loopNodeID;
                    const size_t k1=(k+1)<nodePos.size()? k+1 : 0;
                    long int existinNodeID(-1);
                    for (const auto& netNode : configIO.nodes())
                    {
                        if ((netNode.P-nodePos[k]).squaredNorm()<FLT_EPSILON)
                        {
                            existinNodeID=netNode.sID;
                            break;
                        }
                    }
                    if (existinNodeID>=0)
                    {
                        configIO.loopNodes().emplace_back(loopNodeID + k, loopID, loopNodePos[k], existinNodeID, loopNodeShifts[k], periodicEdgeIDs[k]);
                    oldNodeUsed++;
                    }
                    else
                    {
                        configIO.nodes().emplace_back(nodeID + k-oldNodeUsed, nodePos[k], Eigen::Matrix<double, 1, 3>::Zero(), 1.0, 0);
                        configIO.loopNodes().emplace_back(loopNodeID + k, loopID, loopNodePos[k], nodeID + k-oldNodeUsed, loopNodeShifts[k], periodicEdgeIDs[k]);
                        newNodeAdded++;
                    }

                    configIO.loopLinks().emplace_back(loopID,loopNodeID+k,nextLoopNodeID,(loopNodePos[k]-loopNodePos[k1]).norm()>FLT_EPSILON,0);
                }
                configIO.loops().emplace_back(loopID, b,unitNormal,P0,grainID,loopType);
                nodeID+=newNodeAdded;
                loopNodeID+=loopNodePos.size();
                loopID+=1;
//                snID+=1;
                return true;
            }
            else
            {
                return false;
            }
        }

        /**********************************************************************/
        void MicrostructureGenerator::addStraightDislocations()
        {
            if(targetStraightDislocationDensity>0.0)
            {
                // init counters
                double density=0.0;
                double sessileDensity=0.0;

                std::cout<<magentaBoldColor<<"Generating straight dislocations"<<defaultColor<<std::endl;
                while(density<targetStraightDislocationDensity)
                {
                    
                    const std::pair<LatticeVector<MicrostructureGenerator::dim>,int> rp=randomPointInMesh();
                    const int& grainID=rp.second;   // random grain ID
                    const LatticeVector<MicrostructureGenerator::dim>& L0=rp.first; // random lattice position in the grain
                    const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                    const int rSS=distribution(generator); // a random SlipSystem ID
                    const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    const VectorDimD b=slipSystem.s.cartesian();    // Burgers vector


                    VectorDimD n=slipSystem.unitNormal; // slip plane normal
                    std::uniform_real_distribution<> dis(0.0, 2.0*M_PI);
                    const double theta=dis(generator); // random angle of the dislocation line in the plane from screw orientation.
                    VectorDimD d=Eigen::AngleAxisd(theta, n)*b.normalized();

                    bool isSessile=false;
                    if(sessileDensity/targetStraightDislocationDensity<fractionSessileStraightDislocationDensity)
                    {
                        n=b.normalized();
                        isSessile=true;
                        d=Eigen::AngleAxisd(theta, n)*n.cross(VectorDimD::Random()).normalized();
                    }

                    MeshPlane<MicrostructureGenerator::dim> plane(mesh,grainID,P0,n);
                    const std::vector<VectorDimD> nodePos(DislocationInjector<MicrostructureGenerator::dim>::straightLineBoundaryClosure(P0,d,plane,mesh));

                    const double lineLength=(nodePos.back()-nodePos.front()).norm();
                    //                nodePos.push_back(nodePos[nodePos.size()-1]+1.0/3.0*(nodePos[0]-nodePos[nodePos.size()-1]));
                    //                nodePos.push_back(nodePos[nodePos.size()-1]+2.0/3.0*(nodePos[0]-nodePos[nodePos.size()-1]));

                    double dh=0.0;
                    if(enforceMonotonicHelicity)
                    {
                        dh=deltaHelicity(nodePos,b);
                    }

                    // Write files
                    if(   nodePos.size()>=3
                       && ((fabs(helicity+dh)>fabs(helicity) && fabs(dh)>FLT_EPSILON ) || helicity==0.0 || !enforceMonotonicHelicity))
                    {
                        // write node and edge file
//                        for(size_t k=0;k<nodePos.size();++k)
//                        {
//                            configIO.nodes().emplace_back(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,0);
//                            const int nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
//                            configIO.loopLinks().emplace_back(loopID,nodeID+k,nextNodeID,0);
//                        }
//                        nodeID+=nodePos.size();

                        // write loop file
//                        configIO.loops().emplace_back(loopID+0, b,n,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP);

                        addSingleLoop(false,nodePos,nodePos, b,n,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP,std::vector<VectorDimD>(nodePos.size(),VectorDimD::Zero()),std::vector<short int>(nodePos.size(),-1));
                        
                        if(enforceMonotonicHelicity)
                        {
                            loopPoints.push_back(nodePos);
                            loopBurgers.push_back(b);
                        }

//                        loopID+=1;
//                        snID+=1;
                        density += lineLength/mesh.volume()/std::pow(poly.b_SI,2);
                        if(isSessile)
                        {
                            sessileDensity += lineLength/mesh.volume()/std::pow(poly.b_SI,2);
                        }
                        std::cout<<"theta="<<theta*180.0/M_PI<<", density="<<density<<" (sessileDensity="<<sessileDensity<<")"<<std::endl;
                        if(enforceMonotonicHelicity)
                        {
                            loopPoints.push_back(nodePos);
                            loopBurgers.push_back(b);
                            helicity+=dh;
                            std::cout<<"helicity="<<helicity<<std::endl;
                        }
                        
                    }
                }
            }
        }

        /**********************************************************************/
        void MicrostructureGenerator::addFrankReadSources()
        {
            if(targetFrankReadDislocationDensity>0.0)
            {
                std::cout<<magentaBoldColor<<"Generating Frank-Read sources"<<defaultColor<<std::endl;

                double density=0.0;
                double edgeDensity=0.0;

                const double fractionEdge=1.0; // TEMPORARY

                std::normal_distribution<double> sizeDistribution(FrankReadSizeMean/poly.b_SI,FrankReadSizeStd/poly.b_SI);
                std::normal_distribution<double> aspectRatioDistribution(FrankReadAspectRatioMean,FrankReadAspectRatioStd);

                while(density<targetFrankReadDislocationDensity)
                {
                    const std::pair<LatticeVector<MicrostructureGenerator::dim>,int> rp=randomPointInMesh();
                    const LatticeVector<MicrostructureGenerator::dim> L0=rp.first;
                    const int grainID=rp.second;

                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);

                    const int rSS=distribution(generator); // a random SlipSystem

                    const auto& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    const VectorDimD b(slipSystem.s.cartesian());

                    // Compute the ReciprocalLatticeDirection corresponding to s
                    ReciprocalLatticeDirection<3> sr(poly.grain(grainID).reciprocalLatticeDirection(b));

                    bool isEdge=true;



                    LatticeDirection<3> d1(LatticeVector<MicrostructureGenerator::dim>(sr.cross(slipSystem.n)*randomSign()));
                    double d1cNorm(d1.cartesian().norm());
                    //                    const double size = distribution(generator)*inclusionsDiameterLognormalDistribution_A[f]/poly.b_SI;

                    int a1=sizeDistribution(generator)/d1cNorm;
                    if(a1>0)
                    {
                        LatticeVector<MicrostructureGenerator::dim> L1=L0+d1*a1;


                        if(edgeDensity>fractionEdge*density) // overwrite with screw dislocaiton
                        {
                            isEdge=false;
                            d1cNorm=b.norm();
                            a1=randomSize()/d1cNorm;
                            L1=L0+slipSystem.s.dir*a1;
                        }

                        // Compute the LatticeDireciton corresponding to -n
                        LatticeDirection<3> d2(poly.grain(grainID).latticeDirection(-slipSystem.n.cartesian()*randomSign()));
                        double d2cNorm(d2.cartesian().norm());

                        const int a2=aspectRatioDistribution(generator)*a1; // aspect ratio of double FR source
                        if(a2>0)
                        {
                            LatticeVector<MicrostructureGenerator::dim> L2=L1+d2*a2;
                            LatticeVector<MicrostructureGenerator::dim> L3=L0+d2*a2;

                            const VectorDimD P0=L0.cartesian();
                            const VectorDimD P1=L1.cartesian();
                            const VectorDimD P2=L2.cartesian();
                            const VectorDimD P3=L3.cartesian();
                            const auto search1(mesh.search(P1));
                            const auto search2(mesh.search(P2));
                            const auto search3(mesh.search(P3));

                            double dh=0.0;
                            std::vector<VectorDimD> nodePos;
                            if(enforceMonotonicHelicity)
                            {
                                nodePos.push_back(P0);
                                nodePos.push_back(P1);
                                nodePos.push_back(P2);
                                nodePos.push_back(P3);
                                dh=deltaHelicity(nodePos,-b); // central loop is opposite direction
                            }

                            if(   search1.first && search1.second->region->regionID==grainID
                               && search2.first && search2.second->region->regionID==grainID
                               && search3.first && search3.second->region->regionID==grainID
                               && ((fabs(helicity+dh)>fabs(helicity) && fabs(dh)>FLT_EPSILON ) || helicity==0.0 || !enforceMonotonicHelicity)
                               )
                            {
                                density += 2.0*(d1cNorm*a1 + d2cNorm*a2)/mesh.volume()/std::pow(poly.b_SI,2);
                                if(isEdge)
                                {
                                    edgeDensity+=2.0*(d1cNorm*a1 + d2cNorm*a2)/mesh.volume()/std::pow(poly.b_SI,2);
                                }
                                std::cout<<"density="<<density<<" (edge density="<<edgeDensity<<")"<<std::endl;

                                const VectorDimD P4=0.5*(P0+P1);
                                const VectorDimD P5=0.5*(P2+P3);

                                const VectorDimD n1=slipSystem.unitNormal;
                                const VectorDimD n2=d2.cross(d1).cartesian().normalized();

                                configIO.nodes().emplace_back(nodeID+0,P0,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                                configIO.nodes().emplace_back(nodeID+1,P1,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                                configIO.nodes().emplace_back(nodeID+2,P2,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                                configIO.nodes().emplace_back(nodeID+3,P3,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                                configIO.nodes().emplace_back(nodeID+4,P4,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                                configIO.nodes().emplace_back(nodeID+5,P5,Eigen::Matrix<double,1,3>::Zero(),1.0,0);

                               
                               
                                configIO.loopNodes().emplace_back(loopNodeID+0,loopID+0,P0,nodeID+0,VectorDimD::Zero(),-1);
                                configIO.loopNodes().emplace_back(loopNodeID+1,loopID+0,P1,nodeID+1,VectorDimD::Zero(),-1);
                                configIO.loopNodes().emplace_back(loopNodeID+2,loopID+0,P4,nodeID+4,VectorDimD::Zero(),-1);

                                configIO.loopNodes().emplace_back(loopNodeID+3,loopID+1,P0,nodeID+0,VectorDimD::Zero(),-1);
                                configIO.loopNodes().emplace_back(loopNodeID+4,loopID+1,P3,nodeID+3,VectorDimD::Zero(),-1);
                                configIO.loopNodes().emplace_back(loopNodeID+5,loopID+1,P2,nodeID+2,VectorDimD::Zero(),-1);
                                configIO.loopNodes().emplace_back(loopNodeID+6,loopID+1,P1,nodeID+1,VectorDimD::Zero(),-1);

                                configIO.loopNodes().emplace_back(loopNodeID+7,loopID+2,P3,nodeID+3,VectorDimD::Zero(),-1);
                                configIO.loopNodes().emplace_back(loopNodeID+8,loopID+2,P5,nodeID+5,VectorDimD::Zero(),-1);
                                configIO.loopNodes().emplace_back(loopNodeID+9,loopID+2,P2,nodeID+2,VectorDimD::Zero(),-1);



                                configIO.loops().emplace_back(loopID+0,b,n1,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP);
                                configIO.loops().emplace_back(loopID+1,b,n2,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::SESSILELOOP);
                                configIO.loops().emplace_back(loopID+2,b,n1,P3,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP);

                                // configIO.loopLinks().emplace_back(loopID+0,nodeID+0,nodeID+1,true,0);
                                // configIO.loopLinks().emplace_back(loopID+0,nodeID+1,nodeID+4,true,0);
                                // configIO.loopLinks().emplace_back(loopID+0,nodeID+4,nodeID+0,true,0);

                                // configIO.loopLinks().emplace_back(loopID+1,nodeID+0,nodeID+3,true,0);
                                // configIO.loopLinks().emplace_back(loopID+1,nodeID+3,nodeID+2,true,0);
                                // configIO.loopLinks().emplace_back(loopID+1,nodeID+2,nodeID+1,true,0);
                                // configIO.loopLinks().emplace_back(loopID+1,nodeID+1,nodeID+0,true,0);

                                // configIO.loopLinks().emplace_back(loopID+2,nodeID+3,nodeID+5,true,0);
                                // configIO.loopLinks().emplace_back(loopID+2,nodeID+5,nodeID+2,true,0);
                                // configIO.loopLinks().emplace_back(loopID+2,nodeID+2,nodeID+3,true,0);


                                configIO.loopLinks().emplace_back(loopID+0,loopNodeID+0,loopNodeID+1,true,0);
                                configIO.loopLinks().emplace_back(loopID+0,loopNodeID+1,loopNodeID+2,true,0);
                                configIO.loopLinks().emplace_back(loopID+0,loopNodeID+2,loopNodeID+0,true,0);

                                configIO.loopLinks().emplace_back(loopID+1,loopNodeID+3,loopNodeID+4,true,0);
                                configIO.loopLinks().emplace_back(loopID+1,loopNodeID+4,loopNodeID+5,true,0);
                                configIO.loopLinks().emplace_back(loopID+1,loopNodeID+5,loopNodeID+6,true,0);
                                configIO.loopLinks().emplace_back(loopID+1,loopNodeID+6,loopNodeID+3,true,0);

                                configIO.loopLinks().emplace_back(loopID+2,loopNodeID+7,loopNodeID+8,true,0);
                                configIO.loopLinks().emplace_back(loopID+2,loopNodeID+8,loopNodeID+9,true,0);
                                configIO.loopLinks().emplace_back(loopID+2,loopNodeID+9,loopNodeID+7,true,0);

                                nodeID+=6;
                                loopNodeID+=10;
                                loopID+=3;
//                                snID++;

                                if(enforceMonotonicHelicity)
                                {
                                    loopPoints.push_back(nodePos);
                                    loopBurgers.push_back(-b);
                                    helicity+=dh;
                                    std::cout<<"helicity="<<helicity<<std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }

        /**********************************************************************/
        void MicrostructureGenerator::addSingleArmDislocations()
        {
            if(targetSingleArmDislocationDensity>0.0)
            {
                std::cout<<magentaBoldColor<<"Generating single-arm sources"<<defaultColor<<std::endl;

                //                double fractionEdge=1.0; // TEMPORARY

                double density=0.0;
                //                double edgeDensity=0.0;

                while(density<targetSingleArmDislocationDensity)
                {
                    const std::pair<LatticeVector<MicrostructureGenerator::dim>,int> rp=randomPointInMesh();
                    const LatticeVector<MicrostructureGenerator::dim> L0=rp.first;
                    const int grainID=rp.second;

                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);

                    const int rSS=distribution(generator); // a random SlipSystem

                    const auto& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);

                    // Compute the ReciprocalLatticeDirection corresponding to s
                    ReciprocalLatticeDirection<3> sr(poly.grain(grainID).reciprocalLatticeDirection(slipSystem.s.cartesian()));

                    //                    bool isEdge=true;

                    LatticeDirection<3> d1(LatticeVector<MicrostructureGenerator::dim>(sr.cross(slipSystem.n)*randomSign()));  // this is the projection direction
                    //     if(edgeDensity>fractionEdge*density) // overwrite with screw dislocaiton
                    //      {
                    //           d1=slipSystem.s;
                    //      }
                    // Compute the LatticeDireciton corresponding to -n
                    LatticeDirection<3> d2(poly.grain(grainID).latticeDirection(-slipSystem.n.cartesian()*randomSign()));
                    double d2cNorm(d2.cartesian().norm());
                    const int a2=randomSize()/d2cNorm;
                    LatticeVector<MicrostructureGenerator::dim> L3=L0+d2*a2;


                    const auto search2(mesh.search(L3.cartesian()));

                    if(enforceMonotonicHelicity)
                    {
                        assert(false && "enforceMonotonicHelicity NOT YET IMPLEMENTED FOR SINGLE-ARM SOURCES");
                    }

                    if(  search2.first && search2.second->region->regionID==grainID)
                    {
                        const VectorDimD P0=L0.cartesian();
                        const VectorDimD P3=L3.cartesian();


                        const VectorDimD n1=slipSystem.unitNormal;
                        const VectorDimD n2=d2.cross(d1).cartesian().normalized();

                        MeshPlane<MicrostructureGenerator::dim> plane01(mesh,grainID,P0,n2);
//                        PlaneMeshIntersectionContainerType pmi01=PlaneMeshIntersection<MicrostructureGenerator::dim>(mesh,P0,n2,grainID);
                        const VectorDimD P1=boundaryProjection(P0,d1.cartesian(),plane01.meshIntersections).second;
                        const VectorDimD P2=boundaryProjection(P3,d1.cartesian(),plane01.meshIntersections).second;
                        const std::map<double,VectorDimD> P12=boundaryProjection(P0,P3,d1.cartesian(),plane01.meshIntersections);

                        // const VectorDimD P1=P12.begin().second;
                        // const VectorDimD P2=P12.rbegin().second;
                        const VectorDimD P4=(P0+P1)/2.0;
                        const VectorDimD P5=(P3+P2)/2.0;

                        if ((P1-P0).norm()>a2*0.5 && (P3-P2).norm()>a2*0.5)  //not too small arm
                        {
                            density+=((P1-P0).norm()+(P0-P3).norm()+(P3-P2).norm())/mesh.volume()/std::pow(poly.b_SI,2);
                            /*! Vertex file format is:
                             * ID Px Py Pz Vx Vy Vz velReducCoeff snID meshLocation grainID
                             */
                            const size_t refNodeID=nodeID;
                            configIO.nodes().emplace_back(refNodeID+0,P0,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                            configIO.nodes().emplace_back(refNodeID+1,P3,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                            configIO.nodes().emplace_back(refNodeID+2,P4,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                            configIO.nodes().emplace_back(refNodeID+3,P5,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                            configIO.nodes().emplace_back(refNodeID+4,P1,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                            configIO.nodes().emplace_back(refNodeID+5,P2,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                            nodeID+=6;

                            if (P12.size()==0)
                            {
                                configIO.loopLinks().emplace_back(loopID,refNodeID+4,refNodeID+5,true,0);
                            }
                            else if (P12.size()==1)
                            {
                                configIO.nodes().emplace_back(nodeID,P12.begin()->second,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                                configIO.loopLinks().emplace_back(loopID,refNodeID+4,nodeID,true,0);
                                configIO.loopLinks().emplace_back(loopID,nodeID,refNodeID+5,true,0);
                                nodeID++;
                            }
                            else
                            {
                                for(const auto pair : P12)
                                {
                                    configIO.nodes().emplace_back(nodeID,pair.second,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                                    nodeID++;
                                    if (pair.first==P12.begin()->first)
                                    {
                                        configIO.loopLinks().emplace_back(loopID,refNodeID+4,nodeID,true,0);
                                        configIO.loopLinks().emplace_back(loopID,nodeID,nodeID+1,true,0);
                                    }
                                    else if (pair.first==P12.rbegin()->first )
                                    {
                                        configIO.loopLinks().emplace_back(loopID,nodeID,refNodeID+5,true,0);
                                    }
                                    else
                                    {
                                        configIO.loopLinks().emplace_back(loopID,nodeID,nodeID+1,true,0);
                                    }
                                }
                            }
                            configIO.loops().emplace_back(loopID,slipSystem.s.cartesian(),n2,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP);
                            configIO.loopLinks().emplace_back(loopID,refNodeID+5,refNodeID+1,true,0);
                            configIO.loopLinks().emplace_back(loopID,refNodeID+1,refNodeID+0,true,0);
                            configIO.loopLinks().emplace_back(loopID,refNodeID+0,refNodeID+4,true,0);

                            configIO.loops().emplace_back(loopID+1,slipSystem.s.cartesian(),n1,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::SESSILELOOP);
                            configIO.loopLinks().emplace_back(loopID+1,refNodeID+0,refNodeID+2,true,0);
                            configIO.loopLinks().emplace_back(loopID+1,refNodeID+2,refNodeID+4,true,0);
                            configIO.loopLinks().emplace_back(loopID+1,refNodeID+4,refNodeID+0,true,0);

                            configIO.loops().emplace_back(loopID+2,slipSystem.s.cartesian(),n1,P3,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP);
                            configIO.loopLinks().emplace_back(loopID+2,refNodeID+1,refNodeID+5,true,0);
                            configIO.loopLinks().emplace_back(loopID+2,refNodeID+5,refNodeID+3,true,0);
                            configIO.loopLinks().emplace_back(loopID+2,refNodeID+3,refNodeID+1,true,0);

                            loopID+=3;
//                            snID++;

                            std::cout<<"density="<<density<<std::endl;
                        }
                    }
                }
            }
        }

       
        /**********************************************************************/
        void MicrostructureGenerator::addPrismaticLoops()
        {
            if (targetPrismaticLoopDensity > 0.0)
            {
                std::cout << magentaBoldColor << "Generating prismatic loops" << defaultColor << std::endl;

                double density = 0.0;

                while (density < targetPrismaticLoopDensity)
                {

                    const std::pair<LatticeVector<MicrostructureGenerator::dim>, int> rp = randomPointInMesh();
                    const LatticeVector<MicrostructureGenerator::dim> L0 = rp.first;
                    const int grainID = rp.second;

                    std::uniform_int_distribution<> distribution(0, poly.grain(grainID).slipSystems().size() - 1);
                    const int rSS = distribution(generator); // a random SlipSystem
                    const auto &slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    const VectorDimD b(slipSystem.s.cartesian());

                    // Compute the ReciprocalLatticeDirection corresponding to s
                    ReciprocalLatticeDirection<3> sr(poly.grain(grainID).reciprocalLatticeDirection(b));

                    // find prismatic planes
                    std::vector<int> normalIDs;
                    for (size_t k = 0; k < poly.grain(grainID).planeNormals().size(); ++k)
                    {
                        if (slipSystem.s.dot(poly.grain(grainID).planeNormals()[k]) == 0)
                        {
                            normalIDs.push_back(k);
                        }
                    }
                    if (normalIDs.size() < 2)
                    {
                        std::cout << "normalIDs.size()=" << normalIDs.size() << std::endl;
                        std::cout << "Cannot generate a prismatic loop with less than 2 planes. EXITING." << std::endl;
                        exit(EXIT_FAILURE);
                    }

                    std::vector<LatticeDirection<MicrostructureGenerator::dim>> dirVector;
                    std::vector<int> sizeVector;
                    for (const int &normalID : normalIDs)
                    {
                        dirVector.emplace_back(sr, poly.grain(grainID).planeNormals()[normalID]);
                        sizeVector.emplace_back(randomSize() * randomSign() / dirVector.back().cartesian().norm());
                        //                        sizeVector.emplace_back(10000);
                    }

                    std::vector<VectorDimD> posVector;
                    std::vector<VectorDimD> normalsVector;

                    posVector.push_back(L0.cartesian());
                    for (size_t k = 0; k < dirVector.size(); ++k)
                    {
                        posVector.push_back(posVector[k] + dirVector[k].cartesian() * sizeVector[k]);
                        normalsVector.push_back(poly.grain(grainID).planeNormals()[normalIDs[k]].cartesian().normalized());
                    }
                    for (size_t k = 0; k < dirVector.size(); ++k)
                    {
                        posVector.push_back(posVector[k + dirVector.size()] - dirVector[k].cartesian() * sizeVector[k]);
                        normalsVector.push_back(-poly.grain(grainID).planeNormals()[normalIDs[k]].cartesian().normalized());
                    }
                    posVector.pop_back();
                    //                    std::cout<<"posVector.size()="<<posVector.size()<<std::endl;
                    //                    for(const auto& pos : posVector)
                    //                    {
                    //                        std::cout<<pos.cartesian().transpose()<<std::endl;
                    //                    }

                    bool allInside = true;
                    for (const auto &pos : posVector)
                    {
                        allInside *= mesh.searchRegion(grainID, pos).first;
                        if (!allInside)
                        {
                            break;
                        }
                    }

                    double dh = 0.0;
                    if (enforceMonotonicHelicity)
                    {
                        dh = deltaHelicity(posVector, -b);
                    }

                    if (allInside && ((fabs(helicity + dh) > fabs(helicity) && fabs(dh) > FLT_EPSILON) || helicity == 0.0 || !enforceMonotonicHelicity))
                    {

                        // Add nodes (two for-loops are needed)
                        for (size_t k = 0; k < posVector.size(); ++k)
                        {
                            configIO.nodes().emplace_back(nodeID + k, posVector[k], Eigen::Matrix<double, 1, 3>::Zero(), 1.0, 0);
                        }
                        for (size_t k = 0; k < posVector.size(); ++k)
                        {
                            configIO.nodes().emplace_back(nodeID + k + posVector.size(), posVector[k], Eigen::Matrix<double, 1, 3>::Zero(), 1.0, 0);
                        }

                        // Add lateral loops
                        size_t countLoopNode(0);
                        for (size_t k = 0; k < posVector.size(); ++k)
                        {
                            const int nextNodeID = (k + 1) < posVector.size() ? nodeID + k + 1 : nodeID;

                            configIO.loopNodes().emplace_back(loopNodeID + k +countLoopNode, loopID + k, posVector[nodeID + k- nodeID], nodeID + k, VectorDimD::Zero(), -1);
                            configIO.loopNodes().emplace_back(loopNodeID + k + 1+countLoopNode, loopID + k, posVector[nextNodeID- nodeID], nextNodeID, VectorDimD::Zero(), -1);
                            configIO.loopNodes().emplace_back(loopNodeID + k + 2+countLoopNode, loopID + k, posVector[nextNodeID- nodeID], nextNodeID + posVector.size(), VectorDimD::Zero(), -1);
                            configIO.loopNodes().emplace_back(loopNodeID + k + 3+countLoopNode, loopID + k, posVector[nodeID + k- nodeID], nodeID + k + posVector.size(), VectorDimD::Zero(), -1);
                            // configIO.loopNodes().emplace_back(loopNodeID + k + 4, loopID + k, posVector[nodeID + k], nodeID + k , VectorDimD::Zero(), -1);
                            configIO.loopLinks().emplace_back(loopID + k, loopNodeID+countLoopNode + k, loopNodeID+countLoopNode + k + 1 , true, 0);
                            configIO.loopLinks().emplace_back(loopID + k, loopNodeID+countLoopNode + k + 1, loopNodeID+countLoopNode + k + 2, true, 0);
                            configIO.loopLinks().emplace_back(loopID + k, loopNodeID+countLoopNode + k + 2, loopNodeID+countLoopNode + k + 3, true, 0);
                            configIO.loopLinks().emplace_back(loopID + k, loopNodeID+countLoopNode + k + 3, loopNodeID+countLoopNode + k , true, 0);
                            countLoopNode += 3;

                            configIO.loops().emplace_back(loopID + k, b, normalsVector[k], posVector[k], grainID, DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP);
                        }

                        // Add back loop (sessile)
                        for (size_t k = 0; k < posVector.size(); ++k)
                        {
                            const size_t nextLoopNodeID = (k + 1) < posVector.size() ? loopNodeID + k + 4*posVector.size() + 1 : loopNodeID + 4*posVector.size();
                            configIO.loopNodes().emplace_back(loopNodeID + 4*posVector.size()+k, loopID + posVector.size(), posVector[nodeID + k- nodeID], nodeID + k, VectorDimD::Zero(), -1);

                            configIO.loopLinks().emplace_back(loopID + posVector.size(), loopNodeID + 4*posVector.size() + k, nextLoopNodeID, true, 0);
                        }
                        configIO.loops().emplace_back(loopID + posVector.size(), -b, b.normalized(), posVector[0], grainID, DislocationLoopIO<MicrostructureGenerator::dim>::SESSILELOOP);

                        for (size_t k = 0; k < dirVector.size(); ++k)
                        {
                            density += 2.0 * (dirVector[k] * sizeVector[k]).cartesian().norm() / mesh.volume() / std::pow(poly.b_SI, 2);
                        }

                        std::cout << "density=" << density << std::endl;

                        nodeID += 2 * posVector.size();
                        //                        snID+=1;
                        loopID += posVector.size() + 1;
                        loopNodeID+=4 * posVector.size()+4;


                        if (enforceMonotonicHelicity)
                        {
                            loopPoints.push_back(posVector);
                            loopBurgers.push_back(-b);
                            helicity += dh;
                            std::cout << "helicity=" << helicity << std::endl;
                        }
                    }
                }
            }
        }

        /**********************************************************************/
        void MicrostructureGenerator::addIndividualStraightDislocations()
        {
            if(straightDislocationsSlipSystemIDs.size())
            {
                std::cout<<magentaBoldColor<<"Generating individual straight dislocations"<<defaultColor<<std::endl;
                if(straightDislocationsSlipSystemIDs.size()!=straightDislocationsAngleFromScrewOrientation.size())
                {
                    std::cout<<"straightDislocationsSlipSystemIDs.size()="<<straightDislocationsSlipSystemIDs.size()<<std::endl;
                    std::cout<<"straightDislocationsAngleFromScrewOrientation.size()="<<straightDislocationsAngleFromScrewOrientation.size()<<std::endl;
                    std::cout<<"You must provide one angle for each dislocation. EXITING."<<std::endl;
                    exit(EXIT_FAILURE);
                }
                if(int(straightDislocationsSlipSystemIDs.size())!=pointsAlongStraightDislocations.rows())
                {
                    std::cout<<"straightDislocationsSlipSystemIDs.size()="<<straightDislocationsSlipSystemIDs.size()<<std::endl;
                    std::cout<<"pointsAlongStraightDislocations.rows()="<<pointsAlongStraightDislocations.rows()<<std::endl;
                    std::cout<<"You must provide one point for each dislocation. Each point is a row of the matrix pointsAlongStraightDislocations. EXITING."<<std::endl;
                    exit(EXIT_FAILURE);
                }

                for(size_t k=0;k<straightDislocationsSlipSystemIDs.size();++k)
                {
                    const int& rSS(straightDislocationsSlipSystemIDs[k]);


                    if(rSS>=0)
                    {
                        std::pair<bool,const Simplex<MicrostructureGenerator::dim,MicrostructureGenerator::dim>*> found=mesh.search(pointsAlongStraightDislocations.row(k));
                        if(!found.first)
                        {
                            std::cout<<"Point "<<pointsAlongStraightDislocations.row(k)<<" is outside mesh. EXITING."<<std::endl;
                            exit(EXIT_FAILURE);
                        }

                        int grainID=found.second->region->regionID;



                        std::cout<<"generating individual straight dislocation "<<k<<defaultColor<<std::endl;


                        if(rSS>=int(poly.grain(grainID).slipSystems().size()))
                        {
                            std::cout<<"requested slip system ID="<<rSS<<std::endl;
                            std::cout<<"# of slip systems ="<<poly.grain(grainID).slipSystems().size()<<std::endl;
                            std::cout<<"Requested slip system does not exist. EXITING."<<std::endl;
                            exit(EXIT_FAILURE);
                        }

                        const auto& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                        const std::pair<bool,long int> heightPair=LatticePlane::computeHeight(slipSystem.n,pointsAlongStraightDislocations.row(k));

                        const VectorDimD P0=pointsAlongStraightDislocations.row(k).transpose()-pointsAlongStraightDislocations.row(k).dot(slipSystem.unitNormal)*slipSystem.unitNormal+slipSystem.unitNormal*slipSystem.n.planeSpacing()*heightPair.second;

                        const double theta(straightDislocationsAngleFromScrewOrientation[k]*M_PI/180.0);
                        //const VectorDimD& n(slipSystem.unitNormal);
                        const VectorDimD b(slipSystem.s.cartesian());

                        const VectorDimD d=Eigen::AngleAxisd(theta,slipSystem.unitNormal)*b.normalized();
                        const std::vector<VectorDimD> nodePos(DislocationInjector<MicrostructureGenerator::dim>::straightLineBoundaryClosure(P0,d,slipSystem.unitNormal,grainID,mesh));

                        const double lineLength=(nodePos.back()-nodePos.front()).norm();

                        


                        if(nodePos.size()>=3)
                        {// Write files

                            
                            
                            
//                            for(size_t k=0;k<nodePos.size();++k)
//                            {// write node and edge file
//                                configIO.nodes().emplace_back(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,0);
//                                const int nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
//                                configIO.loopLinks().emplace_back(loopID,nodeID+k,nextNodeID,0);
//                            }
//                            configIO.loops().emplace_back(loopID+0, b,n,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP);  // write loop file
//
//                            nodeID+=nodePos.size();
//                            loopID+=1;
//                            snID+=1;
                            std::cout<<"["<<b.transpose()<<"]("<<slipSystem.unitNormal.transpose()<<") dislocation. Line dir="<<d.transpose()<<". Length="<<lineLength<<std::endl;
                            if(lineLength<FLT_EPSILON)
                            {
                                std::cout<<"Line too short. EXITING."<<std::endl;
                                exit(EXIT_FAILURE);
                            }
                            
                            if(addSingleLoop(true,nodePos,nodePos,b,slipSystem.unitNormal,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP,std::vector<VectorDimD>(nodePos.size(),VectorDimD::Zero()),std::vector<short int>(nodePos.size(),-1)))
                            {
//                                addSingleLoop(true,nodePos,b,VectorDimD::Zero(),P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::SESSILELOOP)
                                std::cout<<"["<<b.transpose()<<"]("<<slipSystem.unitNormal.transpose()<<") dislocation. Line dir="<<d.transpose()<<". Length="<<lineLength<<std::endl;

                            }
                            
                        }
                        else
                        {
                            std::cout<<"nodePos.size="<<nodePos.size()<<std::endl;
                            assert(false && "LOOP DOES NOT HAVE ENOUGH POINTS");
                        }
                    }
                    else
                    {
                        std::cout<<"negative slip system ID. Skipping entries."<<std::endl;
                    }
                }
            }
        }

        /**********************************************************************/
        void MicrostructureGenerator::addNonPlanarLoops()
        {
            if(targetNonPlanarLoopDensity>0.0)
            {
                std::cout<<magentaBoldColor<<"Generating non-planar sessile loops"<<defaultColor<<std::endl;
                double density=0.0;
                while(density<targetNonPlanarLoopDensity)
                {
                    const std::pair<LatticeVector<MicrostructureGenerator::dim>,int> rp=randomPointInMesh();
                    const LatticeVector<MicrostructureGenerator::dim> L0=rp.first;
                    const VectorDimD P0(L0.cartesian());
                    const int grainID=rp.second;
                    
                    std::uniform_int_distribution<> distribution(0,poly.grain(grainID).slipSystems().size()-1);
                    const int rSS=distribution(generator); // a random SlipSystem
                    const auto& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    const VectorDimD b(poly.grain(grainID).latticeDirection(slipSystem.n.cartesian()).cartesian()); // Frank loop
                    const VectorDimD sessileAxis(b.normalized());
                    
                    std::normal_distribution<double> sizeDistribution(nonPlanarLoopRadiusMean/poly.b_SI,nonPlanarLoopRadiusStd/poly.b_SI);
                    const double radius(sizeDistribution(generator));
                    
                    std::normal_distribution<double> heightDistribution(nonPlanarLoopRadiusMean/poly.b_SI,0.5*nonPlanarLoopRadiusMean/poly.b_SI);

                    
                    std::vector<VectorDimD> nodePos;
                    for(int k=0;k<nonPlanarLoopSides;++k)
                    {
                        const double height(heightDistribution(generator));
                        nodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/nonPlanarLoopSides, sessileAxis)*slipSystem.s.cartesian().normalized()*radius+height*sessileAxis);
                    }
                    
                    if(addSingleLoop(true,nodePos,nodePos,b,VectorDimD::Zero(),P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::SESSILELOOP,std::vector<VectorDimD>(nodePos.size(),VectorDimD::Zero()),std::vector<short int>(nodePos.size(),-1)))
                    {
                        density += 2.0*radius*sin(M_PI/nonPlanarLoopSides)/mesh.volume()/std::pow(poly.b_SI,2);
                        std::cout<<"non-planar loop density="<<density<<std::endl;
                    }
                    
//                    if(allPointsInGrain(nodePos,grainID))
//                    {
//                        for(size_t k=0;k<nonPlanarLoopSides;++k)
//                        {
//                            const size_t nextNodeID=(k+1)<nodePos.size()? nodeID+k+1 : nodeID;
//                            configIO.nodes().emplace_back(nodeID+k,nodePos[k],Eigen::Matrix<double,1,3>::Zero(),1.0,0);
//                            configIO.loopLinks().emplace_back(loopID,nodeID+k,nextNodeID,0);
//                        }
//                        configIO.loops().emplace_back(loopID+0, b,VectorDimD::Zero(),P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::SESSILELOOP);  // write loop file
//                        nodeID+=nodePos.size();
//                        loopID+=1;
//                        snID+=1;
//                        density += 2.0*radius*sin(M_PI/nonPlanarLoopSides)/mesh.volume()/std::pow(poly.b_SI,2);
//                        std::cout<<"non-planar loop density="<<density<<std::endl;
//                    }
                    
                }
                
            }
            
        }
        
        
        /**********************************************************************/
        void MicrostructureGenerator::addStatisticallyHomegeneousPeriodicLoops()
        {

            const size_t statisticalHomogeneous = TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<size_t>("statisticalHomogeneous", true);
            if (targetPeriodicLoopDensity > 0.0)
            {
                std::cout << magentaBoldColor << "Generating periodic loops" << defaultColor << std::endl;
                assert(poly.grains().size() == 1 && "PeriodicLoops only supported in single crystals.");
                double density = 0.0;
                size_t periodicNodeID(0);
                //                size_t periodicLoopID(0);
                while (density < targetPeriodicLoopDensity)
                {
                    std::vector<std::pair<LatticeVector<MicrostructureGenerator::dim>, int>> rps;
                    std::vector<LatticeVector<MicrostructureGenerator::dim>> L0s;
                    rps.emplace_back(randomPointInMesh());
                    L0s.emplace_back(rps.begin()->first);
                    rps.emplace_back(randomPointInMesh());
                    L0s.emplace_back(rps.rbegin()->first);

                    const int grainID = rps.begin()->second;
                    std::uniform_int_distribution<> distribution(0, poly.grain(grainID).slipSystems().size() - 1);
                    const int rSS = distribution(generator); // a random SlipSystem
                    const auto &slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    const VectorDimD b(slipSystem.s.cartesian());
                    std::vector<VectorDimD> bs;
                    bs.push_back(b);
                    bs.push_back(-b);
                    size_t numLoopinOne = (statisticalHomogeneous ? 2 : 1);
                    for (size_t i = 0; i < numLoopinOne; i++)
                    {
                        const VectorDimD P0(L0s[i].cartesian());
                        std::normal_distribution<double> sizeDistribution(periodicLoopRadiusMean / poly.b_SI, periodicLoopRadiusStd / poly.b_SI);
                        const double radius(sizeDistribution(generator));

                        std::vector<PolyPoint> dummyPolyPoints;
                        std::vector<std::pair<VectorDimD, const PolyPoint *const>> loopNodePosTemp;
                        for (int k = 0; k < periodicLoopSides; ++k)
                        {
                            dummyPolyPoints.push_back(PolyPoint());
                            loopNodePosTemp.push_back(std::make_pair(P0 + Eigen::AngleAxisd(k * 2.0 * M_PI / periodicLoopSides, slipSystem.unitNormal) * slipSystem.s.cartesian().normalized() * radius, &dummyPolyPoints.back()));
                        }

                        PeriodicGlidePlaneFactory<MicrostructureGenerator::dim> pgpf(poly, glidePlaneFactory);
                        GlidePlaneKey<3> glidePlaneKey(P0, slipSystem.n);
                        std::shared_ptr<PeriodicGlidePlane<3>> periodicGlidePlane(pgpf.get(glidePlaneKey));

                        const auto ppi(periodicGlidePlane->polygonPatchIntersection(loopNodePosTemp));

                        std::vector<VectorDimD> loopNodePos;
                        std::vector<VectorDimD> networkNodePos;
                        std::vector<VectorDimD> loopNodeShifts;
                        std::vector<short int> edgeIDs;

                        for (const auto &tup : ppi)
                        {
                            const VectorDimD gblP(periodicGlidePlane->referencePlane->globalPosition(std::get<0>(tup)));
                            loopNodePos.push_back(gblP);
                            networkNodePos.push_back(gblP + std::get<1>(tup));
                            loopNodeShifts.push_back(std::get<1>(tup));
                            edgeIDs.push_back(std::get<2>(tup));
                        }

                        //
                        if (addSingleLoop(false, networkNodePos, loopNodePos, bs[i], slipSystem.unitNormal, P0, grainID, DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP, loopNodeShifts, edgeIDs))
                        {
                            density += 2.0 *periodicLoopSides * radius * sin(M_PI / periodicLoopSides) / mesh.volume() / std::pow(poly.b_SI, 2);

                            std::cout << "periodicLoop density=" << density << std::endl;
                        }
                    }
                }
            }
        }


        /**********************************************************************/
        void MicrostructureGenerator::addFrankLoops()
        {
            if(targetFrankLoopDensity>0.0)
            {
                std::cout<<magentaBoldColor<<"Generating Frank loops"<<defaultColor<<std::endl;
                double density=0.0;
                while(density<targetFrankLoopDensity)
                {
                    const std::pair<LatticeVector<MicrostructureGenerator::dim>,int> rp=randomPointInMesh();
                    const LatticeVector<MicrostructureGenerator::dim> L0=rp.first;
                    const VectorDimD P0(L0.cartesian());
                    const int grainID=rp.second;

                    std::uniform_int_distribution<> planesDistribution(0,poly.grain(grainID).planeNormals().size()-1);
                    const int rSS=planesDistribution(generator); // a random SlipSystem
                    const auto& n(poly.grain(grainID).planeNormals()[rSS]);
                    const VectorDimD unitNormal(n.cartesian().normalized());
                    const VectorDimD b(n.planeSpacing()*unitNormal); // Frank loop is one plane of vanancies/interstitials

                    std::normal_distribution<double> sizeDistribution(frankLoopRadiusMean/poly.b_SI,frankLoopRadiusStd/poly.b_SI);
                    const double radius(sizeDistribution(generator));
                    const VectorDimD R(randomOrthogonalUnitVector(b)*radius);

                    std::vector<VectorDimD> nodePos;
                    for(int k=0;k<frankLoopSides;++k)
                    {
                        nodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/frankLoopSides, unitNormal)*R);
                    }
                    
                    if(addSingleLoop(true,nodePos,nodePos,b,unitNormal,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::SESSILELOOP,std::vector<VectorDimD>(nodePos.size(),VectorDimD::Zero()),std::vector<short int>(nodePos.size(),-1)))
                    {
                        density += 2.0*radius*sin(M_PI/frankLoopSides)/mesh.volume()/std::pow(poly.b_SI,2);
                        std::cout<<"Frank loop density="<<density<<std::endl;
                    }
                }
            }

        }

        /**********************************************************************/
        void MicrostructureGenerator::addEshelbyInclusions()
        {

            if(targetInclusionDensities.size())
            {
                double totalDensity(0.0);
                for(const double& val : targetInclusionDensities)
                {
                    if(val>0.0)
                    {
                        totalDensity+=val;
                    }
                }

                if(totalDensity)
                {
                    assert(targetInclusionDensities.size()==inclusionsDiameterLognormalDistribution_M.size());
                    assert(targetInclusionDensities.size()==inclusionsDiameterLognormalDistribution_S.size());
                    assert(targetInclusionDensities.size()==inclusionsDiameterLognormalDistribution_A.size());
                    assert(int(targetInclusionDensities.size())==inclusionsTransformationStrains.rows());
                    assert(int(targetInclusionDensities.size())==inclusionsPatterns.rows());


                    std::cout<<magentaBoldColor<<"Generating Inclusions"<<defaultColor<<std::endl;

                    std::ofstream inclusionsfile("E/E_0.txt");
                    std::deque<std::pair<double,VectorDimD>> existingPrecipitates;

                    size_t inclusionID=0;
                    for(size_t f=0;f<targetInclusionDensities.size();++f)
                    {
                        if(   targetInclusionDensities[f]>0.0)
                        {
                            assert(inclusionsDiameterLognormalDistribution_M[f]>0.0);
                            assert(inclusionsDiameterLognormalDistribution_S[f]>0.0);
                            assert(inclusionsDiameterLognormalDistribution_A[f]>0.0);

                            const VectorDimD currentPattern(inclusionsPatterns.row(f)/poly.b_SI); // normalize to length units
                            const double patternHeigth(currentPattern.norm());
                            const bool applyPattern(patternHeigth>0.0);
                            const VectorDimD patternDir(applyPattern? (currentPattern/patternHeigth).eval() : VectorDimD::Zero());

                            double numberDensity=0.0;

                            std::lognormal_distribution<double> distribution(log(inclusionsDiameterLognormalDistribution_M[f]/inclusionsDiameterLognormalDistribution_A[f]),inclusionsDiameterLognormalDistribution_S[f]);

                            while(numberDensity<targetInclusionDensities[f])
                            {
                                const double diameter = distribution(generator)*inclusionsDiameterLognormalDistribution_A[f]/poly.b_SI;
                                const double radius(0.5*diameter);
                                std::pair<LatticeVector<MicrostructureGenerator::dim>,int> pointPair=randomPointInMesh();
                                VectorDimD P=pointPair.first.cartesian();
                                const int& grainID(pointPair.second);

                                if(applyPattern)
                                {
                                    const VectorDimD globalVector(poly.grain(grainID).C2G*currentPattern);
                                    const VectorDimD globalDir(poly.grain(grainID).C2G*patternDir);
                                    const long long pointHeigth=std::round(P.dot(globalDir)/patternHeigth);
                                    const VectorDimD O(pointHeigth*globalVector);
                                    P-=(P-O).dot(globalDir)*globalDir;
                                }

                                bool isGoodPosition=mesh.searchRegion(grainID,P).first;
                                for(const auto& pair : existingPrecipitates)
                                {
                                    isGoodPosition *= (P-pair.second).norm()>pair.first+radius;
                                    if(!isGoodPosition)
                                    {
                                        break;
                                    }
                                }

                                if(isGoodPosition)
                                {
                                    inclusionsfile<<inclusionID
                                    /*          */<<" "<<P.transpose()
                                    /*          */<<" "<<radius
                                    /*          */<<" "<<inclusionsTransformationStrains.row(f)
                                    /*          */<<" "<<f
                                    /*          */<<"\n";

                                    numberDensity+=1.0/mesh.volume()/std::pow(poly.b_SI,3);
                                    inclusionID++;
                                    existingPrecipitates.emplace_back(radius,P);

                                    std::cout<<"inclusions density="<<numberDensity<<std::endl;
                                }
                            }
                        }
                    }
                    inclusionsfile.close();
                }
            }
        }

        

        bool MicrostructureGenerator::isInclusionsUsed() const
        {
            bool temp(false);
            for(const auto& density : targetInclusionDensities)
            {
                if(density>0.0)
                {
                    temp=true;
                    break;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
       MicrostructureGenerator::MicrostructureGenerator(int argc, char* argv[]) :
        /* init*/ generator(std::chrono::system_clock::now().time_since_epoch().count())
        /* init*/,nodeID(0)
       /* init*/,loopNodeID(0)
        /* init*/,loopID(0)
        /* init*/,enforceMonotonicHelicity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("enforceMonotonicHelicity",true))
        /* init*/,helicity(0.0)
        /* init*/,outputBinary(TextFileParser("./inputFiles/DD.txt").readScalar<int>("outputBinary",true))
//        /* init*/,meshID(TextFileParser("./inputFiles/DD.txt").readScalar<int>("meshID",true))
        /* init*/,meshFilename(TextFileParser("./inputFiles/polycrystal.txt").readString("meshFile",true))
//        /* init*/,mesh(meshFilename)
        /* init */,mesh(meshFilename,TextFileParser("./inputFiles/polycrystal.txt").readMatrix<double>("A",3,3,true),TextFileParser("./inputFiles/polycrystal.txt").readMatrix<double>("x0",1,3,true).transpose())
        /* init*/,minSize(0.1*min(mesh.xMax(0)-mesh.xMin(0),min(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2))))
        /* init*/,maxSize(max(mesh.xMax(0)-mesh.xMin(0),max(mesh.xMax(1)-mesh.xMin(1),mesh.xMax(2)-mesh.xMin(2))))
        /* init*/,poly("./inputFiles/polycrystal.txt",mesh)
        /* init*/,glidePlaneFactory(poly)
        /* Straight Dislocations */
        /* init*/,targetStraightDislocationDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetStraightDislocationDensity",true))
        /* init*/,fractionSessileStraightDislocationDensity(targetStraightDislocationDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("fractionSessileStraightDislocationDensity",true) : 0.0)
        /* Frank-Read sources */
        /* init*/,targetFrankReadDislocationDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetFrankReadDislocationDensity",true))
        /* init*/,FrankReadSizeMean(targetFrankReadDislocationDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("FrankReadSizeMean",true) : 0.0)
        /* init*/,FrankReadSizeStd(targetFrankReadDislocationDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("FrankReadSizeStd",true) : 0.0)
        /* init*/,FrankReadAspectRatioMean(targetFrankReadDislocationDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("FrankReadAspectRatioMean",true) : 0.0)
        /* init*/,FrankReadAspectRatioStd(targetFrankReadDislocationDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("FrankReadAspectRatioStd",true) : 0.0)
        /* Single-arm sources */
        /* init*/,targetSingleArmDislocationDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetSingleArmDislocationDensity",true))
        /* Prismatic loops */
        /* init*/,targetPrismaticLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetPrismaticLoopDensity",true))
        /* Indivial straight dislocations */
        /* init*/,straightDislocationsSlipSystemIDs(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<int>("straightDislocationsSlipSystemIDs",true))
        /* init*/,straightDislocationsAngleFromScrewOrientation(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("straightDislocationsAngleFromScrewOrientation",true))
        /* init*/,pointsAlongStraightDislocations(TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("pointsAlongStraightDislocations",straightDislocationsSlipSystemIDs.size(),MicrostructureGenerator::dim,true))
        /* Circular Loops */
        /* init*/,targetFrankLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetFrankLoopDensity",true))
        /* init*/,frankLoopRadiusMean(targetFrankLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("frankLoopRadiusMean",true) : 0.0)
        /* init*/,frankLoopRadiusStd(targetFrankLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("frankLoopRadiusStd",true) : 0.0)
        /* init*/,frankLoopSides(targetFrankLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("frankLoopSides",true) : 0.0)
        /* NonPlanar Loops */
        /* init*/,targetNonPlanarLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetNonPlanarLoopDensity",true))
        /* init*/,nonPlanarLoopRadiusMean(targetNonPlanarLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("nonPlanarLoopRadiusMean",true) : 0.0)
        /* init*/,nonPlanarLoopRadiusStd(targetNonPlanarLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("nonPlanarLoopRadiusStd",true) : 0.0)
        /* init*/,nonPlanarLoopSides(targetNonPlanarLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("nonPlanarLoopSides",true) : 0.0)
        /* PeriodicLoops */
        /* init*/,targetPeriodicLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetPeriodicLoopDensity",true))
        /* init*/,periodicLoopRadiusMean(targetPeriodicLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("periodicLoopRadiusMean",true) : 0.0)
        /* init*/,periodicLoopRadiusStd(targetPeriodicLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("periodicLoopRadiusStd",true) : 0.0)
        /* init*/,periodicLoopSides(targetPeriodicLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("periodicLoopSides",true) : 0.0)
//        /* JunctionLoops */
//        /* init*/,targetJunctionLoops(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<size_t>("targetJunctionLoops",true))
//        /* init*/,targetJunctionLoopsSize(targetJunctionLoops>0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetJunctionLoopsSize",true) : 0.0)
//        /* Planar Dipolar Loops */
//        /* init*/,targetPlanarDipolarLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetPlanarDipolarLoopDensity",true))
//        /* init*/,planarDipolarLoopMean(targetPlanarDipolarLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("planarDipolarLoopMean",true) : 0.0)
//        /* init*/,planarDipolarLoopStd(targetPlanarDipolarLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("planarDipolarLoopStd",true) : 0.0)
//        /* init*/,planarDipolarLoopAspectRatioMean(targetPlanarDipolarLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("planarDipolarLoopAspectRatioMean",true) : 0.0)
//        /* init*/,planarDipolarLoopAspectRatioStd(targetPlanarDipolarLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("planarDipolarLoopAspectRatioStd",true) : 0.0)

        /* Irradiation Loops */
        /* init*/,targetIrradiationLoopDensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetIrradiationLoopDensity",true))
        /* init*/,irradiationLoopsDiameterLognormalDistribution_M(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("irradiationLoopsDiameterLognormalDistribution_M",true) : 0.0)
        /* init*/,irradiationLoopsDiameterLognormalDistribution_S(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("irradiationLoopsDiameterLognormalDistribution_S",true) : 0.0)
        /* init*/,irradiationLoopsDiameterLognormalDistribution_A(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("irradiationLoopsDiameterLognormalDistribution_A",true) : 0.0)
        /* init*/,fraction111Loops(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("fraction111Loops",true) : 0.0)
        /* init*/,mobile111Loops(targetIrradiationLoopDensity>0.0? TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("mobile111Loops",true) : 0.0)
        /* SFTs */
        /* init*/,targetSFTdensity(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<double>("targetSFTdensity",true))
        /* Inclusions */
        /* init*/,targetInclusionDensities(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("targetInclusionDensities",true))
        /* init*/,inclusionsDiameterLognormalDistribution_M(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDiameterLognormalDistribution_M",true) : std::vector<double>())
        /* init*/,inclusionsDiameterLognormalDistribution_S(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDiameterLognormalDistribution_S",true) : std::vector<double>())
        /* init*/,inclusionsDiameterLognormalDistribution_A(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsDiameterLognormalDistribution_A",true) : std::vector<double>())
        /* init*/,inclusionsTransformationStrains(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("inclusionsTransformationStrains",targetInclusionDensities.size(),MicrostructureGenerator::dim*MicrostructureGenerator::dim,true) : Eigen::Matrix<double,Eigen::Dynamic,MicrostructureGenerator::dim*MicrostructureGenerator::dim>::Zero(1,MicrostructureGenerator::dim*MicrostructureGenerator::dim))
        /* init*/,inclusionsPatterns(isInclusionsUsed()? TextFileParser("./inputFiles/initialMicrostructure.txt").readMatrix<double>("inclusionsPatterns",targetInclusionDensities.size(),MicrostructureGenerator::dim,true) : Eigen::Matrix<double,Eigen::Dynamic,MicrostructureGenerator::dim>::Zero(1,MicrostructureGenerator::dim))
        //        /* init*/,inclusionsMobilityReduction(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsMobilityReduction",true))
        {

            

            
            // Some sanity checks
            if(mesh.volume()==0.0)
            {
                std::cout<<"mesh "<<meshFilename<<" is empty. MicrostructureGenerator cannot run. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }

            std::cout<<greenBoldColor<<"Generating initial microstructure"<<defaultColor<<std::endl;
            // Call individual generators
            addStraightDislocations();
            addFrankReadSources();
            addSingleArmDislocations();
            addPrismaticLoops();
            addIndividualStraightDislocations();
            addFrankLoops();
            addNonPlanarLoops();
            // addPeriodicLoops();
            addStatisticallyHomegeneousPeriodicLoops();
//            addStatisticallyHomegeneousPlanarDipolarLoops();
//            addPeriodicJunctionLoops();
            addIrradiationLoops();
            addStackingFaultTetrahedra();
            addEshelbyInclusions();
            writeConfigFiles(0);
            
        }
        
        /**********************************************************************/
        void MicrostructureGenerator::writeConfigFiles(const size_t& fileID)
        {

            auxIO.setGlidePlaneBoundaries(glidePlaneFactory); // change this function to take a GlidePlaneFactory during write
            
            if(outputBinary)
            {
                std::cout<<greenBoldColor<<"Writing configuration to "<<configIO.getBinFilename(fileID)<<defaultColor<<std::endl;
                configIO.writeBin(fileID);
                auxIO.writeBin(fileID);
            }
            else
            {
                std::cout<<greenBoldColor<<"Writing configuration to "<<configIO.getTxtFilename(fileID)<<defaultColor<<std::endl;
                configIO.writeTxt(fileID);
                auxIO.writeTxt(fileID);
            }
        }

        /**********************************************************************/
        void MicrostructureGenerator::addStackingFaultTetrahedra()
        {

            if(targetSFTdensity>0.0)
            {
                std::cout<<magentaBoldColor<<"Generating Stacking Fault Tetrahedra"<<defaultColor<<std::endl;


                if(poly.crystalStructure=="FCC")
                {
                    size_t ndefects=0;
                    double defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);

                    const std::pair<LatticeVector<MicrostructureGenerator::dim>,int> rp=randomPointInMesh();
                    const int& grainID=rp.second;   // random grain ID
                    const LatticeVector<MicrostructureGenerator::dim>& L0=rp.first; // random lattice position in the grain


                    LatticeVector<3> a1(VectorDimI(1,0,0),poly.grain(grainID)); // [011]
                    LatticeVector<3> a2(VectorDimI(0,1,0),poly.grain(grainID)); // [101]
                    LatticeVector<3> a3(VectorDimI(0,0,1),poly.grain(grainID)); // [110]
                    LatticeVector<3> a12(a2-a1); // [1-10]
                    LatticeVector<3> a23(a3-a2); // [01-1]
                    LatticeVector<3> a31(a1-a3); // [-101]

                    int Li=100;
                    const VectorDimD P0(L0.cartesian());
                    const VectorDimD P1(P0+a1.cartesian()*Li);
                    const VectorDimD P2(P0+a2.cartesian()*Li);
                    const VectorDimD P3(P0+a3.cartesian()*Li);

                    if(   mesh.searchRegion(grainID,P1).first
                       && mesh.searchRegion(grainID,P1).first
                       && mesh.searchRegion(grainID,P2).first)
                    {
                        
                        assert(false && "FINISH IMPLEMENTATION< CHECK ALL BURGERS VECTORS");
                        
                        VectorDimD b=L0.cartesian();

                        configIO.nodes().emplace_back(nodeID+0,P0,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                        configIO.nodes().emplace_back(nodeID+1,P1,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                        configIO.nodes().emplace_back(nodeID+2,P2,Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                        configIO.nodes().emplace_back(nodeID+3,P3,Eigen::Matrix<double,1,3>::Zero(),1.0,0);

                        configIO.loopLinks().emplace_back(loopID+0,nodeID+0,nodeID+2,true,0);
                        configIO.loopLinks().emplace_back(loopID+0,nodeID+2,nodeID+1,true,0);
                        configIO.loopLinks().emplace_back(loopID+0,nodeID+1,nodeID+0,true,0);
                        configIO.loops().emplace_back(loopID+0, b,a2.cross(a1).cartesian(),P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP);

                        configIO.loopLinks().emplace_back(loopID+1,nodeID+0,nodeID+1,true,0);
                        configIO.loopLinks().emplace_back(loopID+1,nodeID+1,nodeID+3,true,0);
                        configIO.loopLinks().emplace_back(loopID+1,nodeID+3,nodeID+0,true,0);
                        configIO.loops().emplace_back(loopID+1, b,a1.cross(a3).cartesian(),P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP);

                        configIO.loopLinks().emplace_back(loopID+2,nodeID+0,nodeID+3,true,0);
                        configIO.loopLinks().emplace_back(loopID+2,nodeID+3,nodeID+2,true,0);
                        configIO.loopLinks().emplace_back(loopID+2,nodeID+2,nodeID+0,true,0);
                        configIO.loops().emplace_back(loopID+2, b,a3.cross(a2).cartesian(),P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP);

                        configIO.loopLinks().emplace_back(loopID+3,nodeID+1,nodeID+2,true,0);
                        configIO.loopLinks().emplace_back(loopID+3,nodeID+2,nodeID+3,true,0);
                        configIO.loopLinks().emplace_back(loopID+3,nodeID+3,nodeID+1,true,0);
                        configIO.loops().emplace_back(loopID+3, b,a31.cross(a12).cartesian(),P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP);


//                        snID++;
                        loopID+=4;
                        nodeID+=4;
                        ndefects++;
                        defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                        std::cout<<"SFT density="<<defectsDensity<<std::endl;


                    }

                }
                else
                {
                    std::cout<<"SFTs can only be generated for FCC crystals"<<std::endl;
                }
            }

        }

        /**********************************************************************/
        bool MicrostructureGenerator::allPointsInGrain(const std::vector<VectorDimD>& points,const int& grainID)
        {
            bool temp=true;
            for(const auto& point : points)
            {
                temp*=mesh.searchRegion(grainID,point).first;
                if(!temp)
                {
                    break;
                }
            }
            return temp;
        }

        /**********************************************************************/
        void MicrostructureGenerator::addIrradiationLoopsFCC()
        {// Irradiation loops in FCC are Frank loops

            const size_t irradiationLoopsNumberOfNodes(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("irradiationLoopsNumberOfNodes",true));
            std::lognormal_distribution<double> sizeDistribution(log(irradiationLoopsDiameterLognormalDistribution_M/irradiationLoopsDiameterLognormalDistribution_A),(irradiationLoopsDiameterLognormalDistribution_S));

            size_t ndefects=0;
            double defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
            while(defectsDensity<targetIrradiationLoopDensity)
            {
                const std::pair<LatticeVector<MicrostructureGenerator::dim>,int> rp=randomPointInMesh();
                const int& grainID=rp.second;   // random grain ID
                const LatticeVector<MicrostructureGenerator::dim>& L0=rp.first; // random lattice position in the grain
                const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                
                std::uniform_int_distribution<> planesDistribution(0,poly.grain(grainID).planeNormals().size()-1);
                const int rSS=planesDistribution(generator);
                const auto& n(poly.grain(grainID).planeNormals()[rSS]); // a random {111} plane
                const VectorDimD unitNormal(n.cartesian().normalized());
                const VectorDimD b(n.planeSpacing()*unitNormal); // Frank loops are plateletslane of vanancies/interstitials

                const double diameter_SI = sizeDistribution(generator)*irradiationLoopsDiameterLognormalDistribution_A;
                const double radius(0.5*diameter_SI/poly.b_SI);
                const VectorDimD R(randomOrthogonalUnitVector(b)*radius);

                std::vector<VectorDimD> nodePos;
                for(size_t k=0;k<irradiationLoopsNumberOfNodes;++k)
                {
                    nodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/irradiationLoopsNumberOfNodes, unitNormal)*R);
                }
                
                if(addSingleLoop(true,nodePos,nodePos,b,unitNormal,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::SESSILELOOP,std::vector<VectorDimD>(nodePos.size(),VectorDimD::Zero()),std::vector<short int>(nodePos.size(),-1)))
                {
                    ndefects++;
                    defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                    std::cout<<"irradiation loops density="<<defectsDensity<<std::endl;
                }
            }
        }
        
        /**********************************************************************/
        void MicrostructureGenerator::addIrradiationLoopsHCP()
        {// Irradiation loops in FCC are Frank loops.
         // TO DO: actaully Frank loops can be unfaulted and become prismatic (exagonal) loops bounded by pairs of {111} planes and {001} planes. See P. B. Hirsch , J. Silcox , R. E. Smallman & K. H. Westmacott
            
            const size_t irradiationLoopsNumberOfNodes(TextFileParser("./inputFiles/initialMicrostructure.txt").readScalar<int>("irradiationLoopsNumberOfNodes",true));
            std::lognormal_distribution<double> sizeDistribution(log(irradiationLoopsDiameterLognormalDistribution_M/irradiationLoopsDiameterLognormalDistribution_A),(irradiationLoopsDiameterLognormalDistribution_S));
            
            size_t ndefects=0;
            double defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
            while(defectsDensity<targetIrradiationLoopDensity)
            {
                const std::pair<LatticeVector<MicrostructureGenerator::dim>,int> rp=randomPointInMesh();
                const int& grainID=rp.second;   // random grain ID
                const LatticeVector<MicrostructureGenerator::dim>& L0=rp.first; // random lattice position in the grain
                const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                
                std::uniform_int_distribution<> slipSystemDistribution(0,poly.grain(grainID).slipSystems().size()-1);
                const int rSS=slipSystemDistribution(generator); // a random SlipSystem ID
                const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);

                const double diameter_SI = sizeDistribution(generator)*irradiationLoopsDiameterLognormalDistribution_A;
                const double radius(0.5*diameter_SI/poly.b_SI);
                const VectorDimD R(slipSystem.s.cartesian().normalized()*radius);

                
                std::vector<VectorDimD> nodePos;
                for(size_t k=0;k<irradiationLoopsNumberOfNodes;++k)
                {
                    nodePos.push_back(P0+Eigen::AngleAxisd(k*2.0*M_PI/irradiationLoopsNumberOfNodes, slipSystem.unitNormal)*R);
                }
                
                const double planeSpacing(slipSystem.n.planeSpacing());
                if(fabs(planeSpacing-sqrt(3.0)/2.0)<FLT_EPSILON)
                {// prismatic plane spacing
                    const VectorDimD e=slipSystem.s.cartesian().cross(slipSystem.n.cartesian()).normalized();                 // "edge" direction, along prism axis
                    const VectorDimD b=Eigen::AngleAxisd(randomSign()*M_PI/3.0, e)*(slipSystem.s.cartesian());                  // rotate slip direction out of the plane by 60 deg
                    if(addSingleLoop(true,nodePos,nodePos,b,slipSystem.unitNormal,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::SESSILELOOP,std::vector<VectorDimD>(nodePos.size(),VectorDimD::Zero()),std::vector<short int>(nodePos.size(),-1)))
                    {
                        ndefects++;
                        defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                        std::cout<<"irradiation loops density="<<defectsDensity<<std::endl;
                    }
                }
                else if(fabs(planeSpacing-sqrt(8.0/3.0))<FLT_EPSILON)
                {// basal plane spacing
                    const VectorDimD b(0.5*slipSystem.n.planeSpacing()*slipSystem.n.cartesian().normalized()); // 1/2 c-type loop
                    if(addSingleLoop(true,nodePos,nodePos,b,slipSystem.unitNormal,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::SESSILELOOP,std::vector<VectorDimD>(nodePos.size(),VectorDimD::Zero()),std::vector<short int>(nodePos.size(),-1)))
                    {
                        ndefects++;
                        defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                        std::cout<<"irradiation loops density="<<defectsDensity<<std::endl;
                    }
                }
                else
                {
                    
                }
            }
        }
        

        /**********************************************************************/
        void MicrostructureGenerator::addIrradiationLoopsBCC()
        {
            size_t ndefects=0;
            double defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
            //                    int NP=6;
            std::lognormal_distribution<double> distribution(log(irradiationLoopsDiameterLognormalDistribution_M/irradiationLoopsDiameterLognormalDistribution_A),(irradiationLoopsDiameterLognormalDistribution_S));
            
            while(defectsDensity<targetIrradiationLoopDensity)
            {
                const std::pair<LatticeVector<MicrostructureGenerator::dim>,int> rp=randomPointInMesh();
                const int& grainID=rp.second;   // random grain ID
                const LatticeVector<MicrostructureGenerator::dim>& L0=rp.first; // random lattice position in the grain
                const VectorDimD P0(L0.cartesian());   // cartesian position of L0
                
                if (defectsDensity<targetIrradiationLoopDensity*fraction111Loops)
                {// add [111] loops
                    std::uniform_int_distribution<> slipSystemDistribution(0,poly.grain(grainID).slipSystems().size()-1);
                    const int rSS=slipSystemDistribution(generator); // a random SlipSystem ID
                    const SlipSystem& slipSystem(*poly.grain(grainID).slipSystems()[rSS]);
                    const VectorDimD b=slipSystem.s.cartesian();    // Burgers vector
                    const VectorDimD a=b.normalized();
                    
                    const double diameter_SI = distribution(generator)*irradiationLoopsDiameterLognormalDistribution_A;
                    const double radius_SI(0.5*diameter_SI);
                    const double radius(std::round(radius_SI/poly.b_SI/(2.0*sqrt(2.0)/3.0))*(2.0*sqrt(2.0)/3.0)); // radius must be an integer multiple of 2.0*sqrt(2.0)/3.0, or lateral sides will not be on glide planes
                    
                    if(radius>0.0)
                    {
                        std::vector<VectorDimD> points;
                        for(int k=0;k<6;++k)
                        {
                            points.push_back(P0+Eigen::AngleAxis<double>(k*2.0*M_PI/6-M_PI/6,a)*slipSystem.unitNormal*radius);
                        }
                        
                        if(allPointsInGrain(points,grainID))
                        {
                            
                            for(int k=0;k<6;++k)
                            {// inser the back loop
                                configIO.nodes().emplace_back(nodeID+k,points[k],Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                                const int nextNodeID=(k+1)<6? nodeID+k+1 : nodeID;
                                configIO.loopLinks().emplace_back(loopID,nodeID+k,nextNodeID,true,0);
                            }
                            configIO.loops().emplace_back(loopID+0, b,a,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::SESSILELOOP);
                            loopID++;
                            
                            
                            if(mobile111Loops)
                            {// insert loops on the six sides
                                for(int k=0;k<6;++k)
                                {// inser lateral loops
                                    configIO.nodes().emplace_back(nodeID+k+6,points[k],Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                                    const int nextNodeID=(k+1)<6? nodeID+k+1 : nodeID;
                                    configIO.loopLinks().emplace_back(loopID,nodeID+k,nodeID+k+6,true,0);
                                    configIO.loopLinks().emplace_back(loopID,nodeID+k+6,nextNodeID+6,true,0);
                                    configIO.loopLinks().emplace_back(loopID,nextNodeID+6,nextNodeID,true,0);
                                    configIO.loopLinks().emplace_back(loopID,nextNodeID,nodeID+k,true,0);
                                    
                                    configIO.loops().emplace_back(loopID+0, b,Eigen::AngleAxis<double>(k*2.0*M_PI/6,a)*slipSystem.unitNormal,points[k],grainID,DislocationLoopIO<MicrostructureGenerator::dim>::GLISSILELOOP);
                                    loopID++;
                                }
                                nodeID+=2*points.size();
                            }
                            else
                            {
                                nodeID+=points.size();
                            }
//                            snID++;
                            ndefects++;
                            defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                            std::cout<<"irradiation defects density="<<defectsDensity<<std::endl;
                        }
                    }
                }
                else
                {// add [100] loop (by Yinan Cui)
                    std::vector<LatticeDirectionType> sessileb;
                    sessileb.emplace_back(LatticeVector<MicrostructureGenerator::dim>(VectorDimI(0,1,1),poly.grain(grainID).lattice())); // is ( 1, 0, 0) in cartesian
                    sessileb.emplace_back(LatticeVector<MicrostructureGenerator::dim>(VectorDimI(1,0,1),poly.grain(grainID).lattice())); // is ( 0, 1, 0) in cartesian
                    sessileb.emplace_back(LatticeVector<MicrostructureGenerator::dim>(VectorDimI(1,1,0),poly.grain(grainID).lattice())); // is ( 0, 0, 1) in cartesian
                    
                    std::uniform_int_distribution<> dist(0,2);
                    typedef Eigen::Matrix<long int,2,1> Vector2I;
                    std::vector<Vector2I> sslinedirection;
                    sslinedirection.emplace_back(Vector2I(1,2));
                    sslinedirection.emplace_back(Vector2I(0,2));
                    sslinedirection.emplace_back(Vector2I(0,1));
                    
                    //NP=4;
                    const int rSS_sessile=dist(generator); // a random sessile plane
                    LatticeDirection<3> d1(sessileb[sslinedirection[rSS_sessile][0]]);
                    LatticeDirection<3> d2(sessileb[sslinedirection[rSS_sessile][1]]);
                    
                    const double diameter_SI = distribution(generator)*irradiationLoopsDiameterLognormalDistribution_A;
                    double a1=diameter_SI/poly.b_SI;
                    double a2=diameter_SI/poly.b_SI;
                    
                    std::vector<VectorDimD> points;
                    points.push_back(P0);
                    points.push_back(P0+d1.cartesian().normalized()*a1);
                    points.push_back(P0+d1.cartesian().normalized()*a1+d2.cartesian().normalized()*a2);
                    points.push_back(P0+d2.cartesian().normalized()*a2);
                    
                    // make sure it is interstitial
                    VectorDimD Cycle_plane=d1.cartesian().cross(d2.cartesian());
                    VectorDimD b=sessileb[rSS_sessile].cartesian();
                    if (b.dot(Cycle_plane)<0)
                    {
                        b*=-1.0;
                    }
                    const VectorDimD a=(b.normalized());
                    
                    if(allPointsInGrain(points,grainID))
                    {
                        
                        for(int k=0;k<4;++k)
                        {
                            configIO.nodes().emplace_back(nodeID+k,points[k],Eigen::Matrix<double,1,3>::Zero(),1.0,0);
                            
                            const int nextNodeID=(k+1)<4? nodeID+k+1 : nodeID;
                            configIO.loopLinks().emplace_back(loopID,nodeID+k,nextNodeID,true,0);
                            
                        }
                        
                        configIO.loops().emplace_back(loopID+0, b,a,P0,grainID,DislocationLoopIO<MicrostructureGenerator::dim>::SESSILELOOP);
                        
                        
//                        snID++;
                        loopID++;
                        nodeID+=4;
                        ndefects++;
                        defectsDensity=ndefects/mesh.volume()/std::pow(poly.b_SI,3);
                        std::cout<<"irradiation defects density="<<defectsDensity<<std::endl;
                    }
                }
            }
        }
        
        /**********************************************************************/
        void MicrostructureGenerator::addIrradiationLoops()
        {
            if(targetIrradiationLoopDensity>0.0)
            {
                std::cout<<magentaBoldColor<<"Generating Irradiation Loops"<<defaultColor<<std::endl;

                if(poly.crystalStructure=="BCC")
                {
                    addIrradiationLoopsBCC();
                }
                else if(poly.crystalStructure=="FCC")
                {
                    addIrradiationLoopsFCC();
                }
                else if(poly.crystalStructure=="HEX")
                {
                    addIrradiationLoopsHCP();
                }
                else
                {
                    std::cout<<"irradiation loops supported for FCC/BCC/HEX crystals only"<<std::endl;
                }

            }
        }

        /**********************************************************************/
        double MicrostructureGenerator::deltaHelicity(const std::vector<VectorDimD>& newPoints,
                             const VectorDimD& newBurgers) const
        {

            double h(0.0);
            assert(loopPoints.size()==loopBurgers.size());
            for(size_t k=0;k<loopPoints.size(); ++k)
            {
                h+=LinkingNumber<MicrostructureGenerator::dim>::loopPairHelicity(loopPoints[k],loopBurgers[k],newPoints,newBurgers);
            }

            return h;
        }

        /**********************************************************************/
        std::pair<LatticeVector<MicrostructureGenerator::dim>,int> MicrostructureGenerator::randomPointInMesh() const
        {
            return poly.randomLatticePointInMesh();
        }

        /**********************************************************************/
        double MicrostructureGenerator::randomSize()
        {
            std::uniform_real_distribution<double> dist(minSize,maxSize);
            return dist(generator);
        }

        /**********************************************************************/
        int MicrostructureGenerator::randomSign()
        {
            std::uniform_int_distribution<> dis(0,1);
            return  dis(generator)*2-1;
        }

        /**********************************************************************/
        std::map<double,typename MicrostructureGenerator::VectorDimD> MicrostructureGenerator::boundaryProjection(const VectorDimD& P0,
                                                              const VectorDimD& P1,
                                                              const VectorDimD& D,
//                                                              const PlaneMeshIntersectionContainerType& pp
                                                              const MeshBoundaryContainerType& pp)
        {



            const double dNorm(D.norm());
            assert(dNorm>FLT_EPSILON);
            const VectorDimD dir=D/dNorm;
            // Let a point v on the boundary be written as v=P0+u1*(P1-P0)+u2*d
            // then we have [P1-P0 d]*[u1 u2]^T=v-P0

            Eigen::Matrix<double,3,2> A;
            A.col(0)=P1-P0;
            A.col(1)=dir;
            const Eigen::LLT<Eigen::Matrix<double,2,2>> llt(A.transpose()*A);
            assert(llt.info()==Eigen::Success);

            std::map<double,VectorDimD> temp; // keep points sorted by parameter u1
            for(size_t m=0;m<pp.size();++m)
            {
                const Eigen::Matrix<double,2,1> x=llt.solve(A.transpose()*(pp[m]->P0-P0));
                if(x(0)>FLT_EPSILON && x(0)<1.0-FLT_EPSILON && x(1)>FLT_EPSILON)
                {
                    temp.emplace(x(0),pp[m]->P0);
                }

                const Eigen::Matrix<double,2,1> y=llt.solve(A.transpose()*(pp[m]->P1-P0));
                if(y(0)>FLT_EPSILON && y(0)<1.0-FLT_EPSILON && y(1)>FLT_EPSILON)
                {
                    temp.emplace(y(0),pp[m]->P1);
                }

            }

            return temp;


        }

        /**********************************************************************/
        std::pair<int,typename MicrostructureGenerator::VectorDimD> MicrostructureGenerator::boundaryProjection(const VectorDimD& P,
                                                            const VectorDimD& D,
//                                                            const PlaneMeshIntersectionContainerType& pp,
                                                            const MeshBoundaryContainerType& pp)
        {
            const double dNorm(D.norm());
            assert(dNorm>FLT_EPSILON);
            const VectorDimD dir=D/dNorm;
            // line1 is P+u1*dir, u>0

            bool success=false;
            std::pair<int,VectorDimD> temp=std::make_pair(-1,VectorDimD::Zero());

            for(size_t k=0;k<pp.size();++k)
            {
//                const size_t k1 = ((k==(pp.size()-1))? 0 : k+1);
                const VectorDimD& v0=pp[k]->P0;
                const VectorDimD& v1=pp[k]->P1;
                // line2 is v0+u2*(v1-v0), 0<=u2<=1

                //P+u1*dir=v0+u2*(v1-v0)
                // [dir -(v1-v0)] [u1 u2] = [v0-P]
                // In least square sense
                // [dir -(v1-v0)]^T*[dir -(v1-v0)]* [u1 u2] = [v0-P]

                Eigen::Matrix<double,3,2> A;
                A.col(0)=dir;
                A.col(1)=-(v1-v0);

                const Eigen::Matrix<double,3,1> b=v0-P;

                const Eigen::LLT<Eigen::Matrix<double,2,2>> llt(A.transpose()*A);
                //                std::cout<<"DO NOT USE LLT TO SEE IF SYSTEM HAS SOLUTION. See https://eigen.tuxfamily.org/dox/classEigen_1_1LDLT.html#a858dc77b65dd48248299bb6a6a758abf"<<std::endl;


                if(llt.info()==Eigen::Success)
                {

                    const Eigen::Matrix<double,2,1> x=llt.solve(A.transpose()*b);

                    if(x(0)>=0.0 && x(1)>=0.0 && x(1)<=1.0)
                    {
                        success=true;
                        temp=std::make_pair(k,v0+x(1)*(v1-v0));
                        break;
                    }
                }

            }

            assert(success);
            return temp;
        }

}
#endif
