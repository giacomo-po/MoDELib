
//double diameter = 2127;
//double zSurf = 3*diameter;
//double load = 0.0;
//double disp = 0.0;
//size_t counter =0;
//
//Eigen::Matrix<double,dim,dim> temp(DN.plasticStrainRate());
//
//for (typename Domain::NodeContainerType::const_iterator iter=DN.shared.domain.nodeContainer.begin();iter!=DN.shared.domain.nodeContainer.end();++iter){
//    if (!iter->isBoundaryNode) continue;
//    if ( iter->P(2)!=zSurf) continue;
//    //    std::cout<<"u = "<<iter->u.transpose()+iter->uInf.transpose()<<std::endl;
//    Eigen::Matrix<double,dim,1> UT( iter->u.transpose() + iter->uInf.transpose() );
//    counter++;
//    disp +=  UT(2);
//}
//
//
//for (unsigned int it=0; it < DN.shared.domain.triContainer.size(); it++ ){
//    bvpfe::Triangle* pTri = DN.shared.domain.triContainer[it];
//    if (pTri->eleNodes[0]->P(2)!=zSurf || pTri->eleNodes[1]->P(2)!=zSurf || pTri->eleNodes[2]->P(2)!=zSurf) continue;
//    unsigned int iTet = pTri->neighTetIndx;
//    Eigen::Matrix<double,dim,1> forceVec(pTri->area() * ( DN.shared.domain.tetContainer[iTet].getStress() * pTri->outNormal + pTri->getTriInfiniteTraction<3>(&DN) ))  ;
//    load +=  forceVec(2);
//}
//
//
//UniqueOutputFile<'S'> standard_output;
//standard_output<<counter<<"  "<<disp<<"  "<<load<<"  "<<temp(0,0)<<"  "<<temp(0,1)<<"  "<<temp(0,2)<<"  "<<temp(1,1)<<"  "<<temp(1,2)
//<<"  "<<temp(2,2)<<"  "<<DN.get_dt()<<std::endl;
//
//
//
//SequentialOutputFile<'D',1>::set_increment(outputFrequency); // Displacement_file;
//SequentialOutputFile<'D',1>::set_count(runID); // Displacement_file;
//SequentialOutputFile<'D',1> Dfile;
//for (typename Domain::NodeContainerType::const_iterator iter=DN.shared.domain.nodeContainer.begin();iter!=DN.shared.domain.nodeContainer.end();++iter){
//    if (iter->isBoundaryNode){
//        Dfile<< iter->sID <<" "<< (iter->u + iter->uInf).transpose()<<"\n";
//    }
//}


/*
 // Example of a SequentialOutputFile that is saved in P/P_n.txt (n is the timestep) and to which we append during a single time step.
 SequentialOutputFile<'P',1>::set_increment(outputFrequency); // Edges_file;
 SequentialOutputFile<'P',1>::set_count(runID); // Edges_file;
 SequentialOutputFile<'P',1> p_file;
 // let's output the node velocity
 for (typename NetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
 p_file << nodeIter->second->sID<<"\t"<< std::setprecision(15)<<std::scientific<<nodeIter->second->velocity.transpose()<<"\n";
 }
 */
