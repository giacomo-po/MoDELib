

double zSurf = 16.0e03;
double load = 0.0;
double disp = 0.0;
size_t counter =0;

Eigen::Matrix<double,dim,dim> PSR(plasticStrainRate());    

for (typename bvpfe::Domain::NodeContainerType::const_iterator iter=shared.domain.nodeContainer.begin();iter!=shared.domain.nodeContainer.end();++iter){
    if (!iter->isBoundaryNode) continue;
    if ( iter->P(2)!=zSurf) continue;
    //    std::cout<<"u = "<<iter->u.transpose()+iter->uInf.transpose()<<std::endl;
    Eigen::Matrix<double,dim,1> UT( iter->u.transpose() + iter->uInf.transpose() );
    counter++;
    disp +=  UT(2);
}


for (unsigned int it=0; it < shared.domain.triContainer.size(); it++ ){
    bvpfe::Triangle* pTri = shared.domain.triContainer[it];
    if (pTri->eleNodes[0]->P(2)!=zSurf || pTri->eleNodes[1]->P(2)!=zSurf || pTri->eleNodes[2]->P(2)!=zSurf) continue;
    unsigned int iTet = pTri->neighTetIndx;
    Eigen::Matrix<double,dim,1> force_correction(pTri->area() * ( shared.domain.tetContainer[iTet].getStress() * pTri->outNormal ));
    Eigen::Matrix<double,dim,1> force_infinite(pTri->forceInfinite<4,false,DislocationNetworkType>(this));

    load +=  force_correction(2)+force_infinite(2);        
}




UniqueOutputFile<'S'> standard_output;
standard_output<<counter<<"  "<<disp<<"  "<<load<<"  "<<PSR(0,0)<<"  "<<PSR(0,1)<<"  "<<PSR(0,2)<<"  "<<PSR(1,1)<<"  "<<PSR(1,2)
               <<"  "<<PSR(2,2)<<"  "<<dt<<std::endl;



SequentialOutputFile<'D',1>::set_increment(outputFrequency); // Displacement_file;
SequentialOutputFile<'D',1>::set_count(runID); // Displacement_file;
SequentialOutputFile<'D',1> Dfile;
for (typename bvpfe::Domain::NodeContainerType::const_iterator iter=shared.domain.nodeContainer.begin();iter!=shared.domain.nodeContainer.end();++iter){
    if (iter->isBoundaryNode){		
	Dfile<< iter->sID <<" "<< (iter->u + iter->uInf).transpose()<<"\n";		
	} 
}


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
