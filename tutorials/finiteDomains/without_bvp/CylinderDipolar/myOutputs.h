// Example of a UniqueOutputFile that is saved in S/S_0.txt and to which we append data at each output event.
//UniqueOutputFile<'S'> standard_output;
//standard_output<< runID << "  "  << "  " << network_length()*<<std::endl;

// Example of a SequentialOutputFile that is saved in P/P_n.txt (n is the timestep) and to which we append during a single time step.
SequentialOutputFile<'P',1>::set_increment(outputFrequency); // Edges_file;
SequentialOutputFile<'P',1>::set_count(runID); // Edges_file;
SequentialOutputFile<'P',1> p_file;
// let's output the node velocity
for (typename NetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){				
	p_file << nodeIter->second->sID<<"\t"<< std::setprecision(15)<<std::scientific<<(nodeIter->second->get_V().norm())<<"\n";	
}
