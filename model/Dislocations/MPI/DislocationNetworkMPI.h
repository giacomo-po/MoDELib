
public:
//////////////////////////////////////////////////////////////
// Non-Default Constructor
DislocationNetwork(int argc, char** argv){
	
	
	std::cout<<"MPI Constructor"<<std::endl;
	
	EpetraCols=3*dim+2;
	
	//first, set up our MPI environment...
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &localProc);	// get the local processor ID
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);	// get the total N of processor
	
}

//////////////////////////////////////////////////////////////
// Destructor for MPI
~DislocationNetwork(){
	
	MPI_Finalize();
	
}


private:

int numProcs;
int localProc;
int EpetraCols; // dimensions + num of variable quantities needed

void MPIstep(){
	
	
	
	
	//create an Epetra_Map with rows spread un-evenly over
	//processors.
	Epetra_MpiComm comm(MPI_COMM_WORLD);
	
	
	//NetworkBaseType::speak<int>(3);
	
	//this->NaddressesX<ootndd::DislocationSegment<dim,corder,alpha,InterpolationType,MaterialType>	>();
	//this->ABbegin<LinkType>();
	//this->ABend<LinkType>();
	
	//step();
	
	int global_num_rows = this->linkOrder()*qOrder;
	int num_procs_res = global_num_rows%numProcs;
	int local_num_rows = global_num_rows/numProcs;
	
	if ( (num_procs_res!=0) && (localProc==(numProcs-1)) )
		local_num_rows += num_procs_res;
	
	//now we're ready to create a row-map.
	//Teuchos::RCP<Epetra_Map> rowmap = Teuchos::rcp(new Epetra_Map(global_num_rows, local_num_rows, 0, comm));
	Teuchos::RCP<Epetra_Map> rowMap = Teuchos::rcp(new Epetra_Map(global_num_rows, 0, comm));
	
	
	if (localProc == 0) {
		std::cout << " creating Epetra_MultiVector with un-even distribution..."
		<< std::endl;
	}
	
	/*
	 //create an Epetra_CrsGraph with rows spread un-evenly over
	 //processors.
	 Epetra_MpiComm comm(MPI_COMM_WORLD);
	 int global_num_rows = 26;
	 int num_procs_res = global_num_rows%numProcs;
	 int local_num_rows = global_num_rows/numProcs;
	 
	 if ( (num_procs_res!=0) && (localProc==(numProcs-1)) )
	 local_num_rows += num_procs_res;
	 
	 //now we're ready to create a row-map.
	 Epetra_Map rowmap(global_num_rows, local_num_rows, 0, comm);
	 */
	
	
	//now we create a group of position points.
	//const int dim = 2;
	//	const double cod = 5.0; // cut-off distance
	//	const std::vector<double> len( dim, 20.0 );// len[0]=lx, len[1]=ly, len[2]=lz
	//	const std::vector<int> nne ( dim, static_cast<int>(len[0]/cod) ); // nne[0]=nx, nne[1]=ny, nne[2]=nz
	//	const std::vector<bool> pbcs( dim, false ); // periodic boundary conditions for the block?
	
	//! 1 - Create a Multi Vector
	Teuchos::RCP<Epetra_MultiVector> quadraturePointsData = Teuchos::rcp(new Epetra_MultiVector(*rowMap, EpetraCols));
	fillMultiVector(rowMap,quadraturePointsData);


	
	
	std::cout<<"THIS IS THE DATA"<<std::endl<<*quadraturePointsData<<std::endl;
	
	
	//! 2 - Create a Hyper Graph Matrix
	Teuchos::RCP<Epetra_CrsMatrix> hypGraphMatrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *rowMap, 0));
	
	fillHyperGraphMatrix(rowMap,quadraturePointsData,hypGraphMatrix);
	std::cout <<"THIS IS THE HYPER"<<std::endl<< *hypGraphMatrix<<std::endl;
	
	
	//! 3 - Create gVtxInCol and rowmapCmp
	std::map<int, std::set<int> > gVtxInCol;
	
	Teuchos::RCP<Epetra_Map> rowmapCmp;
	fillRowmapCmp(*hypGraphMatrix, gVtxInCol,rowmapCmp);
	
	//! 4 - create the relationship between rowmapCmp (workingMap) and rowMap (original map)
	Teuchos::RCP<Epetra_Import> imp = Teuchos::rcp(new Epetra_Import(*rowmapCmp, *rowMap));
	
	Teuchos::RCP<Epetra_MultiVector> dataCmp = Teuchos::rcp(new Epetra_MultiVector(*rowmapCmp, EpetraCols));
	dataCmp->Import(*quadraturePointsData, *imp, Add);
	
	std::cout<<"THIS IS THE COMP DATA"<<std::endl<<*dataCmp<<std::endl;	
	
	
	
}



/////////////////////////////////////////////////////////
// fillMultiVector
void fillMultiVector(const Teuchos::RCP<Epetra_Map> & rowMap, 
/*      output    */ Teuchos::RCP<Epetra_MultiVector> & quadraturePointsData){
	
	int global_row;
	int sequentialLinkID;
	int localQuadratureID;
	
	linkIterType linkIter;// = this->linkBegin(); 
	
	for (int i=0; i<quadraturePointsData->MyLength(); ++i) { 
		
		global_row = rowMap->GID(i);
		
		//		
		linkIter = this->linkBegin();			// FIND A SMARTER WAY!!!!!
		sequentialLinkID=global_row/qOrder;
		for (int k=0;k<sequentialLinkID;++k){
			++linkIter;
		}
		
		
		localQuadratureID=global_row%qOrder;
	//	std::cout<<"Processor ID ="<<localProc
	//	<<". link sID ="<<linkIter->second->sID
	//	<<" .localQuadratureID = "<<localQuadratureID
		//<<" rgauss "<<linkIter->second->get_rgauss()
	//	<<std::endl;
		
		for (int d=0;d<dim;++d){
			quadraturePointsData->ReplaceGlobalValue(global_row,d,		linkIter->second->get_rgauss (d,localQuadratureID));	// Rx, Ry, Rz
			quadraturePointsData->ReplaceGlobalValue(global_row,d+dim,	linkIter->second->get_rugauss(d,localQuadratureID));	// Tx, Ty, Tz
			quadraturePointsData->ReplaceGlobalValue(global_row,d+2*dim,linkIter->second->get_Burgers(d                  ));	// Bx, By, Bz
		}
		
		quadraturePointsData->ReplaceGlobalValue(global_row,3*dim,	linkIter->second->get_jgauss(localQuadratureID));	// jacobian
		quadraturePointsData->ReplaceGlobalValue(global_row,3*dim+1,linkIter->second->weight	(localQuadratureID));	// quadrature weight
		
	}


}




/////////////////////////////////////////////////////////
// fillHyperGraphMatrix
void fillHyperGraphMatrix(const Teuchos::RCP<Epetra_Map> & rowMap, 
/*      input          */ const Teuchos::RCP<Epetra_MultiVector> & quadraturePointsData, 
/*      output         */ Teuchos::RCP<Epetra_CrsMatrix> & hypGraphMatrix){
	
	double* systemCoordinateMax = new double[EpetraCols]; 
	double* systemCoordinateMin = new double[EpetraCols]; 
	
	
	// Find the center of the simulation block
	quadraturePointsData->MaxValue(systemCoordinateMax);
	quadraturePointsData->MinValue(systemCoordinateMin);
	
	VectorDimD coordinateMax = VectorDimD::Map(&systemCoordinateMax[0],dim);
	VectorDimD coordinateMin = VectorDimD::Map(&systemCoordinateMin[0],dim);
	
	
	
	double cutOffDistance = 10.0;
	
	ootndd::SpaceDecomposition<dim> SD(coordinateMin,coordinateMax,cutOffDistance);
	
	

	
	int err = 0;
	
	
	VectorDimD pos;
	std::vector<int>  closedNeighborIDs;
	std::vector<double>  coeff;
	
	int global_row;
	for (int i=0; i<quadraturePointsData->MyLength(); ++i) { 
		global_row = rowMap->GID(i);
		
		
		for (int d=0;d<dim;++d){
			pos(d)=(*quadraturePointsData)[d][i];		// you access the multivector by columns (each column is a vector)
		}

		SD.neighborEdgeIDs(pos-coordinateMin,closedNeighborIDs,coeff);

		
		err = hypGraphMatrix->InsertGlobalValues(global_row, closedNeighborIDs.size(),
                                         &coeff[0], &closedNeighborIDs[0]);
		if (err < 0) {
			err = hypGraphMatrix->ReplaceGlobalValues(global_row, closedNeighborIDs.size(),
											  &coeff[0], &closedNeighborIDs[0]);
			if (err < 0) {
				throw Isorropia::Exception("create_epetra_matrix: error inserting matrix values.");
			}
		}
		
		
		
		
	}

	
	err = hypGraphMatrix->FillComplete();
	if (err != 0) {
		throw Isorropia::Exception("create_epetra_matrix: error in matrix.FillComplete()");
	}
	
	
}



/////////////////////////////////////////////////////////
// fillRowmapCmp
void fillRowmapCmp(const Epetra_CrsMatrix &crsmat, std::map<int, std::set<int> > &gVtxInCol, Teuchos::RCP<Epetra_Map> & rowmapCmp) // save global vertices in each edge 
{
	//    const Epetra_Map &rowmap = crsmat.RowMap();
    const Epetra_Map &colmap = crsmat.ColMap();
	//
    const Epetra_Comm &comm = crsmat.Comm();
#ifdef HAVE_MPI
    const Epetra_MpiComm* mpiComm = dynamic_cast<const Epetra_MpiComm*>(&comm);
    MPI_Comm mcomm = mpiComm->Comm();
    //MPI::Intracomm mcomm = mpiComm->Comm();
#endif
	
	
	//    int nProcs      = comm.NumProc();
	//    int myProc      = comm.MyPID();
	//    int numMyRows   = rowmap.NumMyElements();
	//    int numMyCols   = colmap.NumMyElements();
	
    //int nProcs = crsmat.Comm().NumProc();
    int nProcs = comm.NumProc();
    int myProc = comm.MyPID();
    int numMyRows = crsmat.NumMyRows();
    int numNzRow = crsmat.MaxNumEntries();
    int numGzCol = crsmat.NumGlobalCols(); 
	
	
    if (myProc == 0) {
        std::cout << " creating Epetra_Map for computation with overlapping vertices..."
		<< std::endl;
    }
	
    int **vtxpe;
    AllocateArray(vtxpe,nProcs,numGzCol);
    int **srpe;
    AllocateArray(srpe, nProcs,numGzCol); // note: NumGlobalCols = NumGlobalRows, 
	//       should be optimized... !!!!!
    for (int i=0; i<numGzCol; ++i) {
        srpe[myProc][i] = 0;
        vtxpe[myProc][i] = 0;
    }
	
	
    std::map<int, std::vector<int> > vtxInCol;
    std::map<int, std::vector<int> >::iterator gidIter;
    std::map<int, std::set<int> >::iterator gidIterSet;
	
    int* maxCol = new int[nProcs];
    maxCol[myProc]=0;
	
    for (int i=0; i<numMyRows; ++i) {
        int gidRow   = crsmat.GRID(i);
        int *indices = new int[numNzRow];
        double *vals = new double[numNzRow];
		
        int numEntries;
        crsmat.ExtractGlobalRowCopy(gidRow, numNzRow, numEntries, vals, indices); 
		
        for (int j=0; j<numEntries; ++j) {
            int gidCol = indices[j]; //colmap.GID(j);
            if (gidCol>maxCol[myProc]) maxCol[myProc]=gidCol; // test..
            gidIter = vtxInCol.find(gidCol);
            if (gidIter != vtxInCol.end()) {
                vtxInCol[gidCol].push_back(gidRow);
            } else {
                std::vector<int> tmp;
                tmp.push_back(gidRow);
                vtxInCol.insert(make_pair(gidCol,tmp));
            }
            if (vals[j]==1.0) {
                srpe[myProc][gidCol] = 1;
            } else if (vals[j]==-1.0 && srpe[myProc][gidCol]!=1) {
                srpe[myProc][gidCol] = -1;
            } 
        }
		
        delete[] indices;
        delete[] vals;
    }
	
    for (int iEdge = 0; iEdge<numGzCol; ++iEdge) {
        gidIter = vtxInCol.find(iEdge);
        if (gidIter != vtxInCol.end()) {
            vtxpe[myProc][iEdge] = gidIter->second.size();
        } else {
            vtxpe[myProc][iEdge] = 0;
        }
    }
	
    MPI_Allgather(&vtxpe[myProc][0], numGzCol, MPI_INT, 
                  &vtxpe[0][0],      numGzCol, MPI_INT, mcomm);
    MPI_Allgather(&srpe[myProc][0],  numGzCol, MPI_INT, 
                  &srpe[0][0],       numGzCol, MPI_INT, mcomm);
	
    std::set<int> myGlobalElements;
    MPI_Status  status;
	
    for (int iEdge=0; iEdge<numGzCol; ++iEdge) {
		//for (int iEdge=0; iEdge<6; ++iEdge) {
        //MPI_Request** r;
        //AllocateArray(r, nProcs, nProcs);
        //std::vector<MPI_Request> r;
		
        //int sumVtx(0);
        //for (int ip=0; ip<nProcs; ++ip) {
        //    sumVtx += vtxpe[ip][iEdge];
        //}
        //myGlobalElements.resize(sumVtx);
		
        std::set<int> rProc;
        std::set<int> sProc;
        std::set<int>::iterator pIter;
        for (int ip=0; ip<nProcs; ++ip) {
            if ( srpe[ip][iEdge]==1 ) {
                rProc.insert(ip);
            } else if (srpe[ip][iEdge]==-1 ) {
                sProc.insert(ip);
            }
        }
        if (rProc.size()>1) {
            for (int ip=0; ip<nProcs; ++ip) {
                if ( srpe[ip][iEdge]==1 ) {
                    sProc.insert(ip);
                }
            }
        }
		
        int nReq(0);
        MPI_Request *rr;
		
        pIter = sProc.find(myProc);
        if ( pIter != sProc.end() ) {
            nReq = rProc.size();
            rr=new MPI_Request[nReq];
            int ir=0;
            for (std::set<int>::iterator it=rProc.begin(); it!=rProc.end(); ++it) {
                int pDest = *it;
                int pSour = *pIter;
				
                if ( pDest!=pSour ) {
					
                    MPI_Isend(&vtxInCol[iEdge][0], vtxpe[pSour][iEdge], MPI_INT, pDest, 0, mcomm, &rr[ir]);
                    ++ir;
                    //MPI_Send(&vtxInCol[iEdge][0], vtxpe[pSour][iEdge], MPI_INT, pDest, 0, mcomm);
                } else {
                    --nReq;
                }
            }
			
        }
		
        pIter = rProc.find(myProc);
        if ( pIter != rProc.end() ){
            //int ige(0);
            for (int iv=0; iv<vtxpe[*pIter][iEdge]; ++iv) {
                myGlobalElements.insert( vtxInCol[iEdge][iv] );
                //myGlobalElements[ige]= vtxInCol[iEdge][iv] ;
                //++ige;
				
                int gidCol = iEdge;
                gidIterSet = gVtxInCol.find(gidCol);
                if ( gidIterSet != gVtxInCol.end() ) {
                    gVtxInCol[gidCol].insert( vtxInCol[iEdge][iv] );
                } else {
                    std::set<int> tmp;
                    tmp.insert( vtxInCol[iEdge][iv] );
                    gVtxInCol.insert(std::make_pair(gidCol,tmp));
                }
                
            }
            //nReq = sProc.size();
            //r=new MPI_Request[nReq];
            //int ir=0;
            for (std::set<int>::iterator it=sProc.begin(); it!= sProc.end(); ++it) {
                int pSour = *it;
                int pDest = *pIter;
				
                if ( pSour!=pDest ) {
					
                    int *vtmp = new int[ vtxpe[*it][iEdge] ];
                    //MPI_Irecv(&myGlobalElements[ige], vtxpe[*it][iEdge], MPI_INT, pSour, 0, mcomm, r);
                    //++ir;
                    MPI_Recv(vtmp, vtxpe[pSour][iEdge], MPI_INT, pSour, 0, mcomm, &status);
					
                    for (int iv=0; iv<vtxpe[*it][iEdge]; ++iv) {
                        myGlobalElements.insert( vtmp[iv] );
						
                        int gidCol = iEdge;
                        gidIterSet = gVtxInCol.find(gidCol);
                        if ( gidIterSet != gVtxInCol.end() ) {
                            gVtxInCol[gidCol].insert( vtmp[iv] );
                        } else {
                            std::set<int> tmp;
                            tmp.insert( vtmp[iv] );
                            gVtxInCol.insert(std::make_pair(gidCol,tmp));
                        }
						
                    }
                    delete [] vtmp;
                }
            }
        }
        
        
        if (nReq>0){
            //MPI_Wait(&r, MPI_STATUSES_IGNORE);
            MPI_Waitall(nReq, rr, MPI_STATUSES_IGNORE);
            delete [] rr;
        }
        
    }
    
    /*
	 if (myProc==0) {
	 std::cout << "myGlobalElements = ";
	 for (std::set<int>::iterator it=myGlobalElements.begin(); it!=myGlobalElements.end(); ++it) {
	 std::cout << *it << ", ";
	 }
	 std::cout<<endl;
	 std::cout << "myGlobalElements in col[4] =";
	 for (std::set<int>::iterator it=gVtxInCol[4].begin(); it!=gVtxInCol[4].end(); ++it) {
	 std::cout << *it << ", ";
	 }
	 std::cout<<endl;
	 }
	 
     
	 if (myProc==0) {
	 std::cout << "vtxpe: ";
	 for (int i=0; i<numGzCol; ++i) 
	 std::cout << vtxpe[1][i] << " ";
	 std::cout << endl;
	 std::cout << "srpe: ";
	 for (int i=0; i<numGzCol; ++i) 
	 std::cout << srpe[1][i] << " ";
	 std::cout << endl;
	 }
	 */
	
    FreeArray(vtxpe);
    FreeArray(srpe);
	
	
    int numMyElements(myGlobalElements.size());
    int *myGlobalElementsArray = new int[numMyElements];
	
    int ige(0);
    for (std::set<int>::iterator it=myGlobalElements.begin(); it!=myGlobalElements.end(); ++it) {
        myGlobalElementsArray[ige++] = *it;
    }
	
    rowmapCmp = Teuchos::rcp(new Epetra_Map(numMyElements, numMyElements, myGlobalElementsArray, 0, comm));
	
    /*
	 if (myProc !=4) {
	 std::cout << "Processor " << myProc << " has columns: ";
	 for (int i=0; i<rmap->NumMyElements(); ++i) {
	 std::cout << rmap->MyGlobalElements()[i] << " ";
	 }
	 std::cout << endl;
	 std::cout << "Processor " << myProc << " has elementsize= " << rmap->ElementSize() << endl;
	 }
	 */
	
	
    
}




private:







