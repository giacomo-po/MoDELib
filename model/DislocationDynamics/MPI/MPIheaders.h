#ifdef HAVE_MPI
#include <mpi.h>
#endif

///////////////////////////////////////////////////////////////////////////
// Parallel
//Include Isorropia_Exception.hpp only because the helper functions at
//the bottom of this file (which create the epetra objects) can
//potentially throw exceptions.
#include <Isorropia_Exception.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

//The Isorropia user-interface functions being demonstrated are declared
//in Isorropia_Epetra.hpp.
#include <Isorropia_Epetra.hpp>

#include <Teuchos_ParameterList.hpp>

#ifdef HAVE_EPETRA
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "Epetra_Comm.h"
#include <Epetra_Map.h>
#include "Epetra_BlockMap.h"
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#endif

#include "ispatest_lbeval_utils.hpp"

//Declarations for helper-functions that create epetra objects. These
//functions are implemented at the bottom of this file.

inline int edgeID(int i, int j, int nx) {
    return i + j*nx;
}

inline int edgeID(int i, int j, int k, int nx, int ny) {
    return i + j*nx + k*nx*ny;
}

#ifdef HAVE_EPETRA
Teuchos::RCP<Epetra_MultiVector>
create_epetra_multivec(const Epetra_Map &rowmap, 
					   int numProcs, int localProc, int dim);

Teuchos::RCP<Epetra_CrsMatrix>
create_epetra_matrix(const Epetra_MultiVector &pos, 
					 const Epetra_Map &rowmap, 
					 int numProcs, int localProc);

Teuchos::RCP<Epetra_Map> 
create_epetra_mapcmp(const Epetra_CrsMatrix &crsmatrix, 
					 std::map<int, std::set<int> > &gVtxInCol);

Teuchos::RCP<Epetra_MultiVector>
nn_epetra_solver(const Epetra_MultiVector &posCmp, 
				 const Epetra_Map &rowmapCmp, 
				 std::map<int, std::set<int> > &gVtxInCol, 
				 const Epetra_CrsMatrix & crsmatrix);
#endif

///////////////////////////////////////////////////////////////////////////