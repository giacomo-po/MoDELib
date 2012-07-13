//==============================================================================
//auxilary tool to generate initial file for the traction and displacement 
// boundary conditions on surface nodes
//
//        Mamdouh Mohamed , Sept 2011, mamdouh.s.mohamed@gmail.com
//============================================================================

#include <vector>
#include <stdio.h>
#include <iostream>

std::vector<std::vector<float> > P;                // matrix nNodes x 3 contains the nodes coordinates
std::vector<std::vector<unsigned int> > Tri;        // matrix nTriangles x 3 contains the indexes of the triangle nodes

    //=======================================================================
    // function to read the mesh nodes file
    //=======================================================================
    
    unsigned int readNodes(){
      
      int nNodes , d1, d2 , d3, in;
      
      FILE *fp =fopen("mesh.node", "r");
			
      fscanf (fp, "%d%d%d%d", &nNodes , &d1, &d2, &d3);      
      P.resize(nNodes, std::vector<float>(3));
	  
      for (unsigned int i = 0; i< nNodes; i++){
	fscanf (fp, "%d%f%f%f", &in , &P[i][0], &P[i][1], &P[i][2]);
      }
      
      fclose(fp); 
      
      return nNodes;
    }
    
    
    //=======================================================================
    // function to read the face mesh file
    //=======================================================================
    
    unsigned int readFaceMesh(){
      
      unsigned int nTris , di , iFc, iTet , ti;
		  
      FILE *fTri =fopen("mesh.face", "r");
		  
      fscanf (fTri, "%d%d", &nTris , &di);
           
      Tri.resize(nTris, std::vector<unsigned int>(3));
	  
      for (unsigned int i = 0; i< nTris; i++){
	fscanf (fTri, "%d%d%d%d%d%d", &ti , &Tri[i][0], &Tri[i][1], &Tri[i][2], &iTet  , &iFc );
      }
      
      fclose(fTri); 
      
      return nTris;
    }

//=======================================================================
//=======================================================================
int main (int argc, char * const argv[]) {

    unsigned int nNodes = readNodes();   
    
    unsigned int nTris = readFaceMesh();
    
    //---------------- set BCs for selected nodes --------------------
    FILE *fout =fopen("BCs_0.txt", "w");
    
    for (unsigned int i=0; i<nNodes; i++){
      if (P[i][2] == 0.0) {
	//  node_index , bool u_x , bool u_y, bool u_z , u_x  , u_y , u_z , traction_x , traction_y , traction_z
	fprintf (fout, "%8d %2d %2d %2d %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f \n", i , 1, 1, 1 , 0.0  , 0.0 , 0.0 , 0.0  , 0.0 , 0.0);
      }
      /*
      else if (P[i][1] == 0.0) {
	//  node_index , bool u_x , bool u_y, bool u_z , u_x  , u_y , u_z , traction_x , traction_y , traction_z
	fprintf (fout, "%8d %2d %2d %2d %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f \n", i , 0, 1, 0 , 0.0  , 0.0 , 0.0 , 0.0  , 0.0 , 0.0);
      }
      
      else if (P[i][2] == 0.0) {
	//  node_index , bool u_x , bool u_y, bool u_z , u_x  , u_y , u_z , traction_x , traction_y , traction_z
	fprintf (fout, "%8d %2d %2d %2d %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f \n", i , 0, 0, 1 , 0.0  , 0.0 , 0.0 , 0.0  , 0.0 , 0.0);
      }
      */
      //else if (P[i][1] == 1000.0){
	//fprintf (fout, "%8d %2d %2d %2d %15.8f %15.8f %15.8f %15.8f %15.8f %15.8f \n", i , 0, 0, 0 , 0.0  , 0.0 , 0.0 , -0.06 , 0.0 , -0.06);
    // }
      
    }
        
    fclose(fout); 
    
    return 0;
} 



    