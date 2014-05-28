/* This file is part of model, the Mechanics of Defects Evolution Library.
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>,
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_VIRTUALBOUNDARYSLIPCONTAINER_H_
#define  model_VIRTUALBOUNDARYSLIPCONTAINER_H_

#include <boost/ptr_container/ptr_vector.hpp>
#include <Eigen/Dense>

#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/Utilities/CompareVectorsByComponent.h>
#include <model/DislocationDynamics/VirtualBoundarySlipSurface.h>
#include <model/Utilities/SequentialOutputFile.h>
#include <model/DislocationDynamics/DislocationSharedObjects.h>

#include <model/BVP/Domain.h>

#include <stdio.h>
#include <string>


namespace model {
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename DislocationSegmentType>
    class VirtualBoundarySlipContainer : public boost::ptr_vector<VirtualBoundarySlipSurface<DislocationSegmentType> > {
        
        enum{dim=TypeTraits<DislocationSegmentType>::dim};
        
        typedef VirtualBoundarySlipSurface<DislocationSegmentType> VirtualBoundarySlipSurfaceType;
        
        
        typedef boost::ptr_vector<VirtualBoundarySlipSurfaceType> BaseContainerType;
        
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
        //    typedef std::map<Eigen::Matrix<double,dim+1,1> , std::map<Eigen::Matrix<double,dim,1>, std::auto_ptr<VirtualBoundarySlipSurfaceType > , model::CompareVectorsByComponent<double,dim,float> > ,
        //					               model::CompareVectorsByComponent<double,dim+1,float> > radialSegmentsContainerType;
        
        typedef std::map< Eigen::Matrix<double,dim+1,1> , std::vector<VirtualBoundarySlipSurfaceType* > ,
        model::CompareVectorsByComponent<double,dim+1,float>,
        Eigen::aligned_allocator<std::pair<const Eigen::Matrix<double,dim+1,1>,std::vector<VirtualBoundarySlipSurfaceType* > > > > radialSegmentsContainerType;
        
        typedef typename radialSegmentsContainerType::iterator radialSegmentsContainerIterator;
        typedef typename radialSegmentsContainerType::const_iterator radialSegmentsContainerConstIterator;
        
        //------------ map, with the key to be the glide-plane-key, that has pointers for the VirtualBoundarySlipSurface entities that has radial segments with non-zero Burgers (used for stress calculations)
        radialSegmentsContainerType  radialSegmentsContainer;
        
    public:
        
        //=======================================================================
        // add a new entity to the container
        //======================================================================
        //      template <typename DislocationSegmentType>
        void add (const DislocationSegmentType& ds) {
			
            std::auto_ptr<VirtualBoundarySlipSurfaceType > pVBSS (new VirtualBoundarySlipSurfaceType (ds) );
            this->push_back(pVBSS);
            
            //---------- check radial segments that are required for stress calculation
            const Eigen::Matrix<double,dim+1,1> gpKey((Eigen::Matrix<double,dim+1,1>()<<ds.pGlidePlane->planeNormal.normalized() , ds.pGlidePlane->height).finished());
            //gpKey << ds.pGlidePlane->planeNormal.normalized() , ds.pGlidePlane->height;
            
            radialSegmentsContainerIterator gpIter = radialSegmentsContainer.find(gpKey);
            
            if (gpIter == radialSegmentsContainer.end()) {     // no segments from the same glide plane were stored before, so just push_back
                
                std::vector<VirtualBoundarySlipSurfaceType* > tempV;   tempV.push_back(&(*(this->rbegin())));
                radialSegmentsContainer.insert(std::make_pair (gpKey , tempV ) );
                
                this->rbegin()->addRadialSegments ();
            }
            
            //-- else, check the segments repeatness --------
            else {
                
                //--------- WARNING : removing the bool variables, and placing the functions call directly inside the if-statement
                //--------- may cause to have only the first function excuted if it returns TRUE, and so the call for the 2nd function will be neglected
                bool source_repeated = checkSegmentsRepeatness (gpIter, ds.source->get_P(), ds.Burgers, 0);
                bool sink_repeated   = checkSegmentsRepeatness (gpIter, ds.sink->get_P()  , ds.Burgers, 1);
                if ( source_repeated || sink_repeated )  gpIter->second.push_back(&(*(this->rbegin())));
                
            }
            
            /*
             for (gpIter = radialSegmentsContainer.begin(); gpIter!=radialSegmentsContainer.end(); gpIter++) {
             std::cout << "Container size " << gpIter->second.size() << std::endl;
             }
             std::cout << std::endl;
             */
            
            
        }
        
        
        
        /* initializeVirtualSegments ********************************************/
        template<typename DislocationNetworkType>
        void initializeVirtualSegments (DislocationNetworkType& DN) {
            
            VectorDim sourceP, sinkP, Burgers;
            
            std::stringstream filename;
            filename << "B/B_" << DN.runningID() << ".txt";
            
            unsigned int sourceID, sinkID;
            unsigned int ii=0;
            
            FILE *fp =fopen(filename.str().c_str(), "r");
            
            if (fp != NULL) {
                while (!feof(fp)) {
                    if(fscanf(fp, "%le%le%le%le%le%le%le%le%le", &sourceP(0),&sourceP(1),&sourceP(2),&sinkP(0),&sinkP(1),&sinkP(2),&Burgers(0),&Burgers(1),&Burgers(2))==9){
                        ii++;
                        
                        sourceID = DN.insertVertex(sourceP);
                        if (DN.node(sourceID).second->meshLocation()!=onMeshBoundary) {
                            std::cout << "Error: Source node of Virtual dislocation no. " << ii << " was not recognized as boundary node. Check its coordinates" << std::endl;
                            assert(0);
                        }
                        
                        sinkID   = DN.insertVertex(sinkP);
                        
                        if (DN.node(sinkID).second->meshLocation()!=onMeshBoundary) {
                            std::cout << "Error: Sink node of Virtual dislocation no. " << ii << " was not recognized as boundary node. Check its coordinates" << std::endl;
                            assert(0);
                        }
                        if ((sinkP-sourceP).norm()>1.0e-5)
                        {
                            DN.connect(sourceID,sinkID,Burgers); // create a dislocation segment
                            DN.template disconnect<true>(sourceID,sinkID); // destroy the dislocation segment in order to create the boundary segment. true=remove isolated nodes
                        }
                    }
                }
                
                fclose(fp);
            }
            
        }
        
        
        
        
        
        //==========================================================================
        // function to check and consider the radial segments repeatness
        //==========================================================================
        
        bool checkSegmentsRepeatness (const radialSegmentsContainerIterator& gpIter, const VectorDim& P, const VectorDim& disBurg, const unsigned int& addIndex) {
            
            double tol = 1.0e-7;
            
            int toBeRemoved_r = -1;
            int toBeRemoved_v = -1;
            
            int sign = 1;      if (addIndex==1) sign = -1;         // sign of the Burgers vector at this node (if source -> +ve, if sink -> -ve)
            
            //---------- loop over all virtual segments on this glide plane ----
            for (unsigned int iv=0; iv<gpIter->second.size(); iv++) {
                
                //--------------- loop over all existing radial segments that belong to this virtual segment ---
                for (unsigned int ir=0; ir<gpIter->second[iv]->radialSegmentsVector.size(); ir++) {
                    
                    if( (P-gpIter->second[iv]->radialSegmentsVector[ir].sourceP).norm() < tol && (sign*disBurg + gpIter->second[iv]->radialSegmentsVector[ir].Burg).norm() < tol ) {
                        //--- segment is repeated, so remove it from "gpIter->second[iv]->radialSegmentsVector"
                        toBeRemoved_r = ir;
                        toBeRemoved_v = iv;
                        break;
                    }
                }
                
                if (toBeRemoved_r >= 0) break;
            }
            
            if (toBeRemoved_r >= 0) {
                gpIter->second[toBeRemoved_v]->radialSegmentsVector.erase(gpIter->second[toBeRemoved_v]->radialSegmentsVector.begin()+toBeRemoved_r);
                //std::cout<< "Repeated " << addIndex << "  " << gpIter->second[toBeRemoved_v]->radialSegmentsVector.size() << std::endl;
                
                //------- if a virtual dislocation didn't have any more radial segments, remove it from the list "radialSegmentsContainer"
                if (gpIter->second[toBeRemoved_v]->radialSegmentsVector.size() == 0) gpIter->second.erase(gpIter->second.begin()+toBeRemoved_v);
            }
            else	this->rbegin()->addRadialSegments (addIndex);
            
            return (toBeRemoved_r == -1);
        }
        
        /* stress *************************************************************/
        MatrixDim stress(const VectorDim& Rfield) const {
            /*! the stress field induced by boundary segments stored in this container
             */
            MatrixDim temp(MatrixDim::Zero());
            for(typename BaseContainerType::const_iterator sIter=this->begin();sIter!=this->end();++sIter){
                if(sIter->radialSegmentsVector.size()){
                    temp+=sIter->stress(Rfield);
                }
            }
            
            return temp;
        }
        
        
        /* stressFromGlidePlane ***********************************************/
        MatrixDim stressFromGlidePlane(const Eigen::Matrix<double,dim+1,1>& gpKey, const VectorDim& Rfield) const {
            /*! the stress field induced by boundary segments that belongs to a given glide plane
             */
            MatrixDim temp(MatrixDim::Zero());
            radialSegmentsContainerConstIterator gpIter = radialSegmentsContainer.find(gpKey);
            if( gpIter != radialSegmentsContainer.end() ) {
                for (unsigned int iv=0; iv<gpIter->second.size(); iv++) {
                    temp+=gpIter->second[iv]->stress(Rfield);
                }
            }
            return temp;
        }
        
        
        /* displacement *******************************************************/
        VectorDim displacement(const VectorDim& Rfield, const VectorDim& S) const {
            /*! the displacement field induced by boundary segments stored in this container
             */
            VectorDim temp(VectorDim::Zero());
            for(typename BaseContainerType::const_iterator sIter=this->begin();sIter!=this->end();++sIter){
                temp+=sIter->displacement(Rfield,S);
            }
            return temp;
        }
        
        
        /* outputVirtualDislocations ******************************************/
        void outputVirtualDislocations (const int& outputFrequency, const unsigned int& runID) const {
            /*! Function to output the Virtual Dislocations structure
             */
            model::SequentialOutputFile<'B',true>::set_increment(outputFrequency);
            model::SequentialOutputFile<'B',true>::set_count(runID);
            model::SequentialOutputFile<'B',true> BSFile;
            for(typename BaseContainerType::const_iterator sIter=this->begin();sIter!=this->end();++sIter){
                BSFile << std::setprecision(15)<<std::scientific << sIter->sourceP.transpose() << "  " <<
                sIter->sinkP.transpose()   << "  " <<
                sIter->Burgers.transpose() << "  " << std::endl;   // coreL shouldn't be outputted
            }
        }
        
    };
    /**************************************************************************/
    /**************************************************************************/
} // namespace model
#endif



//=======================================================================
// add a new entity to the container (from a restart file)
//======================================================================
//      template <typename SharedType>
//      void add (const double gp_H, const VectorDim gp_N, const VectorDim sourceP, const VectorDim sourceN, const VectorDim sinkP, const VectorDim sinkN,
//	        const VectorDim Burgers, const double coreL, const SharedType* sharedPtr) {

//	std::auto_ptr<VirtualBoundarySlipSurfaceType > pVBSS (new VirtualBoundarySlipSurfaceType (gp_H,gp_N,sourceP,sourceN,sinkP,sinkN,Burgers,coreL,sharedPtr) );
//	this->push_back(pVBSS);

//	//---------- check radial segments that are required for stress calculation
//	Eigen::Matrix<double,dim+1,1> gpKey;        gpKey << gp_N , gp_H;

//	radialSegmentsContainerIterator gpIter = radialSegmentsContainer.find(gpKey);

//	if (gpIter == radialSegmentsContainer.end()) {     // no segments from the same glide plane were stored before, so just push_back

//	  std::vector<VirtualBoundarySlipSurfaceType* > tempV;   tempV.push_back(&(*(this->rbegin())));
//	  radialSegmentsContainer.insert(std::make_pair (gpKey , tempV ) );

//	  this->rbegin()->addRadialSegments();
//	}

//	//-- else, check the segments repeatness --------
//	else {

//	  //--------- WARNING : removing the bool variables, and placing the functions call directly inside the if-statement
//	  //--------- may cause to have only the first function excuted if it returns TRUE, and so the call for the 2nd function will be neglected
//	  bool source_repeated = checkSegmentsRepeatness (gpIter, sourceP, Burgers, 0);
//	  bool sink_repeated   = checkSegmentsRepeatness (gpIter, sinkP  , Burgers, 1);
//	  if ( source_repeated || sink_repeated )  gpIter->second.push_back(&(*(this->rbegin())));

//	}

//	/*
//	for (gpIter = radialSegmentsContainer.begin(); gpIter!=radialSegmentsContainer.end(); gpIter++) {
//	  std::cout << "Container size " << gpIter->second.size() << std::endl;
//	}
//	std::cout << std::endl;

//	*/

//      }


/*
 //================================================================================================
 // function to read the input B_ files, that contain the virtual dislocations data
 //================================================================================================
 template <typename SharedType>
 void read(const unsigned int runID , const SharedType* sharedPtr) {
 
 VectorDim gp_N, sourceP, sourceN, sinkP, sinkN, Burgers;
 double gp_H, coreL;
 
 //std::cout<< " Node container size " << sharedPtr->nodeContainer.size() << std::endl;
 
 std::stringstream filename;
 filename << "B/B_" << runID << ".txt";
 
 FILE *fp =fopen(filename.str().c_str(), "r");
 
 //std::cout << "Reading Virtual dislocations file " <<runID << "   " << (fp == NULL) << std::endl;
 
 if (fp != NULL) {
 while (!feof(fp)) {
 if(fscanf(fp, "%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le",
 &gp_H, &gp_N(0) , &gp_N(1), &gp_N(2), &sourceP(0),&sourceP(1),&sourceP(2),&sourceN(0),&sourceN(1),&sourceN(2),
 &sinkP(0),&sinkP(1),&sinkP(2),&sinkN(0),&sinkN(1),&sinkN(2),
 &Burgers(0),&Burgers(1),&Burgers(2), &coreL)==20){
 
 add (gp_H, gp_N, sourceP, sourceN, sinkP, sinkN, Burgers, coreL ,sharedPtr );
 }
 }
 
 fclose(fp);
 }
 
 }
 
 */

/*
 //================================================================================================
 // function to read the input B_ files, that contain the virtual dislocations data
 //================================================================================================
 template <typename SharedType>
 void read(const unsigned int runID , const SharedType* sharedPtr) {
 
 VectorDim gp_N, sourceP, sourceN, sinkP, sinkN, Burgers;
 double gp_H, coreL;
 
 //std::cout<< " Node container size " << sharedPtr->nodeContainer.size() << std::endl;
 
 std::stringstream filename;
 filename << "B/B_" << runID << ".txt";
 
 FILE *fp =fopen(filename.str().c_str(), "r");
 
 //std::cout << "Reading Virtual dislocations file " <<runID << "   " << (fp == NULL) << std::endl;
 
 if (fp != NULL) {
 while (!feof(fp)) {
 if(fscanf(fp, "%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le%le",
 &gp_H, &gp_N(0) , &gp_N(1), &gp_N(2), &sourceP(0),&sourceP(1),&sourceP(2),&sourceN(0),&sourceN(1),&sourceN(2),
 &sinkP(0),&sinkP(1),&sinkP(2),&sinkN(0),&sinkN(1),&sinkN(2),
 &Burgers(0),&Burgers(1),&Burgers(2), &coreL)==20){
 
 add (gp_H, gp_N, sourceP, sourceN, sinkP, sinkN, Burgers, coreL ,sharedPtr );
 }
 }
 
 fclose(fp);
 }
 
 }
 
 */

