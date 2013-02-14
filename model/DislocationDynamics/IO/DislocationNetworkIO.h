/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONNETWORKIO_H_
#define model_DISLOCATIONNETWORKIO_H_

#include <string>

#include <model/Utilities/SequentialOutputFile.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>

namespace model {
	
	std::string defaultColor    = "\033[0m";	  // the default color for the console
	std::string redBoldColor    = "\033[1;31m";   // a bold red color
	std::string greenBoldColor  = "\033[1;32m";   // a bold green color
	std::string blueBoldColor   = "\033[1;34m";   // a bold blue color
	std::string magentaColor    = "\033[0;35m";   // a magenta color
	
	/**************************************************************************/
	/**************************************************************************/
    template <typename DislocationNetworkType>
	class DislocationNetworkIO
    {
        
        const DislocationNetworkType& DN; // a const reference to the DislocationNetwork
        
        enum {dim=DislocationNetworkType::dim}; // make dim available outside class

        
    public:
        
        typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
        typedef typename DislocationNetworkType::NetworkNodeContainerType NetworkNodeContainerType;
        typedef typename DislocationNetworkType::GlidePlaneObserverType GlidePlaneObserverType;
        typedef typename DislocationNetworkType::SpaceCellObserverType SpaceCellObserverType;
        typedef typename SpaceCellObserverType::CellMapType CellMapType;
        
        static int outputFrequency;
        static bool outputGlidePlanes;
        static bool outputSpaceCells;
        static bool outputPKforce;
        static bool outputMeshDisplacement;
        
        /* Constructor ******************************/
		DislocationNetworkIO(const DislocationNetworkType& DN_in) :
        /* init list */ DN(DN_in) // initialize DN
        {}

        /* outputTXT **********************************************************/
		void outputTXT(const unsigned int& runID) const
        {/*! Outputs DislocationNetwork data to the following files (x is the runID):
          * ./E/E_x.txt (DislocationSegment(s) are always outputted)
          * ./V/V_x.txt (DislocationNode(s) are always outputted)
          * ./C/C_x.txt (DislocationCell(s) only if outputSpaceCells==true)
          * ./G/G_x.txt (GlidePlane(s) only if outputGlidePlanes==true)
          * ./P/P_x.txt (PK force only if outputPKforce==true)
          * ./D/D_x.txt (mesh displacement only if outputMeshDisplacement==true)
          */
			std::cout<<"		Writing to "<<std::flush;
			
			//! 1- Outputs the Edge informations to file E_*.txt where * is the current simulation step
			SequentialOutputFile<'E',1>::set_increment(outputFrequency); // edgeFile;
			SequentialOutputFile<'E',1>::set_count(runID); // edgeFile;
			SequentialOutputFile<'E',1> edgeFile;
			edgeFile << *(const NetworkLinkContainerType*)(&DN);
			std::cout<<" E/E_"<<edgeFile.sID<<std::flush;

			//! 2- Outputs the Vertex informations to file V_*.txt where * is the current simulation step
			SequentialOutputFile<'V',1>::set_increment(outputFrequency); // vertexFile;
			SequentialOutputFile<'V',1>::set_count(runID); // vertexFile;
			SequentialOutputFile<'V',1> vertexFile;
            vertexFile << *(const NetworkNodeContainerType*)(&DN);
			std::cout<<", V/V_"<<vertexFile.sID<<std::flush;
			
            if(outputSpaceCells){
                //! 3- Outputs the nearest neighbor Cell structures to file C_*.txt where * is the current simulation step
                SequentialOutputFile<'C',1>::set_increment(outputFrequency); // Cell_file;
                SequentialOutputFile<'C',1>::set_count(runID); // Cell_file;
                SequentialOutputFile<'C',1> Cell_file;
                //              SpaceCellObserverType SPC;
                int cID(0);
                for (typename CellMapType::const_iterator cellIter=SpaceCellObserverType::begin();cellIter!=SpaceCellObserverType::end();++cellIter){
                    Cell_file<<cID<<"\t"<<cellIter->second->cellID.transpose()<<"\t"<<cellSize<<std::endl;
                    ++cID;
                }
                std::cout<<", C/C_"<<Cell_file.sID<<std::flush;
            }
			
            if(outputGlidePlanes){
                //! 4- Outputs the glide planes
                SequentialOutputFile<'G',1>::set_increment(outputFrequency); // GlidePlanes_file;
                SequentialOutputFile<'G',1>::set_count(runID); // GlidePlanes_file;
                SequentialOutputFile<'G',1> glide_file;
                glide_file << *dynamic_cast<const GlidePlaneObserverType*>(&DN);
                std::cout<<", G/G_"<<glide_file.sID<<std::flush;
            }
            
            if(outputPKforce){
                SequentialOutputFile<'P',1>::set_increment(outputFrequency); // Edges_file;
                SequentialOutputFile<'P',1>::set_count(runID); // Edges_file;
                SequentialOutputFile<'P',1> p_file;
                // let's output the node velocity
                int ll=0;
                for (typename NetworkLinkContainerType::const_iterator linkIter=DN.linkBegin();linkIter!=DN.linkEnd();++linkIter){
                    const int qOrder(linkIter->second->rgauss.cols());
                    for (int q=0;q<qOrder;++q){
                        p_file << ll*qOrder+q<<" "<< linkIter->second->rgauss.col(q).transpose()<<" "<<linkIter->second->pkGauss.col(q).transpose()<<"\n";
                    }
                    ll++;
                }
                std::cout<<", P/P_"<<p_file.sID<<std::flush;
            }
            
            if (DN.shared.use_bvp){
                if(outputMeshDisplacement){
                    model::SequentialOutputFile<'D',1>::set_increment(outputFrequency); // Vertices_file;
                    model::SequentialOutputFile<'D',1>::set_count(runID); // Vertices_file;
                    model::SequentialOutputFile<'D',true> d_file;
                    for (unsigned int i = 0; i< DN.shared.domain.nodeContainer.size(); i++){
                        d_file<< DN.shared.domain.nodeContainer[i].sID<<"	" << (DN.shared.domain.nodeContainer[i].u+DN.shared.domain.nodeContainer[i].uInf).transpose()<<std::endl;
                    }
                    std::cout<<", D/D_"<<d_file.sID<<std::flush;
                }
                if(DN.shared.boundary_type==1){
                    DN.shared.vbsc.outputVirtualDislocations(outputFrequency,runID);
                    
                }
            }
            
			
#ifdef customUserOutputs
#include customUserOutputs
#endif
			
		}
        
        
        /* outputBIN **********************************************************/
        void outputBIN(const unsigned int& runID) const
        {
            //			typedef std::pair<std::pair<int,int>,Eigen::Matrix<double,1,9> > BinEdgeType;
            //			SequentialBinFile<'E',BinEdgeType>::set_increment(outputFrequency);
            //			SequentialBinFile<'E',BinEdgeType>::set_count(runID);
            //			SequentialBinFile<'E',BinEdgeType> binEdgeFile;
            //			for (typename NetworkLinkContainerType::const_iterator linkIter=this->linkBegin();linkIter!=this->linkEnd();++linkIter){
            //				Eigen::Matrix<double,1,9> temp;
            //				temp<< linkIter->second->flow.transpose(),
            //				/*  */ linkIter->second->glidePlaneNormal.transpose(),
            //				/*  */ linkIter->second->sourceTfactor,
            //				/*  */ linkIter->second->sinkTfactor,
            //				/*  */ linkIter->second->pSN()->sID;
            //				binEdgeFile.write(std::make_pair(linkIter->first,temp));
            //			}
        }
        
    };
    
    // Declare static data
    template <typename DislocationNetworkType>
    int DislocationNetworkIO<DislocationNetworkType>::outputFrequency=1;

    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputGlidePlanes=false;

    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputSpaceCells=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputPKforce=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputMeshDisplacement=false;
    
    /**************************************************************************/
    /**************************************************************************/
} // namespace model
#endif





//			for (typename NetworkNodeContainerType::const_iterator nodeIter=this->nodeBegin();nodeIter!=this->nodeEnd();++nodeIter){
//				vertexFile << nodeIter->second->sID<<"\t"
//				/*         */ << std::setprecision(15)<<std::scientific<<nodeIter->second->get_P().transpose()<<"\t"
//				/*         */ << std::setprecision(15)<<std::scientific<<nodeIter->second->get_T().transpose()<<"\t"
//				/*         */ << nodeIter->second->pSN()->sID<<"\t";
//				if (shared.use_bvp){ //output in deformed configuration
//					vertexFile << std::setprecision(15)<<std::scientific<<nodeIter->second->deformedPosition().transpose()<<"\t";
//				}
//				else{
//					vertexFile<< VectorDimD::Zero().transpose();
//				}
//				vertexFile << std::endl;
//
//			}

