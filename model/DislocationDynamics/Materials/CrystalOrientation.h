/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CRYSTALORIENTATION_H_
#define model_CRYSTALORIENTATION_H_

#include <assert.h>
#include <float.h>
#include <vector>
#include <iterator>     // std::distance
#include <Eigen/Dense>
#include <model/Utilities/TerminalColors.h>
#include <model/DislocationDynamics/Materials/CrystalStructures.h>
#include <model/MPI/MPIcout.h>
#include <model/LatticeMath/LatticeBase.h>



namespace model {
    
    /**************************************************************************/
    /**************************************************************************/
    template <int dim>
    class CrystalOrientation
    {
        
    public:
        typedef Eigen::Matrix<double,dim,1> VectorDim;
//        typedef std::vector<VectorDim> PlaneNormalContainerType;
        typedef std::vector<VectorDim> PlaneNormalContainerType;
        typedef std::vector<typename PlaneNormalContainerType::const_iterator> PlaneNormalIteratorContainerType;
        //        typedef std::vector<const VectorDim> PlaneNormalContainerType;

        
    private:
        static PlaneNormalContainerType planeNormalContainer;
        
        static Eigen::Matrix<double,dim,dim> C2G;

//        static Eigen::Matrix<double,dim,dim> latticeMatrix;
//        static Eigen::Matrix<double,dim,dim> inverseLatticeMatrix;

        
    public:
        
        
        static double tol;
        
        /* rotate *************************************************************/
        template <typename CrystalStructure>
        static void rotate(const Eigen::Matrix<double,dim,dim>& C2G_in)
        {
            
			// make sure that C2G is orthogonal
			assert((C2G_in*C2G_in.transpose()-Eigen::Matrix<double,dim,dim>::Identity()).norm()<2.0*DBL_EPSILON*dim*dim && "CRYSTAL TO GLOBAL ROTATION MATRIX IS NOT ORTHOGONAL.");
            // make sure that C2G is proper
            assert(std::fabs(C2G_in.determinant()-1.0) < FLT_EPSILON && "C2G IS NOT PROPER.");
            // store C2G
            C2G=C2G_in;
            
            const PlaneNormalContainerType tempPN(CrystalStructure::template getPlaneNormals<dim>());
            planeNormalContainer.clear();
            
            for (unsigned int k=0; k<tempPN.size();++k)
            {
                VectorDim temp(C2G*tempPN[k]);
                double tempNorm(temp.norm());
                assert(tempNorm > FLT_EPSILON && "PLANE NORMAL HAS ZERO NORM.");
                planeNormalContainer.push_back(temp/tempNorm);
            }
            
            
//            latticeMatrix=C2G*CrystalStructure::template getLatticeMatrix<dim>();
//            inverseLatticeMatrix=latticeMatrix.inverse();
            LatticeBase<dim>::setLatticeBasis(C2G*CrystalStructure::template getLatticeBasis<dim>());

            
            
            model::cout<<magentaColor<<"Current Crystal Plane Normals are:"<<std::endl;
            for (unsigned int k=0; k<planeNormalContainer.size();++k)
            {
                std::cout<<"    "<<planeNormalContainer[k].transpose()<<std::endl;
            }
            model::cout<<defaultColor<<std::endl;
            
		}
        
        /**********************************************************************/
		static void find_slipSystem(const VectorDim& chord, const VectorDim& Burgers,
                                    PlaneNormalIteratorContainerType& allowedSlipSystems)
        {/*!
          */
			assert(  chord.norm()>tol && "CHORD HAS ZERO NORM");
			assert(Burgers.norm()>tol && "BURGERS HAS ZERO NORM");
			
			
			const VectorDim normalizedChord(chord.normalized());
            const VectorDim normalizedBurgers(Burgers.normalized());
            
			allowedSlipSystems.clear();
            
            // Try to find a plane which has normal orthogonal to both the chord and the Burgers
            //			for (typename std::set<SlipSystem<dim,Nslips> >::iterator iter=slipSystemContainer.begin();iter!=slipSystemContainer.end();++iter){
            for (typename PlaneNormalContainerType::const_iterator iter=planeNormalContainer.begin();iter!=planeNormalContainer.end();++iter)
            {
				if(std::fabs(iter->dot(  normalizedChord))<tol &&
                   std::fabs(iter->dot(normalizedBurgers))<tol)
                {
					//allowedSlipSystems.insert( *iter );
					allowedSlipSystems.push_back( iter );
				}
			}
            
            if(allowedSlipSystems.size()==0)
            {
                for (typename PlaneNormalContainerType::const_iterator iter=planeNormalContainer.begin();iter!=planeNormalContainer.end();++iter)
                {
                    //						std::cout<<"c*n="<< std::fabs( iter->normal.dot(normalizedChord)) << " tol is "<<tol<<std::endl;
                    if(	std::fabs( iter->dot(normalizedChord))<tol )
                    {
                        //allowedSlipSystems.insert( *iter );
                        allowedSlipSystems.push_back( iter );
                    }
                }
                if (allowedSlipSystems.size()<2)
                {
                    std::cout<<"chord="<<chord.transpose()<<std::endl;
                    std::cout<<"Burgers="<<Burgers.transpose()<<std::endl;
                    for (typename PlaneNormalContainerType::const_iterator iter=planeNormalContainer.begin();iter!=planeNormalContainer.end();++iter){
                        
                        std::cout<<"n="<<iter->transpose()<<" |c*n|="<< std::fabs( iter->dot(normalizedChord)) << " tol is "<<tol<<std::endl;
                    }
                    assert(allowedSlipSystems.size()>=2 && "SESSILE SEGMENTS MUST FORM ON THE INTERSECTION OF TWO CRYSTALLOGRAPHIC PLANES.");
                }
            }
			
//			const unsigned int N(allowedSlipSystems.size());
//			switch (N) {
//				case 0: // CHECK FOR SESSILE
//
//					break;
//                    
//				case 1: // OK
//					break;
//                    
//				default: // More than one slip plane found. This must be a screw segment
//					break;
//			}
            
		}
        
        /**********************************************************************/
		static const VectorDim& find_planeNormal(const VectorDim& chord,
                                                 const VectorDim& Burgers)
        {/*!@param[in] chord the chord of a DislocationSegment
          * @param[in] Burgers the Burgers vector of a DislocationSegment
          *\returns A const reference to the first vector in planeNormalContainer 
          * which is orthogonal to both chord and Burgers.
          */
            PlaneNormalIteratorContainerType allowedSlipSystems;
			//shared.material.find_slipSystem(chord,Burgers,allowedSlipSystems);
            //            CrystalBase<dim,Nslips>::find_slipSystem(chord,Burgers,allowedSlipSystems);
            find_slipSystem(chord,Burgers,allowedSlipSystems);
            assert(std::fabs(  chord.dot(**allowedSlipSystems.begin()))<tol && "CHORD AND NORMAL ARE NOT ORTHOGONAL");
			return **allowedSlipSystems.begin(); // RETURNING THE FIRST PLANE FOUND IS SOMEWHAT ARBITRARY
		}
        
        /**********************************************************************/
		static VectorDim get_sessileNormal(const VectorDim& chord, const VectorDim& Burgers)
        {
            //			assert(chord.norm()>FLT_EPSILON && "CHORD TOO SMALL");
            //			assert(Burgers.norm()>FLT_EPSILON && "Burgers TOO SMALL");
            //            std::set<SlipSystem<dim,Nslips> > allowedSlipSystems;
            PlaneNormalIteratorContainerType allowedSlipSystems;
            
            //shared.material.find_slipSystem(chord,Burgers,allowedSlipSystems);
            //            CrystalBase<dim,Nslips>::find_slipSystem(chord,Burgers,allowedSlipSystems);
            find_slipSystem(chord,Burgers,allowedSlipSystems);
            
            VectorDim temp(chord.normalized().cross(Burgers));
			double tempNorm(temp.norm());
			if (tempNorm<FLT_EPSILON) // a screw segment
            { 
				//temp.normalize();
                assert(allowedSlipSystems.size()>=2);
                temp.setZero(); // allow glide on primary plane
			}
			else{ // not a screw segment
                if (allowedSlipSystems.size()>=2){ // a sessile segment
                    //temp=allowedSlipSystems.rbegin()->normal;
                    temp=**allowedSlipSystems.rbegin();
                    
                }
                else{ // a glissile segment
                    temp.setZero();
                }
			}
            
            assert(std::fabs(  chord.normalized().dot(temp))<tol && "CHORD AND NORMAL ARE NOT ORTHOGONAL");
            //assert(std::fabs(Burgers.dot(temp))<FLT_EPSILON && "BURGERS AND NORMAL ARE NOT ORTHOGONAL");
            
            
			return temp;
		}
		
        
        /**********************************************************************/
        static PlaneNormalContainerType conjugatePlaneNormal(const VectorDim& B, const VectorDim& N)
        {
            
			int count(0);
//			VectorDim temp(VectorDim::Zero());
            PlaneNormalContainerType temp;
            assert((std::fabs(B.normalized().dot(N.normalized()))<tol) && "CANNOT DETERMINE CONJUGATE PLANE FOR SESSILE SEGMENT");
            for (typename PlaneNormalContainerType::const_iterator iter=planeNormalContainer.begin();iter!=planeNormalContainer.end();++iter){
//                std::cout<<*iter.transpose()<<std::endl;
                if(	 std::fabs(B.normalized().dot(*iter))<tol && N.normalized().cross(*iter).norm()>tol){
//                    temp=*iter;
                    temp.push_back(*iter);
                    ++count;
                }
			}
			//assert(count==1 && "FOUND MORE THAN ONE CONJUGATE PLANES"); // IN BCC THERE IS MORE THAN ONE PLANE!
			return temp;
		}
        
        /**********************************************************************/
        static const Eigen::Matrix<double,dim,dim>& c2g()
        {
            return C2G;
        }
        
//        static const Eigen::Matrix<double,dim,dim>& lattMat()
//        {
//            return latticeMatrix;
//        }
//        
//        static const Eigen::Matrix<double,dim,dim>& invLattMat()
//        {
//            return inverseLatticeMatrix;
//        }
        
        
        /**********************************************************************/
        static size_t planeID(const VectorDim& planeNormal)
        {
            size_t temp=planeNormalContainer.size();
            
            for (int n=0;n<planeNormalContainer.size();++n)
            {
            if ((planeNormalContainer[n]-planeNormal).norm()<FLT_EPSILON)
            {
                temp=n;
            }
            }
                        assert(temp!=planeNormalContainer.size() && "PLANE NORMAL NOT FOUND");
            return temp;
            
//            typename PlaneNormalContainerType::iterator pIter(std::find(planeNormalContainer.begin(),planeNormalContainer.end(),planeNormal));
//            assert(pIter!=planeNormalContainer.end() && "PLANE NORMAL NOT FOUND");
//            return std::distance(planeNormalContainer.begin(),pIter);
        }
        
        /**********************************************************************/
        static const std::vector<Eigen::Matrix<double,dim,1> >& planeNormals()
        {
            return planeNormalContainer;
        }
        
    };
    
    template <int dim>
    std::vector<Eigen::Matrix<double,dim,1> > CrystalOrientation<dim>::planeNormalContainer=FCC::getPlaneNormals<dim>();
    
    template <int dim>
    Eigen::Matrix<double,dim,dim> CrystalOrientation<dim>::C2G=Eigen::Matrix<double,dim,dim>::Identity();

//    template <int dim>
//    Eigen::Matrix<double,dim,dim> CrystalOrientation<dim>::latticeMatrix=FCC::getLatticeMatrix<dim>();
//
//    template <int dim>
//    Eigen::Matrix<double,dim,dim> CrystalOrientation<dim>::inverseLatticeMatrix=latticeMatrix.inverse();

    
    template <int dim>
    double CrystalOrientation<dim>::tol=FLT_EPSILON;

    
    /**************************************************************************/
} // namespace model
#endif
