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
        //        typedef std::vector<VectorDim> PlaneNormalContainerType;
        typedef           LatticeVector<dim>              LatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        
        typedef std::vector<ReciprocalLatticeDirectionType> PlaneNormalContainerType;
        
        
        //        typedef std::vector<const ReciprocalLatticeDirectionType*> PlaneNormalPointerContainerType;
        typedef std::vector<unsigned int> PlaneNormalIDContainerType;
        
        //        typedef std::vector<const VectorDim> PlaneNormalContainerType;
        
        
    private:
        static PlaneNormalContainerType planeNormalContainer;
        
        static Eigen::Matrix<double,dim,dim> C2G;
        
        
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
            
            // Get crystallographic plane normals from the CrystalStructure
            planeNormalContainer=CrystalStructure::template reciprocalPlaneNormals<dim>();
            
            // Rotate the LatticeBasis of the CrystalStructure using C2G
            LatticeBase<dim>::setLatticeBasis(C2G*CrystalStructure::template getLatticeBasis<dim>());
            
            model::cout<<magentaColor<<"Current Crystal Plane Normals are:"<<std::endl;
            for (unsigned int k=0; k<planeNormalContainer.size();++k)
            {
                std::cout<<"    "<<planeNormalContainer[k].cartesian().transpose()<<std::endl;
            }
            model::cout<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        static PlaneNormalIDContainerType find_slipSystem(const VectorDim& chord,
                                                          const LatticeVectorType& Burgers)
        {/*!
          */
            assert(  chord.norm()>tol && "CHORD HAS ZERO NORM");
            assert(Burgers.squaredNorm()>0 && "BURGERS HAS ZERO NORM");
            
            
            const VectorDim normalizedChord(chord.normalized());
            //            const VectorDim normalizedBurgers(Burgers.normalized());
            
            //
            PlaneNormalIDContainerType allowedSlipSystems;
            
            // Try to find a plane which has normal orthogonal to both the chord and the Burgers
            for (unsigned int k=0;k<planeNormalContainer.size();++k)
            {
                if(std::fabs(planeNormalContainer[k].cartesian().dot(normalizedChord))<tol
                   && planeNormalContainer[k].dot(Burgers)==0)
                {
                    //allowedSlipSystems.insert( *iter );
                    allowedSlipSystems.push_back(k);
                }
            }
            
            // If no planes are found, check only chord to detect possibly sessile segment
            if(allowedSlipSystems.size()==0)
            {
                for (unsigned int k=0;k<planeNormalContainer.size();++k)
                {
                    if(	std::fabs( planeNormalContainer[k].cartesian().dot(normalizedChord))<tol )
                    {
                        allowedSlipSystems.push_back(k);
                    }
                }
                if (allowedSlipSystems.size()<2)
                {
                    std::cout<<"chord="<<chord.transpose()<<std::endl;
                    std::cout<<"Burgers="<<Burgers.cartesian().transpose()<<std::endl;
                    for (const auto& planeNormal : planeNormalContainer)
                    {
                        
                        std::cout<<"n="<<planeNormal.cartesian().transpose()<<" |c*n|="<< std::fabs( planeNormal.cartesian().dot(normalizedChord)) << " tol is "<<tol<<std::endl;
                    }
                    assert(allowedSlipSystems.size()>=2 && "SESSILE SEGMENTS MUST FORM ON THE INTERSECTION OF TWO CRYSTALLOGRAPHIC PLANES.");
                }
            }
        
            return allowedSlipSystems;
        }
        
        /**********************************************************************/
        static const ReciprocalLatticeDirectionType& find_planeNormal(const VectorDim& chord,
                                                                      const LatticeVectorType& Burgers)
        {/*!@param[in] chord the chord of a DislocationSegment
          * @param[in] Burgers the Burgers vector of a DislocationSegment
          *\returns A const reference to the first vector in planeNormalContainer
          * which is orthogonal to both chord and Burgers.
          */
            const PlaneNormalIDContainerType allowedSlipSystems=find_slipSystem(chord,Burgers);
            return planeNormalContainer[allowedSlipSystems[0]]; // RETURNING THE FIRST PLANE FOUND IS SOMEWHAT ARBITRARY
        }
        
        /**********************************************************************/
        static const VectorDim find_sessileNormal(const VectorDim& chord,
                                                  const LatticeVectorType& Burgers)
        {
            
            const PlaneNormalIDContainerType allowedSlipSystems=find_slipSystem(chord,Burgers);
            
            
            VectorDim temp(chord.normalized().cross(Burgers.cartesian())); // CHANGE HERE
            double tempNorm(temp.norm());
            if (tempNorm<FLT_EPSILON) // a screw segment
            {
                //temp.normalize();
                assert(allowedSlipSystems.size()>=2);
                temp.setZero(); // allow glide on primary plane
            }
            else{ // not a screw segment
                if (allowedSlipSystems.size()>=2)
                { // a sessile segment
                    //temp=allowedSlipSystems.rbegin()->normal;
//                    temp=(*allowedSlipSystems.rbegin())->cartesian();
                    temp=planeNormalContainer[allowedSlipSystems[1]].cartesian();
                    
                }
                else
                { // a glissile segment
                    temp.setZero();
                }
            }
            
            assert(std::fabs(  chord.normalized().dot(temp))<tol && "CHORD AND NORMAL ARE NOT ORTHOGONAL");
            return temp;
        }
        
        
        /**********************************************************************/
        static std::vector<VectorDim> conjugatePlaneNormal(const LatticeVectorType& B,
                                                           const VectorDim& N)
        {
            
            int count(0);
            //			VectorDim temp(VectorDim::Zero());
            std::vector<VectorDim> temp;
            assert((std::fabs(B.cartesian().normalized().dot(N.normalized()))<tol) && "CANNOT DETERMINE CONJUGATE PLANE FOR SESSILE SEGMENT");
            for (const auto& planeNormal : planeNormalContainer)
            {
                //                std::cout<<*iter.transpose()<<std::endl;
                if(	 B.dot(planeNormal)==0 && N.normalized().cross(planeNormal.cartesian()).norm()>tol)
                {
                    //                    temp=*iter;
                    temp.push_back(planeNormal.cartesian());
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
        
        /**********************************************************************/
        static size_t planeID(const VectorDim& planeNormal)
        {
            size_t temp=planeNormalContainer.size();
            
            for (int n=0;n<planeNormalContainer.size();++n)
            {
                if ((planeNormalContainer[n].cartesian()-planeNormal).norm()<FLT_EPSILON)
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
        
        //        /**********************************************************************/
        //        static const std::vector<Eigen::Matrix<double,dim,1> >& planeNormals()
        //        {
        //            return planeNormalContainer;
        //        }
        
    };
    
    template <int dim>
    std::vector<ReciprocalLatticeDirection<dim>> CrystalOrientation<dim>::planeNormalContainer=FCC::reciprocalPlaneNormals<dim>();
    
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
