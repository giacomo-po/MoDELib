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
#include <deque>
#include <iterator>     // std::distance
#include <Eigen/Dense>
#include <model/Utilities/TerminalColors.h>
//#include <model/DislocationDynamics/Materials/CrystalStructures.h>
#include <model/MPI/MPIcout.h>
#include <model/LatticeMath/LatticeBase.h>
#include <model/DislocationDynamics/Materials/FCCcrystal.h>
#include <model/DislocationDynamics/Materials/BCCcrystal.h>
#include <model/DislocationDynamics/Materials/SlipSystem.h>



namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <int dim>
    class CrystalOrientation
    {
        
    public:
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef           LatticeVector<dim>              LatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        //        typedef LatticePlaneBase ReciprocalLatticeDirectionType;
        typedef std::vector<LatticePlaneBase> PlaneNormalContainerType;
        typedef std::vector<SlipSystem> SlipSystemContainerType;
        typedef std::vector<unsigned int> PlaneNormalIDContainerType;
        
    private:
        static PlaneNormalContainerType planeNormalContainer;
        static SlipSystemContainerType slipSystemContainer;
        
        static Eigen::Matrix<double,dim,dim> C2G;
        
        
    public:
        
        static double tol;
        
        /* rotate *************************************************************/
        template <typename MaterialType>
        static void rotate(const Eigen::Matrix<double,dim,dim>& C2G_in)
        {
            typedef typename MaterialType::CrystalStructure CrystalStructure;
            
            // make sure that C2G is orthogonal
            assert((C2G_in*C2G_in.transpose()-Eigen::Matrix<double,dim,dim>::Identity()).norm()<2.0*DBL_EPSILON*dim*dim && "CRYSTAL TO GLOBAL ROTATION MATRIX IS NOT ORTHOGONAL.");
            // make sure that C2G is proper
            assert(std::fabs(C2G_in.determinant()-1.0) < FLT_EPSILON && "C2G IS NOT PROPER.");
            // store C2G
            C2G=C2G_in;
            
            // Get crystallographic plane normals from the CrystalStructure
            planeNormalContainer=CrystalStructure::template reciprocalPlaneNormals<dim>();
            slipSystemContainer=CrystalStructure::slipSystems();
            
            // Rotate the LatticeBasis of the CrystalStructure using C2G
            LatticeBase<dim>::setLatticeBasis(C2G*CrystalStructure::template getLatticeBasis<dim,MaterialType>());
            
            model::cout<<magentaColor<<"Current Crystal Plane Normals are:"<<std::endl;
            for (unsigned int k=0; k<planeNormalContainer.size();++k)
            {
                std::cout<<"    "<<planeNormalContainer[k].cartesian().normalized().transpose()<<std::endl;
            }
            model::cout<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        static PlaneNormalIDContainerType find_slipSystem(const LatticeVectorType& chord,
                                                          const LatticeVectorType& Burgers)
        {/*!
          */
            assert(  chord.squaredNorm()>0 && "CHORD HAS ZERO NORM");
            assert(Burgers.squaredNorm()>0 && "BURGERS HAS ZERO NORM");
            
            
            PlaneNormalIDContainerType allowedSlipSystems;
            
            // Try to find a plane which has normal orthogonal to both the chord and the Burgers
            for (unsigned int k=0;k<planeNormalContainer.size();++k)
            {
                if(planeNormalContainer[k].dot(chord)==0
                   && planeNormalContainer[k].dot(Burgers)==0)
                {
                    allowedSlipSystems.push_back(k);
                }
            }
            
            // If no planes are found, check only chord to detect possibly sessile segment
            if(allowedSlipSystems.size()==0)
            {
                for (unsigned int k=0;k<planeNormalContainer.size();++k)
                {
                    if(	planeNormalContainer[k].dot(chord)==0)
                    {
                        allowedSlipSystems.push_back(k);
                    }
                }
                if (allowedSlipSystems.size()<2)
                {
                    std::cout<<"chord="<<chord.cartesian().transpose()<<std::endl;
                    std::cout<<"Burgers="<<Burgers.cartesian().transpose()<<std::endl;
                    for (const auto& planeNormal : planeNormalContainer)
                    {
                        
                        std::cout<<"n="<<planeNormal.cartesian().normalized().transpose()<<" c*n="<< planeNormal.dot(chord)<<" b*n="<< planeNormal.dot(Burgers) <<std::endl;
                    }
                    assert(allowedSlipSystems.size()>=2 && "SESSILE SEGMENTS MUST FORM ON THE INTERSECTION OF TWO CRYSTALLOGRAPHIC PLANES.");
                }
            }
            
            return allowedSlipSystems;
        }
        
        /**********************************************************************/
        static const LatticePlaneBase& find_glidePlane(const LatticeVectorType& chord,
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
        static const LatticePlaneBase& find_sessilePlane(const LatticeVectorType& chord,
                                                         const LatticeVectorType& Burgers)
        {
            const PlaneNormalIDContainerType allowedSlipSystems=find_slipSystem(chord,Burgers);
            
            int planeID(0);
            if (allowedSlipSystems.size()>=2)
            {
                if (chord.cross(Burgers).squaredNorm()!=0) // a sessile segment
                {
                    planeID=1;
                }
                else // a screw segment
                {
                    planeID=0; // allow glide on primary plane
                }
            }
            return planeNormalContainer[allowedSlipSystems[planeID]];
        }
        
        
        /**********************************************************************/
        static std::deque<const LatticePlaneBase*> conjugatePlaneNormal(const LatticeVectorType& B,
                                                                        const ReciprocalLatticeDirectionType& N)
        {
            //            assert(B.dot(N)==0 && "CANNOT DETERMINE CONJUGATE PLANE FOR SESSILE SEGMENT");
            
            
            std::deque<const LatticePlaneBase*> temp;
            if(B.dot(N)==0) // not sessile
            {
                for (const auto& planeNormal : planeNormalContainer)
                {
                    if(	 B.dot(planeNormal)==0 && N.cross(planeNormal).squaredNorm()>0)
                    {
                        temp.push_back(&planeNormal);
                    }
                }
            }
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
            
            for (size_t n=0;n<planeNormalContainer.size();++n)
            {
                if ((planeNormalContainer[n].cartesian().normalized()-planeNormal.normalized()).norm()<FLT_EPSILON)
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
        static const PlaneNormalContainerType& planeNormals()
        {
            return planeNormalContainer;
        }
        
        /**********************************************************************/
        static const SlipSystemContainerType& slipSystems()
        {
            return slipSystemContainer;
        }
        
    };
    
    template <int dim>
    std::vector<LatticePlaneBase> CrystalOrientation<dim>::planeNormalContainer=FCC::reciprocalPlaneNormals<dim>();
    
    template <int dim>
    std::vector<SlipSystem> CrystalOrientation<dim>::slipSystemContainer=FCC::slipSystems();
    
    
    template <int dim>
    Eigen::Matrix<double,dim,dim> CrystalOrientation<dim>::C2G=Eigen::Matrix<double,dim,dim>::Identity();
    
    template <int dim>
    double CrystalOrientation<dim>::tol=FLT_EPSILON;
    
    
    /**************************************************************************/
} // namespace model
#endif
