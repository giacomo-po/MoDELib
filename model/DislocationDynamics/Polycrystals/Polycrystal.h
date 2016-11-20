/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Polycrystal_H_
#define model_Polycrystal_H_

#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>
#include <utility>
#include <tuple>
#include <map>
#include <vector>
#include <deque>
#include <Eigen/Core>
#include <model/MPI/MPIcout.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundary.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <model/Utilities/EigenDataReader.h>
#include <model/DislocationDynamics/StressStraight.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundaryType.h>

namespace model
{
    
    
    
    template <int dim>
    class Polycrystal :
    /* base */ private std::map<size_t,Grain<dim>>,
    /* base */ private std::map<std::pair<size_t,size_t>,GrainBoundary<dim>>,
    /* base */ public std::vector<StressStraight<dim> >
    {
        
        typedef SimplicialMesh<dim> SimplicialMeshType;
        typedef MeshRegion<Simplex<dim,dim> > MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;

        static constexpr PeriodicElement<13,Isotropic> Al=PeriodicElement<13,Isotropic>();
        static constexpr PeriodicElement<28,Isotropic> Ni=PeriodicElement<28,Isotropic>();
        static constexpr PeriodicElement<29,Isotropic> Cu=PeriodicElement<29,Isotropic>();
        static constexpr PeriodicElement<74,Isotropic>  W=PeriodicElement<74,Isotropic>();

        unsigned int materialZ;
        
    public:
        
        const SimplicialMeshType& mesh;
        
        /**********************************************************************/
        Polycrystal(const SimplicialMeshType& mesh_in) :
        /* init */ materialZ(0),
        /* init */ mesh(mesh_in)
        {
            model::cout<<"Creating Polycrystal"<<std::endl;

            
        }
        
        /**********************************************************************/
        void init(const std::string& fullName)
        {
            
            for(const auto& rIter : MeshRegionObserverType::regions())
            {
                grains().emplace(rIter.second->regionID,*(rIter.second));
            }
            
            for(const auto& rgnBnd : mesh.regionBoundaries())
            {
                grainBoundaries().emplace(std::piecewise_construct,
                                          std::forward_as_tuple(rgnBnd.first),
                                          std::forward_as_tuple(rgnBnd.second,grain(rgnBnd.first.first),grain(rgnBnd.first.second)));
                //model::cout<<"mesh region "<<rIter.second->regionID<<" contains "<<rIter.second->size()<<" Simplex<"<<dim<<","<<dim<<">"<<std::endl;
            }
            
            EigenDataReader EDR;

            EDR.readScalarInFile(fullName,"material",materialZ); // material by atomic number Z
            //Material<Isotropic>::select(materialZ);
            
            for(auto& gr : grains())
            {
                gr.second.selectMaterial(materialZ);
                
                Eigen::Matrix<double,dim,dim> C2Gtemp(Eigen::Matrix<double,dim,dim>::Identity());
                EDR.readMatrixInFile(fullName,"C2G"+std::to_string(gr.second.grainID),C2Gtemp); // crystal-to-global orientation
                gr.second.rotate(C2Gtemp);
            }
            

            // Initialize GrainBoundary objects
            grainBoundaryDislocations().clear();
            model::cout<<yellowBoldColor<<"Initializing GrainBoundaries"<<defaultColor<<std::endl;
            for(auto& gb : grainBoundaries())
            {
                gb.second.initializeGrainBoundary(*this);
            }
            
            model::SequentialOutputFile<'L',1>::set_count(0); // Vertices_file;
            model::SequentialOutputFile<'L',1>::set_increment(1); // Vertices_file;
            model::SequentialOutputFile<'L',true> stressStraightFile;
            size_t n=0;
            for(const auto& sStraight : grainBoundaryDislocations())
            {
                stressStraightFile<<n<<"\t"<<sStraight.P0.transpose()<<"\t"<<sStraight.P1.transpose()<<"\t"<<sStraight.b.transpose()<<std::endl;
                n++;
            }
        }
        
        /**********************************************************************/
        const std::deque<GrainBoundaryType<dim>>& grainBoundaryTypes() const
        {
            switch (materialZ)
            {
                case Cu.Z:
                    return PeriodicElement<Cu.Z,Isotropic>::grainBoundaryTypes;
                    break;
                default:
                    assert(0 && "grainBoundaryTypes not implemented.");
                    break;
            }
        }

        
//        /**********************************************************************/
//        void createLatticePlanes()
//        {
//
//        }
        
        /**********************************************************************/
        Grain<dim>& grain(const size_t& k)
        {
            return grains().at(k);
        }
        
        /**********************************************************************/
        const Grain<dim>& grain(const size_t& k) const
        {
            return grains().at(k);
        }
        
        /**********************************************************************/
        const std::map<size_t,Grain<dim>>& grains() const
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<size_t,Grain<dim>>& grains()
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<std::pair<size_t,size_t>,GrainBoundary<dim>>& grainBoundaries()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::map<std::pair<size_t,size_t>,GrainBoundary<dim>>& grainBoundaries() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const GrainBoundary<dim>& grainBoundary(const size_t& i,
                                                const size_t& j) const
        {
            assert(i!=j && "GrainBoundary IDs cannot be the same.");
            return (i<j)? grainBoundaries().at(std::make_pair(i,j)) : grainBoundaries().at(std::make_pair(j,i));
        }
        
        /**********************************************************************/
        GrainBoundary<dim>& grainBoundary(const size_t& i,
                                                const size_t& j)
        {
            assert(i!=j && "GrainBoundary IDs cannot be the same.");
            return (i<j)? grainBoundaries().at(std::make_pair(i,j)) : grainBoundaries().at(std::make_pair(j,i));
        }
        
        /**********************************************************************/
        LatticeVectorType latticeVectorFromPosition(const VectorDimD& p,
                                                    const Simplex<dim,dim>* const guess) const
        {
            const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(p,guess));
            assert(temp.first && "Position not found in mesh");
            return grain(temp.second->region->regionID).latticeVector(p);
        }
        
        /**********************************************************************/
        LatticeVectorType latticeVectorFromPosition(const VectorDimD& p) const
        {
            return latticeVectorFromPosition(p,&(mesh.simplices().begin()->second));
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType reciprocalLatticeVectorFromPosition(const VectorDimD& p,
                                                    const Simplex<dim,dim>* const guess) const
        {
            const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(p,guess));
            assert(temp.first && "Position not found in mesh");
            return grain(temp.second->region->regionID).reciprocalLatticeVector(p);
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType reciprocalLatticeVectorFromPosition(const VectorDimD& p) const
        {
            return reciprocalLatticeVectorFromPosition(p,&(mesh.simplices().begin()->second));
        }
        
        /**********************************************************************/
        std::vector<StressStraight<dim>>& grainBoundaryDislocations()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::vector<StressStraight<dim>>& grainBoundaryDislocations() const
        {
            return *this;
        }
        
    };
    
    
} // end namespace
#endif

