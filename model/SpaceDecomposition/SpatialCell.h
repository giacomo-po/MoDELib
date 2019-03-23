/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SpatialCell_h_
#define model_SpatialCell_h_

#include <math.h>
#include <set>
//#include <boost/utility.hpp>
#include <utility>      // std::pair, std::make_pair
#include <memory> // std::shared_ptr (c++11)
#include <map>
#include <Eigen/Dense>
#include <model/Math/CompileTimeMath/CTM.h>
#include <model/Utilities/CompareVectorsByComponent.h>
//#include <model/Utilities/modelMacros.h> // model_execAssert(
#include <model/Utilities/CRTP.h>
#include <model/SpaceDecomposition/CellShift.h>
#include <model/Utilities/TypeTraits.h>
#include <model/Utilities/NonCopyable.h>
#include <model/Math/CompileTimeMath/PermutationWithRepetition.h>
//#include <stdexcept>      // std::out_of_range


namespace model
{
	
    //class pre-declaration
    template<typename ParticleType, short unsigned int dim>
	struct SpatialCellObserver;
    
	/**************************************************************************/
	/**************************************************************************/
	/*!\brief A dim-dimensional cell occupying the spatial region
     * cellID <= x/cellSize < (cellID+1). SpatialCell is aware off all
     * ParticleType objects present inside it.
	 */
	template<typename ParticleType, short unsigned int dim>
	struct SpatialCell : public NonCopyable,
//    /* inheritance    */ boost::noncopyable,
    /* inheritance    */ public TypeTraits<ParticleType>::SpatialCellProperties
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef Eigen::Matrix<int,dim,1> CellIdType;
        typedef SpatialCell<ParticleType,dim> SpatialCellType;
		typedef std::map<CellIdType, SpatialCellType* const,CompareVectorsByComponent<int,dim> >  CellMapType;
		typedef std::shared_ptr<SpatialCellType> SharedPtrType;
        typedef std::pair<bool,SpatialCellType* const> isCellType;
		typedef SpatialCellObserver<ParticleType,dim> SpatialCellObserverType;
        typedef std::set<const ParticleType*> ParticleContainerType; // PTR COMPARE IS NOT NECESSARY
        
        typedef typename TypeTraits<ParticleType>::SpatialCellProperties SpatialCellProperties;
        
        typedef Eigen::Matrix<double,dim,CTM::pow(2,dim)> MatrixDimVerticesType;

        typedef Eigen::Matrix<double,1,CTM::pow(2,dim)> VectorVerticesType;
        
        enum{neighborLayer=1}; // = (1+2*1)^dim  cells = 27  cells in 3d
//        enum{    nearLayer=2}; // = (1+2*3)^dim  cells = 343 cells in 3d
        
        typedef CellShift<dim,1>    CellShiftType;
//        typedef CellShift<dim,1>     NearShiftType;
        
        
        CellMapType  neighborCells;
        //        CellMapType      nearCells;
        //        CellMapType       farCells;
        
#ifdef _MODEL_MPI_
        int assignedRank;
#else
        const int assignedRank;
#endif
        
        /* isNeighborCell *****************************************************/
        bool isNeighborCell(const CellIdType& otherCellID) const
        {
            bool temp(false);
            for (int k=0;k<neighborCellIDs.cols();++k)
            {
                if(neighborCellIDs.col(k)==otherCellID)
                { // cell is a neighbor
                    temp=true;
                    break;
                }
            }
            return temp;
        }
        
        
	public:
        
		//! The container of pointers to particles in this cell
		ParticleContainerType particleContainer; // TO DO: MAKE THIS A BASE CLASS
        
        //! The ID of this cell, defining the dim-dimensional spatial region cellID<= x/cellSize < (cellID+1).
		const CellIdType cellID;
        
        //! The position vector of the center of this SpatialCell
        const VectorDimD center;
        
        static const MatrixDimVerticesType standardVertices;
        
		//! The cellID(s) of the neighboring SpatialCell(s) (in column)
		const Eigen::Matrix<int,dim, CellShiftType::Nneighbors> neighborCellIDs;
        
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        
		/**********************************************************************/
		SpatialCell(const CellIdType& cellID_in) :
        /* base init */ SpatialCellProperties(TypeTraits<ParticleType>::init()),
        /* init list */ assignedRank(0),
        /* init list */ cellID(cellID_in),
//        /* init list */ center((cellID.template cast<double>().array()+0.5).matrix()*SpatialCellObserverType::cellSize()),
        /* init list */ center(cellID.template cast<double>()*SpatialCellObserverType::cellSize()),
//        /* init list */ vertices(PermutationWithRepetition<dim>::permute((Eigen::Matrix<double,1,2>()<<0.5,-0.5).finished()*SpatialCellObserverType::cellSize()).colwise()+center),
//        /* init list */ vertices((PermutationWithRepetition<dim>::permute((Eigen::Matrix<double,1,2>()<<-0.5,0.5).finished())*SpatialCellObserverType::cellSize()).colwise()+center),
        /* init list */ neighborCellIDs(CellShiftType::neighborIDs(cellID))
        {/*! @param[in] cellID_in The ID of the current cell.
          * The constructor initializes cellID, center, neighborCellIDs, and
          * nearCellIDs. It then populates neighborCells, nearCells, and
          * farCells.
          */
            
			//! 1- Adds this to static SpatialCellObserver::cellMap
            const bool success=SpatialCellObserverType::cellMap.emplace(cellID,this).second;
            assert(success && "CANNOT INSERT Spatial CELL IN STATIC cellMap.");

            //! 2- Populate neighborCells
            for (int c=0;c<neighborCellIDs.cols();++c)
            {
                isCellType isC(SpatialCellObserverType::isCell(neighborCellIDs.col(c)));
                if (isC.first)
                {
                    const bool success1=neighborCells.emplace(isC.second->cellID,isC.second).second;
                    assert(success1 && "CANNOT INSERT CELL IN NEIGHBORCELLS");

                    if(cellID!=isC.second->cellID)
                    {
                        const bool success2=isC.second->neighborCells.emplace(cellID,this).second;
                        assert(success2 && "CANNOT INSERT THIS IN NEIGHBORCELLS");
                    }
                }
            }
            
		}
		
		/**********************************************************************/
		~SpatialCell()
        {/*! Removes this from static SpatialCellObserver::cellMap
          */
            for (typename CellMapType::iterator cIter=neighborCells.begin();cIter!=neighborCells.end();++cIter)
            {
                if(cellID!=cIter->second->cellID)
                {
                    const size_t erased=cIter->second->neighborCells.erase(cellID);
                    assert(erased==1 && "CANNOT ERASE SPATIALCELL FROM NEIGHBOR-CELLMAP.");
                }
            }
            const size_t erased1=SpatialCellObserverType::cellMap.erase(cellID);
            assert(erased1==1 && "CANNOT ERASE SPATIALCELL FROM CELLMAP.");
		}
		
		/**********************************************************************/
		void addParticle(const ParticleType* const pP)
        {/*!@param[in] pP pointer to a ParticleType
          * Adds pP to *this
          */
            const bool success=particleContainer.insert(pP).second;
			assert(success && "CANNOT INSERT PARTICLE IN SpatialCELL");
		}
		
		/**********************************************************************/
		void removeParticle(const ParticleType* const pP)
        {/*!@param[in] pP pointer to a ParticleType
          * Removes pP to *this
          */
            const size_t erased=particleContainer.erase(pP);
			assert(erased==1 && "CANNOT ERASE PARTICLE FROM particleContainer.");
		}
        
		/**********************************************************************/
        typename CellMapType::const_iterator neighborCellsBegin() const
        {/*!\returns A const iterator referring to the first neighbor cell 
          * of *this.
          */
            return neighborCells.begin();
        }
        
		/**********************************************************************/
        typename CellMapType::const_iterator neighborCellsEnd() const
        {/*!\returns A const iterator referring to the past-the-end neighbor 
          * cell of *this.
          */
            return neighborCells.end();
        }
        
		/**********************************************************************/
        typename ParticleContainerType::const_iterator particleBegin() const
        {
            return particleContainer.begin();
        }
        
		/**********************************************************************/
        typename ParticleContainerType::const_iterator particleEnd() const
        {
            return particleContainer.end();
        }
        
		/**********************************************************************/
        size_t size() const
        {
            return particleContainer.size();
        }
        
        /**********************************************************************/
        MatrixDimVerticesType vertices() const
        {
            return (standardVertices*SpatialCellObserverType::cellSize()).colwise()+center;
        }
        
//        /**********************************************************************/
//        VectorVerticesType vertexWeigths(const VectorDimD& P) const
//        {
//            return VectorVerticesType::Constant(1.0)-(vertices().colwise()-P).colwise().prod().array().abs().matrix()/CTM::pow(SpatialCellObserverType::cellSize(),dim);
//        }
        
        /**********************************************************************/
        VectorVerticesType vertexWeigths(const VectorDimD& P) const
        {            
            return (standardVertices.colwise()+(P-center)/SpatialCellObserverType::cellSize()).colwise().prod().array().abs();
        }
        
		/**********************************************************************/
        size_t neighborSize() const
        {
            size_t temp(0);
            for (typename CellMapType::const_iterator cIter =neighborCellsBegin();
                 /*                                */ cIter!=neighborCellsEnd();
                 /*                                */ cIter++)
            {
                temp+=cIter->second->size();
            }
            return temp;
        }
        
		/**********************************************************************/
        size_t n2Weight() const
        {
            return size()*neighborSize();
        }
        
		/**********************************************************************/
        template <typename T=double>
        T n2WeightD() const
        {
            return static_cast<T>(size())*neighborSize();
        }
        
        //        /* size ***************************************************************/
        //        template <size_t N>
        //        ??? getProperty()
        //        {
        //            return FINISH HERE;
        //        }
        
	};
    
    template<typename ParticleType, short unsigned int dim>
    const Eigen::Matrix<double,dim,CTM::pow(2,dim)> SpatialCell<ParticleType,dim>::standardVertices=PermutationWithRepetition<dim>::permute((Eigen::Matrix<double,1,2>()<<-0.5,0.5).finished());
    
}	// close namespace model
#endif

