/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_InteractionBase_h
#define _model_InteractionBase_h



namespace model {
    
    
    
    //template <typename _ParticleType, typename UserSystemProperties = SystemProperties<> >
    template<typename Derived, typename _DataType, int _DataPerParticle, int _DataPerCellMoment> //char...otherChars
    struct InteractionBase
    {
        // note: "Derived" is used to keep resultVector distinct for different interactions
        
        //enum{DataType=MpiDataType<_DataType>::DataType};
        enum{DataPerParticle=_DataPerParticle};
        enum{DataPerCellMoment=_DataPerCellMoment};

        
        
        static std::vector<_DataType> resultVector;
        static std::vector<_DataType> momentVector;

//        /*****************************************/
//        template <typename DerivedParticle,int dim>
//        static ResultType getResult(const SpatialParticle<DerivedParticle,dim>& p){
//            const int rid(p.rID);
//            
//            return (ResultType()<<resultVector[rid*3+0],resultVector[rid*3+1],resultVector[rid*3+2]).finished();
//            
//            //resultVector[rid*3+1]
//            //resultVector[rid*3+2]
//        }
        
        static void resize(const unsigned int&  k, const _DataType& val = _DataType())
        {
            resultVector.resize(k,val);
        }
        
        
        static _DataType* data()
        {
            return resultVector.data();
        }
        
//        operator const char *()
//        {
//            m_pStr[m_nSize] = '\0';
//            return(m_pStr);
//        }
        
        
//        static _DataType& operator[](const unsigned int& k)
//        {
//            return  resultVector[k];
//        }
        
        //        /**********************************************************************/
        //        static ResultType get(const ChargedParticleType& cp1)
        //        {
        //            return ResultType((ResultType()<<this->resultVector[cp1.mpiID*3+0],
        //                               /*         */ this->resultVector[cp1.mpiID*3+1],
        //                               /*         */ this->resultVector[cp1.mpiID*3+2]).finished());
        //        }

        
    };
    
    
    // declare statica data
    template<typename Derived, typename _DataType, int _DataPerParticle, int _DataPerCellMoment> //char...otherChars
    std::vector<_DataType> InteractionBase<Derived,_DataType,_DataPerParticle,_DataPerCellMoment>::resultVector;

    template<typename Derived, typename _DataType, int _DataPerParticle, int _DataPerCellMoment> //char...otherChars
    std::vector<_DataType> InteractionBase<Derived,_DataType,_DataPerParticle,_DataPerCellMoment>::momentVector;

    
} // end namespace
#endif
