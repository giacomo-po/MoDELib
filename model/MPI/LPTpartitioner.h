/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_LPTpartitioner_h_
#define _model_LPTpartitioner_h_

#include <map>
#include <vector>
#include <deque>
#include <algorithm>    // std::min_element, std::max_element
#include <assert.h>
#include <iostream>
#include <iomanip>

namespace model {
    
    
    /*!\brief Temaplte implementation of the Longest Processing Time (LPT) algorithm for
     * task partiotioning.
     */
    template <typename T>
    class LPTpartitioner :
    /* inheritance    */ private std::multimap<double,T*>,
    /* inheritance    */ private std::vector<std::deque<T*> >,
    /* inheritance    */ private std::vector<double>
    {
        
   
        
    public:

        typedef std::multimap<double,T*> BaseMapType;
        typedef std::deque<T*> BaseDequeType;
        typedef std::vector<BaseDequeType> BaseDequeVectorType;
        typedef std::vector<double> BaseWeightVectorType;
        
        
        void clear()
        {
            BaseMapType::clear();
            BaseDequeVectorType::clear();
            BaseWeightVectorType::clear();
        }
        
        /**********************************************************************/
        void insert(const double& w,T* pT)
        {
            assert(w>=0 && "WEIGHT MUST BE POSITIVE");
            BaseMapType::insert(std::make_pair(w,pT));
        }
        
        /**********************************************************************/
        void partition(const int& nProc)
        {
            // Resize
            BaseDequeVectorType::resize(nProc);
            BaseWeightVectorType::resize(nProc,0.0);
            
            // Traverse BaseMapType in descending order (highest to lowest weight)
            for (typename BaseMapType::const_reverse_iterator rIter =BaseMapType::rbegin();
                 /*                                        */ rIter!=BaseMapType::rend();
                 /*                                      */ ++rIter)
            {
                // Compute the index of minimum weight in BaseWeightVectorType
                const size_t argMin(min_element(BaseWeightVectorType::begin(), BaseWeightVectorType::end())-BaseWeightVectorType::begin());
                // Add the current weigth to the bin of minimum weigth
                BaseWeightVectorType::operator[](argMin)+=rIter->first;
                // Add the current pointer to the bin of minimum weight
                BaseDequeVectorType::operator[](argMin).push_back(rIter->second);
            }
            
        }
        
        /**********************************************************************/
        double totalWeight() const
        {
            double temp(0.0);
            for (unsigned int k=0;k<BaseWeightVectorType::size();++k)
            {
                temp+=weight(k);
            }
            return temp;
        }


        /**********************************************************************/
        size_t bins() const
        {
            return BaseDequeVectorType::size();
        }
        
        /**********************************************************************/
        const BaseDequeType& bin(const int& k) const
        {
            return BaseDequeVectorType::operator[](k);
        }

//        /**********************************************************************/
//        const size_t& binSize(const int& k) const
//        {
//            return BaseDequeVectorType::operator[](k).size();
//        }
        
//        /**********************************************************************/
//        size_t binOffset(const int& k) const
//        {
//            size_t n(0);
//            for (int j=0;j<k;++j)
//            {
//                n+=binSize(j);
//            }
//            
//            return n;
//        }
        
        /**********************************************************************/
        const double& weight(const int& k) const
        {
            return BaseWeightVectorType::operator[](k);
        }
        
        /**********************************************************************/
        void show() const
        {
            const double wT(totalWeight());
            for (unsigned int k=0;k<BaseWeightVectorType::size();++k)
            {
                double wr(weight(k)/wT);
                int nC(wr*20*BaseWeightVectorType::size());
                std::cout<<"processor "<<std::setw(3)<<k<<std::setw(5)<<" ("<<wr*100<<"%): ";

                for (int c=0;c<nC;++c)
                {
                    std::cout<<"*";
                }
                std::cout<<"\n";
            }
        }
        
        
    };
} // end namespace
#endif