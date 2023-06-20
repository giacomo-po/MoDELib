/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNetworkRemesh_h_
#define model_DislocationNetworkRemesh_h_

#include <DislocationDynamicsModule.h>

namespace model
{
    
    /*! \brief Class template that handles the nodal remesh of the DislocationNetwork.
     */
    template <typename DislocationNetworkType>
    class DislocationNetworkRemesh
    {
        
        typedef TypeTraits<DislocationNetworkType> TraitsType;
        static constexpr int dim=TraitsType::dim;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef typename TraitsType::VectorLowerDim VectorLowerDim;

        
        DislocationNetworkType& DN;
        
    public:
        
        const double Lmax;
        const double Lmin;
        const double relativeAreaThreshold;
        const short unsigned int remeshFrequency;
        
        DislocationNetworkRemesh(DislocationNetworkType& DN_in);
        void remesh(const long int &runID);
        void remeshByRemoval();
        void remeshByExpansion();
        void contractBoundaryNodes(); //This function contracts boundary nodes if they are at the same position and share atleast one same neightbors
        void contract0chordSegments();
        void remove0AreaLoopAcrossBnd();
//        static double minMeshSize(const SimplicialMesh<dim> &mesh);
    };
    
} // namespace model
#endif

