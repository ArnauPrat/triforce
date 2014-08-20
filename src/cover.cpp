/*Triforce is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Triforce is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  */


#include "cover.h"
#include <algorithm>
#include <cassert>
#include <tuple>
#include <cstring>

namespace triforce {


    Cover* Create( const Graph* graph ) {
        Cover* cover = new Cover();
        cover->m_Graph = const_cast<Graph*>(graph);
        cover->m_MembershipStats = new MembershipStats[graph->GetNumNodes()];
        cover->m_NodeMemberships = new std::set<long>[graph->GetNumNodes()];
        cover->m_Weights = new long[graph->GetNumEdges()*2];
        cover->m_Communities = new std::set<long>[graph->GetNumNodes()];
        return cover;
    }

    void Copy( Cover& dest, const Cover& src ) {
        assert(dest.m_Graph == src.m_Graph);
        memcpy(dest.m_MembershipStats,src.m_MembershipStats, sizeof(MembershipStats)*src.m_Graph->GetNumNodes());
        memcpy(dest.m_Weights,src.m_Weights, sizeof(long)*src.m_Graph->GetNumEdges()*2);
        for( int i = 0; i < src.m_Graph->GetNumNodes(); ++i) {
            dest.m_NodeMemberships[i] = src.m_NodeMemberships[i];
            dest.m_Communities[i] = src.m_Communities[i];
        }
    }

    void   Destroy( Cover* cover ) {
        delete [] cover->m_MembershipStats;
        delete [] cover->m_NodeMemberships;
        delete [] cover->m_Weights;
        delete [] cover->m_Communities;
        delete cover;
    }

}

