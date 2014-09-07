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


#ifndef COVER_H
#define COVER_H

#include "graph.h"
#include <set>
#include <vector>

namespace triforce {
   struct Movement; 
   struct MembershipStats {
       double   m_Score;
       double   m_R;
       double   m_Din;
       double   m_DinPrima;
       double   m_ConnectedAdd;
       double   m_NotConnectedAdd;
       double   m_ConnectedRemove;
       double   m_NotConnectedRemove;
   };

   struct CommunityStats {
       int  m_InDegree;
       int  m_OutDegree;
       int  m_Din;
       int  m_DinPrima;
       int  m_R;
       double  m_Score;
   };


   struct Cover {
       Graph*                           m_Graph;
       MembershipStats*                 m_MembershipStats;
       std::set<long>*                  m_NodeMemberships;
       long*                            m_Weights;
       std::vector<std::set<long>*>     m_Communities;
       CommunityStats*                  m_CommunityStats;
   };

   Cover* Create( const Graph* graph );
   void Copy( Cover& dest, const Cover& src );
   void Destroy( Cover* cover );

}
#endif 
