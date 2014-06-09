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
#include "graph.h"
#include "graph_utils.h"
#include <algorithm>
#include <tuple>

namespace triforce {

    Cover::Community::Community( Cover& cover,  long id ) :
        m_Cover(cover),
        m_Id(id)
    {
    }  

    Cover::Community::~Community() {
    }

    Cover::Cover( const Graph & graph ) : 
        m_Graph(graph),
        m_NodeCommunities(graph.GetNumNodes())
    {
        // Initialize Cover structures.
        for(  long i = 0; i < m_Graph.GetNumNodes(); ++i ) { 
            m_NodeCommunities[i] = new std::set<long>();
        }
        Initialize();

    }

    Cover::~Cover() {
        for(long i = 0; i < m_Graph.GetNumNodes(); ++i ) { 
            delete m_NodeCommunities[i];
        }
    }

    const std::set< long>& Cover::NodeCommunities(  long nodeId ) const {
        return static_cast<std::set< long>&>(*m_NodeCommunities[nodeId]);
    }

    Cover::Community& Cover::GetCommunity(  long communityId ) {
        return static_cast<Community&>(*m_Communities[communityId]);
    }

    const Cover::Community& Cover::GetCommunity(  long communityId ) const {
        return static_cast<Community&>(*m_Communities[communityId]);
    }

    /** @brief Compares two node clustering.
     *  @param a The first tuple to sort.
     *  @param b The second tuple to sort.   
     *  @return -1 if a goes before b. 1 if a goes after b. 0 if a and b are equal.*/
    static bool CompareClusterings(const std::tuple<  long, double,  long >& a, const std::tuple<  long, double,  long > b) {
        return std::tie( std::get<1>(a), std::get<2>(a) ) > std::tie( std::get<1>(b), std::get<2>(b) ) ;
    }

    void Cover::Initialize() {

        std::vector< std::tuple<  long, double,  long > > nC( m_Graph.GetNumNodes() ); 
        for(  long i = 0; i < m_Graph.GetNumNodes(); ++i ) { 
            nC.push_back( std::tuple<  long, double,  long >( i,
                        ClusteringCoefficient( m_Graph, i ),
                        m_Graph.GetDegree(i)) 
                    );  
        }   
        std::sort(nC.begin(), nC.end(), CompareClusterings );                     //Sorts the nC
        std::vector<bool> visited( m_Graph.GetNumNodes(), false );                //Vector to track those visited members.
        long nextLabel = 0;
        for( auto it = nC.begin(); it != nC.end(); ++it ) {    
            long nodeId = std::get<0>(*it);
            if(!visited[nodeId]) {
                visited[nodeId] = true;
                Community* community = new Community( *this, nextLabel );
                m_Communities.push_back(community);
                community->Add(nodeId);
                m_NodeCommunities[nodeId]->insert(community->Id());
                const  long* adjacencies = m_Graph.GetNeighbors(nodeId);
                long degree = m_Graph.GetDegree(nodeId);
                for( long j = 0; j < degree; ++j ) { 
                    visited[adjacencies[j]] = true;
                    m_NodeCommunities[adjacencies[j]]->insert(nextLabel);
                    community->Add(adjacencies[j]);
                }   
                nextLabel++;
            }   
        }   
    }
}

