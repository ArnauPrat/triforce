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

#include "cover_utils.h"
#include "graph.h"
#include <set>

namespace triforce {

    static bool ValidOverlapp( const Cover& cover,  long nodeId ) {
        const std::set< long>& communities = cover.NodeCommunities( nodeId );        
        for( auto it1 = communities.begin(); it1 != communities.end(); ++it1 ) {
            for( auto it2 = communities.begin(); it2 != communities.end(); ++it2 ) {
                if( *it1 != *it2 && Similar(cover.GetCommunity(*it1), cover.GetCommunity(*it2))) {
                    return false; 
                }
            }
        }
        return true;
    }

    bool Similar( const Cover::Community& communityA, const Cover::Community& communityB ) {
        const std::set< long>& nodesA = communityA.Nodes();
        const std::set< long>& nodesB = communityB.Nodes();
        auto it1 = nodesA.begin();
        auto it2 = nodesB.begin();
         long common = 0;
        while( (it1 != nodesA.end()) && (it2 != nodesB.end()) ) { 
            if( *it1 < *it2 ) { 
               ++it1; 
            }   
            else if( *it1 > *it2 ) { 
                ++it2;
            } else {
                common++;
                ++it1;
                ++it2;
            }   
        }   
        if( common / (double)(nodesA.size() + nodesB.size() - common) >= 0.7 ) {
            return true;
        }
        return false;
    }

    double Score( const Cover& cover,  long nodeId ) {
        const Graph& graph = cover.GetGraph();
        if( ValidOverlapp( cover, nodeId ) ) {
            std::set< long> unionSet = CommunityUnion( cover, nodeId );
            const  long* adjacencies = graph.GetNeighbors( nodeId );
             long degree = graph.GetDegree( nodeId );
             long kin = 0;
             long kout = 0;
            for(  long i = 0; i < degree; ++i ) {
                if( unionSet.find(adjacencies[i]) != unionSet.end() ) {
                    kin++;
                } else {
                    kout++;
                }
            }
            return kin /(double)( kout + unionSet.size() - 1);
        }
        return 0.0;
    }

    double Score( const Cover& cover ) {
        const Graph& graph = cover.GetGraph();
         long numNodes = graph.GetNumNodes();
        double sum = 0;
        for(  long i = 0; i < numNodes; ++i ) {
            sum += Score( cover, i );
        }
        return sum / (double)numNodes;
    }

    std::set< long> CommunityUnion( const Cover& cover, 
                                            const  long nodeId ) {
        std::set< long> unionSet;
        const std::set< long>& communities = cover.NodeCommunities( nodeId );        
        for( auto it1 = communities.begin(); it1 != communities.end(); ++it1 ) {
            const std::set< long>& communityNodes = cover.GetCommunity(*it1).Nodes();
            for( auto it2 = communityNodes.begin(); it2 != communityNodes.end(); ++it2) {
                unionSet.insert(*it2);
            }
        }
        return unionSet;
    }
}
