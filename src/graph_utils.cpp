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

#include "graph_utils.h"

namespace triforce {

 long CommonNeighbors( const Graph& graph, 
        const  long nodeId1, 
        const  long nodeId2 ) {

         long i = 0;
         long j = 0;
        const  long* adjacencies1 = graph.GetNeighbors( nodeId1 );
        const  long* adjacencies2 = graph.GetNeighbors( nodeId2 );
         long degree1 = graph.GetDegree(nodeId1);
         long degree2 = graph.GetDegree(nodeId2);
         long common = 0;
        while( i < degree1 && j < degree2 ) { 
            if( adjacencies1[i] < adjacencies2[j] ) { 
                i++;
            }   
            else if( adjacencies1[i] > adjacencies2[j] ) { 
                j++;
            } else {
                common++;
                i++;
                j++;
            }   
        }   
        return common;
    }   


double ClusteringCoefficient( const Graph& graph, const long nodeId ) {

     long degree = graph.GetDegree( nodeId );
    const  long* adjacencies = graph.GetNeighbors( nodeId );
     long possible = degree * (degree - 1);
     long numTriangles = 0;
    for(  long i = 0; i < degree; ++i ) {
       numTriangles += CommonNeighbors( graph, nodeId, adjacencies[i]);
    } 
    return numTriangles/(double)possible;
}

}
