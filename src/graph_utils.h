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


#ifndef GRAPH_UTILS_H
#define GRAPH_UTILS_H

#include "graph.h"

namespace triforce {

/** @brief Computes the number of common neighbros of two nodes.
 *  @param graph The graph the nodes belong to.
 *  @param nodeId1 The first node.
 *  @param nodeId2 The second node.
 *  @return The number of common neighbors.*/
 long CommonNeighbors( const Graph& graph, 
                        const  long nodeId1, 
                        const  long nodeId2 );

/** @brief Computes the clustering coefficient of a node.
 *  @param graph The graph the node belongs to.
 *  @param nodeId The node identifier.
 *  @return The clustering coefficient of the node.*/
double ClusteringCoefficient( const Graph& graph, const  long nodeId );

}

#endif
