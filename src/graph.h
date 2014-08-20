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

#ifndef CGRAPH_H
#define CGRAPH_H

#include <assert.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <set>
#include <time.h>
#include <vector>

namespace triforce {

    static int CompareLong (const void * a, const void * b) {
        if ( *(long*)a > *(long*)b ) return 1;
        if ( *(long*)a < *(long*)b ) return -1;
        if ( *(long*)a == *(long*)b ) return 0;
    }

		/**	@brief This class represents a graph.*/
class Graph {
public:	
	Graph();
	~Graph();

	/** @brief Reads a graph from a file. The file must contain a list of
	*	undirected edges (an edge cannot appear twice). 
	*	The identifiers of the nodes have to bee between
	*	0 and N-1. The edges must be sorted by the first identifier first,
	*	and then the second.
	*	@param[in] fileName The name of the file.
	*	@return 0 if the load was successful. 1 if there were errors.*/
	int 		Load(const char * fileName, const int numThreads);

	/**	@brief Gets the number of nodes in the graph.
	*	@return The number of nodes.*/
	inline long 	GetNumNodes() const {
		return m_NumNodes;
	}

	/**	@brief Gets the number of edges in the graph.
	*	@return The number of edges.*/
	inline  long 	GetNumEdges() const {
			return m_NumEdges;
	}

	/**	@brief Gets the degree of a node.
	*	@param[in] nodeId The identifier of the node.
	*	@return	The degree of the node.*/
	inline  long 	GetDegree( long nodeId) const {
		assert(nodeId<m_NumNodes);
		return m_Nodes[nodeId+1] - m_Nodes[nodeId];
	}

	/**	@brief Gets the neighbors of a node.
	*	@param[in] nodeId The identifier of the node.
	*	@return	The neighbors of the node.*/
	inline const  long*  GetNeighbors( long nodeId) const {
		assert(nodeId<m_NumNodes);
		return &m_Adjacencies[m_Nodes[nodeId]];
	}

	/** @brief Gets the identifier corresponding to the given node.
	*	@param[in] nodeId The identifier of the node.
	*	@return The original identifier.*/
	inline  long 	ReMap( long nodeId) const {
		assert(nodeId<m_NumNodes);
		return m_Map[nodeId];
	}

    inline long GetEdgeIndex( long nodeId, long neighborOffset ) const {
       return m_Nodes[nodeId]+neighborOffset; 
    }


    inline long HasNeighbor( long nodeId, long neighbor ) const {
        const long* adjacencies = GetNeighbors(nodeId);
        long* pos = (long*) std::bsearch(&neighbor, adjacencies,
                GetDegree(nodeId), sizeof(long), CompareLong);
        if(pos) return static_cast<long>(pos - adjacencies);
        return -1;
    }

	/** @brief Returns the map between internal identifiers to external ones. Used by some tools.
	 *  @return The map vector.*/
	inline const  long* GetMap() const {
		return m_Map;
	}

private:
	long	    m_NumNodes; 		/**< @brief The number of nodes in the graph.*/
	long 	    m_NumEdges;			/**< @brief The number of edges in the graph.*/
	long*		m_Nodes;			/**< @brief The array of nodes of the graph.*/
	long* 		m_Adjacencies;		/**< @brief The array of adjacencies of the graph.*/
	long* 		m_Map;				/**< @brief The map of internal identifiers to original identifiers.*/
};

}
#endif
