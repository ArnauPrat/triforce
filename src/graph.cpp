/*SCD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SCD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <algorithm>
#include <assert.h>
#include <list>
#include <string.h>

#ifndef _WINDOWS
#include <sys/time.h>
#endif 

#include <cstdlib>
#include <omp.h>
#include <map>
#include "graph.h"

namespace triforce {

    /**	@brief Compares two  integers.
     * 	@param e1 Void pointer to the first  integer.
     * 	@param e2 Void pointer to the second  integer.
     *	@return -1 if e1 goes before e2. 1 if e1 goes after e2. 0 if e1 and e2 are equal.*/
    static int 	Compare_Ids(const void* e1, const void* e2) {
        long id1 = *(long*)e1;
        long id2 = *(long*)e2;
        if (id1 < id2) return -1;
        if (id2 < id1) return 1;
        return 0;
    }

    Graph::Graph() :
        m_NumNodes(0),
        m_NumEdges(0),
        m_Nodes(NULL),
        m_Adjacencies(NULL),
        m_Map(NULL)
    {

    }

    Graph::~Graph() {

        if (m_Nodes != NULL) {
            delete[] m_Nodes;
            m_Nodes = NULL;
        }

        if (m_Adjacencies != NULL) {
            delete[] m_Adjacencies;
            m_Adjacencies = NULL;
        }

        if (m_Map != NULL) {
            delete[]  m_Map;
            m_Map = NULL;
        }
    }

    int Graph::Load(const char * fileName, const int numThreads) {

        printf("Graph: Loading Graph\n");
        std::ifstream inFile;
        inFile.open((const char *)fileName);
        if (!inFile) {
            printf("Graph: Error Openning Graph File\n");
            return 1;
        }


        printf("Graph: Relabeling nodes ...\n");
#ifndef _WINDOWS
        timeval time;
        gettimeofday(&time, NULL);
        long initTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif
        std::map< long, long>* mapa = new std::map< long, long>();
        if (!mapa) {
            printf("\t Graph: Error allocating mapa\n");
            return 1;
        }
        long index = 0;
        m_NumEdges = 0;
        long node1;
        while (inFile >> node1) {
            long node2;
            inFile >> node2;

            if (!mapa->count(node1)) {
                mapa->insert(std::pair< long, long>(node1, index));
                index++;
            }

            if (!mapa->count(node2)) {
                mapa->insert(std::pair< long, long>(node2, index));
                index++;
            }
            m_NumEdges++;
        }
#ifndef _WINDOWS
        gettimeofday(&time, NULL);
        long endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
        printf("Graph: Nodes relabeled in %lu ms\n", endTime - initTime);
#endif

        printf("Graph: Reading degrees ...\n");
        //We set the file cursor to the beginning.
        inFile.close();
        inFile.open((const char *)fileName);
#ifndef _WINDOWS
        gettimeofday(&time, NULL);
        initTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif
        //Allocate space for nodes and initialize de degree.
        m_NumNodes = index;
        m_Nodes = new  long[m_NumNodes + 1];
        if (!m_Nodes) {
            printf("Graph: Error allocating nodes\n");
            return 1;
        }

        for (long i = 0; i <= m_NumNodes; i++) {
            m_Nodes[i] = 0;
        }

        //Compute the degree of each node.
        while (inFile >> node1) {
            long node2;
            inFile >> node2;
            m_Nodes[(*mapa->find(node1)).second]++;
            m_Nodes[(*mapa->find(node2)).second]++;
            m_Nodes[m_NumNodes] += 2;
        }

        //Computing the adjacency indices, average degree and maximum Degree.
        double averageDegree = 0.0f;
        double maxDegree = 0.0f;
        long currentAdjacencyIndex = 0;
        for (long i = static_cast<long>(m_NumNodes - 1); i >= 0; --i) {
            averageDegree += m_Nodes[i];
            if (m_Nodes[i] > maxDegree) {
                maxDegree = m_Nodes[i];
            }
            m_Nodes[i] = m_Nodes[i + 1] - m_Nodes[i];
        }
        averageDegree /= (m_NumNodes);
#ifndef _WINDOWS 
        gettimeofday(&time, NULL);
        endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
        printf("\t Graph: Degrees read in %lu ms\n", endTime - initTime);
#endif

        printf("Graph: Reading adjacencies ...\n");
#ifndef _WINDOWS
        gettimeofday(&time, NULL);
        initTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif
        //We set the file cursor to the beginning.
        inFile.close();
        inFile.open((const char *)fileName);

        m_Adjacencies = new long[m_NumEdges * 2];
        if (!m_Adjacencies) {
            printf("Graph: Error allocating adjacencies\n");
            return 1;
        }
        long* counters = new  long[m_NumNodes];
        for (long i = 0; i < m_NumNodes; i++) {
            counters[i] = 0;
        }

        //Filling adjacencies
        while (inFile >> node1) {
            long node2;
            inFile >> node2;
            long tail = (*mapa->find(node1)).second;
            long head = (*mapa->find(node2)).second;
            assert(counters[tail] < (m_Nodes[tail + 1] - m_Nodes[tail]));
            assert(counters[head] < (m_Nodes[head + 1] - m_Nodes[head]));
            m_Adjacencies[m_Nodes[tail] + counters[tail]] = head;
            m_Adjacencies[m_Nodes[head] + counters[head]] = tail;
            counters[tail]++;
            counters[head]++;
        }
        delete[] counters;
        inFile.close();
#ifndef _WINDOWS
        gettimeofday(&time, NULL);
        endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
        printf("\t Graph: Adjacencies read in %lu ms\n", endTime - initTime);

        printf( "Graph: Sorting adjacencies ...\n" );
        gettimeofday(&time, NULL);
        initTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif
        //Sorting adjacencies.
        //        #pragma omp parallel for schedule(static, 16)
        for (long i = 0; i < m_NumNodes; i++) {
            qsort(&m_Adjacencies[m_Nodes[i]], m_Nodes[i + 1] - m_Nodes[i], sizeof(long), Compare_Ids);
        }

#ifndef _WINDOWS
        gettimeofday(&time, NULL);
        endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
        printf("\t Graph: Adjacencies sorted in %lu ms\n", endTime - initTime);
#endif

        printf("Graph: Filling map array ...\n");
#ifndef _WINDOWS
        gettimeofday(&time, NULL);
        initTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
#endif
        //Filling the map array
        m_Map = new  long[m_NumNodes];
        for (std::map< long, long>::iterator it = mapa->begin(); it != mapa->end(); it++) {
            m_Map[(*it).second] = (*it).first;
        }

        delete mapa;
#ifndef _WINDOWS
        gettimeofday(&time, NULL);
        endTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
        printf("\t Graph: Map array filled in %lu ms\n", endTime - initTime);
#endif

        printf("Graph: Graph Loaded\n");
        printf("Graph: Number of Nodes: %lu\n", m_NumNodes);
        printf("Graph: Number of Edges: %lu\n", m_NumEdges);
        printf("Graph: Average Degree: %f\n", averageDegree);
        printf("Graph: Maximum Degree: %f\n", maxDegree);
        printf("Graph: Memory \n");
        printf("..............\n");
        long memNodes = static_cast<long>((char*)&m_Nodes[m_NumNodes] - (char*)&m_Nodes[0]);
        long memEdges = static_cast<long>((char*)&m_Adjacencies[m_NumEdges * 2] - (char*)&m_Adjacencies[0]);
        long memMap = static_cast<long>((char*)&m_Map[m_NumNodes - 1] - (char*)&m_Map[0]);
        printf("%-16s %-10lu Bytes\n", "Nodes:", memNodes);
        printf("%-16s %-10lu Bytes\n", "Adjacencies:", memEdges);
        printf("%-16s %-10lu Bytes\n", "Map:", memMap);
        printf("%-16s %-10lu Bytes\n", "Total:", memNodes + memEdges + memMap);
        printf("..............\n");
        return 0;
    }


    static int compareInt(const void * a, const void * b)
    {
        if (*(int*)a > *(int*)b) return 1;
        if (*(int*)a < *(int*)b) return -1;
        if (*(int*)a == *(int*)b) return 0;
    }

}
