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

#include "sparse_matrix.h"
#include "graph.h"
#include <set>
#include <vector>

namespace triforce {
   
   struct Movement; 

    class Cover {
        public: 

            class Community {
                public:
                    ~Community();

                    inline void Print() {
                        std::cout << m_Id << ": ";
                        for( long l : m_Nodes ) {
                            std::cout << l << " ";
                        }
                        std::cout << std::endl;
                    }

                    /** @brief Retrieves the size of the community.
                     *  @return The size of the community.*/
                    inline long Size() const {
                       return m_Nodes.size(); 
                    }

                    /** @brief Adds a node into the community.
                     *  @param nodeId The node to add.*/
                    inline void Add( long nodeId ) {
                        assert(!Contains(nodeId));
                        m_Nodes.insert( nodeId );
                        m_Cover.m_NodeCommunities[nodeId]->insert(m_Id);
                        const Graph& graph = GetGraph();
                        long degree = graph.GetDegree( nodeId );
                        const long* adjacencies = graph.GetNeighbors( nodeId );
                        for( long i = 0; i < degree; ++i ) {
                            if( Contains(adjacencies[i]) ) {
                                m_Kin[adjacencies[i]]++;
                                m_InternalEdges++;
                                m_ExternalEdges--;
                            } else {
                                m_ExternalEdges++;
                            }
                        }
                        m_Kin[nodeId] = m_InternalEdges;
                        /*const std::set<long>& communities = m_Cover.NodeCommunities(nodeId);
                        for(auto it = communities.begin(); it != communities.end(); ++it ) {
                            if( *it != m_Id ) {
                                m_Cover.m_CommunityMatrix.Inc( m_Id, *it );
                            }
                        }
                        */
                        
                    }

                    /** @brief Tests if a node belongs to the community.
                     *  @param nodeId The node to test.
                     *  @return true if the node belongs to the community.*/
                    inline bool Contains( long nodeId ) const {
                        return m_Nodes.find(nodeId) != m_Nodes.end();
                    }

                    /** @brief Removes a node from the community.
                     *  @param nodeId The node to remove.*/
                    inline void  Remove( long nodeId ) {
                        assert(Contains(nodeId));
                        m_Nodes.erase( nodeId );
                        m_Cover.m_NodeCommunities[nodeId]->erase(m_Id);
                        const Graph& graph = GetGraph();
                        long degree = graph.GetDegree( nodeId );
                        const long* adjacencies = graph.GetNeighbors( nodeId );
                        for( long i = 0; i < degree; ++i ) {
                            if( Contains(adjacencies[i]) ) {
                                m_InternalEdges--;
                                m_Kin[adjacencies[i]]--;
                                m_ExternalEdges++;
                            } else {
                                m_ExternalEdges--;
                            }
                        }
                        m_Kin.erase(nodeId);
                        /*const std::set<long>& communities = m_Cover.NodeCommunities(nodeId);
                        for(auto it = communities.begin(); it != communities.end(); ++it ) {
                            if( *it != m_Id ) {
                                m_Cover.m_CommunityMatrix.Dec( m_Id, *it );
                            }
                        }
                        */
                    }

                    /** @brief Gets the Id of the community.
                     *  @return The id of the community.*/
                    inline  long    Id() const {
                       return m_Id; 
                    }

                    /** @brief Gets the graph the community belong.
                     *  @return The graph.*/
                    inline const Graph& GetGraph() const {
                        return m_Cover.GetGraph();
                    }

                    /** @brief Gets the cover the community belongs to.
                     *  @return The cover.*/
                    inline const Cover& GetCover() const {
                        return m_Cover;
                    }

                    inline const std::set< long> Nodes() const {
                        return static_cast<const std::set< long> &>(m_Nodes);
                    }

                    /** @brief Gets the number of internal edges of the community.
                     *  @return The number of internal edges.*/
                    inline long GetInternalEdges() const {
                        return m_InternalEdges;
                    }

                    /** @brief Gets the number of internal edges of a node.
                     *  @param nodeId The node to retrieve the internal edges from.
                     *  @return The number of internal edges of the node.*/
                    inline long GetInternalEdges(long nodeId) const {
                        auto it = m_Kin.find(nodeId);
                        if( it != m_Kin.end() ) {
                            return (*it).second;
                        }
                        assert(0);
                        return -1;
                    }

                    /** @brief Gets the number of external edges of the community.
                     *  @return The number of external edges.*/
                    inline long GetExternalEdges() const {
                        return m_ExternalEdges;
                    }

                private:  
                    friend class Cover;
                    friend void triforce::Initialize( Cover& cover, double alpha, double overlapp );
                    friend void triforce::RefineCommunities( Cover& cover, double alpha , double overlapp );
                    Community( Cover & cover,  long id );
                    std::set<long> m_Nodes;
                    std::map<long, long> m_Kin;
                    Cover&         m_Cover;
                    const long     m_Id;
                    long           m_InternalEdges;
                    long           m_ExternalEdges;
            };

            Cover( const Graph& graph );
            ~Cover();

            /** @brief  Retrieve the community identifiers where a node belongs to.
             *  @params nodeId The node to retrieve the communities from.
             *  @return a set with the community identifiers.*/  
            const std::set<long> NodeCommunities( long nodeId ) const;


            /** @brief  Returns the number of communities a node belongs to.
             *  @params nodeId The node to retrieve the communities from.
             *  @return the number of communities */  
            long NumNodeCommunities( long nodeId ) const;


            /** @brief  Retrieves a community.
             *  @params communityId The community identifier .
             *  @return the community.*/  
            Community&  GetCommunity( long communityId );

            /** @brief  Retrieves a community.
             *  @params communityId The community identifier .
             *  @return the community.*/  
            const Community& GetCommunity( long communityId ) const;

            /** @brief Retrieves the graph of the cover.
             *  @return The graph of the cover.*/
            inline const Graph& GetGraph() const {
                return m_Graph;
            }

            /** @brief Retrieves the number of communities in the cover.
             *  @return The number of communities.*/
            inline long NumCommunities() const {
                return m_Communities.size();
            }


            /** @brief Retrieves the number of communities a node belongs to.
             *  @param nodeId The node.
             *  @return The number of communities a node belongs to.*/
            inline long NumCommunities( long nodeId ) const {
                return m_NodeCommunities[nodeId]->size();
            }

            /** @brief Clears the cover.*/
            void Clear();

            /** @brief Serializes the cover into an array.
             *  @param[out] The size of the created array.
             *  @return An array containign community sizes and communities.*/
            long* Serialize( long& size );

            /** @brief Deserializes an array into a cover.
             *  @param array The serialized cover. 
             *  @param size The size of the array.*/
            void Deserialize( long* array, long size );


        private:

            friend class Community;
            friend bool triforce::Similar( const Cover& cover, const Cover::Community& communityA, const Cover::Community& communityB, double overlapp );
            friend bool triforce::ViolatesOverlappInsert( const long i, const Cover::Community& community);
            friend bool triforce::ViolatesOverlappRemove( const long i, const Cover::Community& community);
            friend void triforce::Initialize( Cover& cover, double alpha, double overlapp );
            friend void triforce::RefineCommunities( Cover& cover, double alpha , double overlapp );
            friend bool PerformBestMovement( Cover& cover, const long nodeId, double alpha, double overlapp, double bestScore, Movement& movement );

            /** @brief Initializes the cover with an initial assignment of nodes to communities.*/
//            void  Initialize();

            const Graph&                                        m_Graph;            /**< @brief The graph of the cover.*/
            std::vector<Community*>                             m_Communities;      /**< @brief The vector of communities in the graph.*/
            std::vector<std::set<long>*>                        m_NodeCommunities;  /**< @brief The communities each node belong.*/
            SparseMatrix<long>                                  m_CommunityMatrix;  /**< @brief The matrix of community overlapps.*/
    };
}
#endif 
