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
#include <set>
#include <vector>

namespace triforce {

    
    class Graph; //forward declaration.

    class Cover {
        public: 

            class Community {
                public:
                    ~Community();

                    /** @brief Retrieves the size of the community.
                     *  @return The size of the community.*/
                    inline long Size() const {
                       return m_Nodes.size(); 
                    }

                    /** @brief Adds a node into the community.
                     *  @param nodeId The node to add.*/
                    inline void Add( long nodeId ) {
                        m_Nodes.insert( nodeId );
                        m_Cover.m_NodeCommunities[nodeId]->insert(m_Id);
                        const std::set<long>& communities = m_Cover.NodeCommunities(nodeId);
                        for(auto it = communities.begin(); it != communities.end(); ++it ) {
                            if( *it != m_Id ) {
                                m_Cover.m_CommunityMatrix.Inc( m_Id, *it );
                            }
                        }
                    }

                    /** @brief Removes a node from the community.
                     *  @param nodeId The node to remove.*/
                    inline void  Remove( long nodeId ) {
                        m_Nodes.erase( nodeId );
                        m_Cover.m_NodeCommunities[nodeId]->erase(m_Id);
                        const std::set<long>& communities = m_Cover.NodeCommunities(nodeId);
                        for(auto it = communities.begin(); it != communities.end(); ++it ) {
                            if( *it != m_Id ) {
                                m_Cover.m_CommunityMatrix.Dec( m_Id, *it );
                            }
                        }
                    }

                    /** @brief Gets the Id of the community.
                     *  @return The id of the community.*/
                    inline  long    Id() {
                       return m_Id; 
                    }

                    inline const std::set< long>& Nodes() const {
                        return static_cast<const std::set< long> &>(m_Nodes);
                    }

                private:  
                    friend class Cover;
                    Community( Cover & cover,  long id );
                    std::set<long> m_Nodes;
                    Cover&         m_Cover;
                    const long     m_Id;
            };

            Cover( const Graph& graph );
            ~Cover();

            /** @brief  Retrieve the community identifiers where a node belongs to.
             *  @params nodeId The node to retrieve the communities from.
             *  @return a set with the community identifiers.*/  
            const std::set< long>& NodeCommunities( long nodeId ) const;


            /** @brief  Retrieves a community.
             *  @params communityId The community identifier .
             *  @return the community.*/  
            Community&  GetCommunity(  long communityId );

            /** @brief  Retrieves a community.
             *  @params communityId The community identifier .
             *  @return the community.*/  
            const Community& GetCommunity(  long communityId ) const;

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

            /** @brief Refines the communities in the cover.*/
            void RefineCommunities();

        private:
            friend class Community;
            /** @brief Initializes the cover with an initial assignment of nodes to communities.*/
            void  Initialize();

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
             
            const Graph&                                        m_Graph;            /**< @brief The graph of the cover.*/
            std::vector<Community*>                             m_Communities;      /**< @brief The vector of communities in the graph.*/
            std::vector<std::set<long>*>                        m_NodeCommunities;  /**< @brief The communities each node belong.*/
            SparseMatrix<long>                                  m_CommunityMatrix;  /**< @brief The matrix of community overlapps.*/
    };
}
#endif 
