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
#include "cover_utils.h"
#include "graph.h"
#include "graph_utils.h"
#include <algorithm>
#include <tuple>

#define LOOKAHEAD 5 

namespace triforce {

    Cover::Community::Community( Cover& cover,  long id ) :
        m_Cover(cover),
        m_Id(id),
        m_InternalEdges(0)
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
        Clear();
        for(long i = 0; i < m_Graph.GetNumNodes(); ++i ) { 
            delete m_NodeCommunities[i];
        }
    }

    const std::set< long> Cover::NodeCommunities(  long nodeId ) const {
        return static_cast<std::set< long>&>(*m_NodeCommunities[nodeId]);
    }

    long Cover::NumNodeCommunities( long nodeId ) const {
        return m_NodeCommunities[nodeId]->size();
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
        //return std::tie( std::get<1>(a), std::get<2>(a) ) < std::tie( std::get<1>(b), std::get<2>(b) ) ;
        //return std::tie( std::get<2>(a), std::get<1>(a) ) < std::tie( std::get<2>(b), std::get<1>(b) ) ;
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
                const  long* adjacencies = m_Graph.GetNeighbors(nodeId);
                long degree = m_Graph.GetDegree(nodeId);
                for( long j = 0; j < degree; ++j ) { 
                    if( !visited[adjacencies[j]] ){
                        visited[adjacencies[j]] = true;
        //                community->Add(adjacencies[j]);
                    }
                }   
                nextLabel++;
            }   
        }   
    }

    void Cover::Clear() {
        for(long i = 0; i < m_Graph.GetNumNodes(); ++i ) { 
            m_NodeCommunities[i]->clear();
        }
        m_Communities.clear();
        m_CommunityMatrix.Clear();
    }

    long* Cover::Serialize( long& arraySize ) {
        arraySize = 0;
        long size = m_Communities.size();
        for( long i = 0; i < size; ++i ) {
            arraySize+=m_Communities[i]->Size()+1;
        }

        long* array = new long[arraySize];
        long j = 0;
        for( long i = 0; i < size; ++i ) {
            array[j++] = m_Communities[i]->Size();
            const std::set<long>& nodes = m_Communities[i]->Nodes();
            for( auto it = nodes.begin(); it != nodes.end(); ++it, ++j ) {
                array[j] = *it;
            }
        }
        return array;
    }

    void Cover::Deserialize( long* array , long size) {
        Clear();
        long nextLabel = 0;
        for( long i = 0; i < size; ) {
            Community* community = new Community(*this, nextLabel++);
            m_Communities.push_back(community);
            long limit = i + array[i] + 1;
            i++;
            for(; i < limit; ++i) {
                community->Add(array[i]);
            }  
        } 
    }

    bool Cover::PerformBestMovement( const long nodeId, double alpha, double overlapp, double bestScore, Movement& movement ) {
        bool movementFound = false;
        double bestImprovement = 0.0;
        std::set<long> nodeCommunities = NodeCommunities( nodeId );
        for(auto it = nodeCommunities.begin(); it != nodeCommunities.end(); ++it) {
            Cover::Community& community = this->GetCommunity(*it);
            community.Remove(nodeId);    
            double newScore = Score(*this, alpha, overlapp );
            community.Add(nodeId);
            double scoreDiff =  newScore - bestScore;
            if( scoreDiff > bestImprovement ) {
                bestImprovement = scoreDiff; 
                movement.m_Type = E_REMOVE;
                movement.m_NodeId = nodeId;
                movement.m_CommunityId1 = community.Id();
                movementFound = true;
            }

            /*double aux = TestRemove( nodeId, community, alpha, overlapp ) / m_Graph.GetNumNodes();
            if( !ValidOverlapp(*this,nodeId, overlapp) || aux > bestImprovement  ) {
                bestImprovement = aux; 
                movement.m_Type = E_REMOVE;
                movement.m_NodeId = nodeId;
                movement.m_CommunityId1 = community.Id();
                movementFound = true;
            }*/
        }

        const long* adjacencies = m_Graph.GetNeighbors(nodeId); 
        long degree = m_Graph.GetDegree(nodeId);
        std::set<long> possibleCommunities;
        for( long i = 0; i < degree; ++i ) {
            long neighbor = adjacencies[i];
            const std::set<long>& neighborCommunities = this->NodeCommunities( neighbor );
            for( auto it = neighborCommunities.begin(); it != neighborCommunities.end(); ++it ) {
                const Community& community = this->GetCommunity(*it);
                if( nodeCommunities.find(*it) == nodeCommunities.end())  {
                    possibleCommunities.insert(*it);
                }
            }
        }

        for( auto it = possibleCommunities.begin(); it != possibleCommunities.end(); ++it ) {
            Cover::Community& community = GetCommunity(*it);
            /*community.Add(nodeId);
            double newScore = Score( *this, alpha, overlapp );
            community.Remove(nodeId);
            double scoreDiff =  newScore - bestScore;
            if( scoreDiff > bestImprovement ) {
                bestImprovement = scoreDiff; 
                movement.m_Type = E_INSERT;
                movement.m_NodeId = nodeId;
                movement.m_CommunityId1 = community.Id();
                movementFound = true;
            }*/
            
            double aux = TestInsert( nodeId, community, alpha, overlapp);
            if( aux > bestImprovement && CanOverlapMore( *this, nodeId, overlapp ) ) {
                bestImprovement = aux; 
                movement.m_Type = E_INSERT;
                movement.m_NodeId = nodeId;
                movement.m_CommunityId1 = community.Id();
                movementFound = true;
            }
        }

        nodeCommunities = NodeCommunities( nodeId );
        for(auto it = nodeCommunities.begin(); it != nodeCommunities.end(); ++it) {
            Cover::Community& community1 = this->GetCommunity(*it);
            for( auto it2 = possibleCommunities.begin(); it2 != possibleCommunities.end(); ++it2 ) {
                Cover::Community& community2 = this->GetCommunity(*it2);
                if( (*it != *it2) && community1.Contains(nodeId) && !community2.Contains(nodeId) ) {
                    community1.Remove(nodeId);
                    community2.Add(nodeId);
                    double newScore = Score( *this, alpha, overlapp );
                    community1.Add(nodeId);
                    community2.Remove(nodeId);
                    double scoreDiff = newScore - bestScore ;
                    if( scoreDiff > bestImprovement) {
                        bestImprovement = newScore - bestScore; 
                        movement.m_Type = E_TRANSFER;
                        movement.m_NodeId = nodeId;
                        movement.m_CommunityId1 = community1.Id();
                        movement.m_CommunityId2 = community2.Id();
                        movementFound = true;
                    }
                    

                /*    double aux = TestInsert( nodeId, community2, alpha, overlapp ) + TestRemove( nodeId, community1, alpha, overlapp );
                    if( aux > bestImprovement && CanOverlapMore(*this,nodeId, overlapp) ) {
                        bestImprovement = aux; 
                        movement.m_Type = E_TRANSFER;
                        movement.m_NodeId = nodeId;
                        movement.m_CommunityId1 = community1.Id();
                        movement.m_CommunityId2 = community2.Id();
                        movementFound = true;
                    }
                    */
                }
            }
        }
        return movementFound;
    }

void Cover::RefineCommunities( double alpha, double overlapp ) {
        long bestCoverSize = 0;
        long* bestCover = Serialize( bestCoverSize );
        double bestScore = Score(*this, alpha, overlapp );
        double currentScore = bestScore;
        std::cout << "Initial Score: " << bestScore << std::endl;
        int lookahead = LOOKAHEAD;
        while(lookahead > 0) {
            lookahead--;
            std::cout << "Starting Iteration" << std::endl;
            //Compute new cover.
            std::vector<Movement> movements;
            for( long i = 0; i < m_Graph.GetNumNodes(); ++i ) {
                Movement movement;
                if(PerformBestMovement(i, alpha, overlapp, currentScore, movement)) {
                    movements.push_back(movement);
                }
            }
            std::cout << "Number of movements performed " << movements.size() << std::endl;
            Community* community1; 
            Community* community2;
            for( long i = 0; i < movements.size(); ++i ) {
                Movement& movement = movements[i];
                switch(movement.m_Type) {
                    case E_REMOVE:
                        community1 = &(GetCommunity(movement.m_CommunityId1));
                        community1->Remove(movement.m_NodeId);
                        break;
                    case E_INSERT:
                        community1 = &(GetCommunity(movement.m_CommunityId1));
                        community1->Add(movement.m_NodeId);
                        break;
                    case E_TRANSFER:
                        community1 = &(GetCommunity(movement.m_CommunityId1));
                        community2 = &(GetCommunity(movement.m_CommunityId2));
                        community1->Remove(movement.m_NodeId);
                        community2->Add(movement.m_NodeId);
                        break;
                };
            }
            currentScore = Score(*this, alpha, overlapp );
            std::cout << "New Score " << currentScore << std::endl;
            if( currentScore > bestScore ) {
                std::cout << "Cover Improved" << std::endl;
                delete [] bestCover;
                bestCover = Serialize( bestCoverSize );
                bestScore = currentScore;
                lookahead = LOOKAHEAD;
            }
        }        
        Deserialize(bestCover, bestCoverSize);
        delete [] bestCover;
    }
}

