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
#include "graph_utils.h"
#include "graph.h"
#include <set>
#include <algorithm>
#include <cmath>

#define LOOKAHEAD 5 

namespace triforce {

    double Score( const Cover& cover, double alpha, double overlap ) {
        double score = 0.0;
        for( long i = 0; i < cover.m_Graph->GetNumNodes(); ++i ) {
            score+=Score(cover,i,alpha,overlap); 
            if((cover.m_NodeMemberships[i].size()-1) / static_cast<double>(cover.m_Graph->GetDegree(i)) > overlap ) {
                std::cout << "ENTRA" << std::endl;
                std::cout << i << " " << cover.m_NodeMemberships[i].size() << " " << cover.m_Graph->GetDegree(i) << std::endl;
                return 0.0;
            }
        }
        return score / cover.m_Graph->GetNumNodes();
    }

    double Score( const Cover& cover,  long nodeId, double alpha, double overlap) {
        return cover.m_MembershipStats[nodeId].m_Score;
    }


    void Print( const Cover& cover, std::ostream& stream ) {
        /*        const Graph& graph = cover.GetGraph();
                  for( long i = 0; i < cover.NumCommunities(); ++i ) {
                  const Cover::Community& community = cover.GetCommunity(i);
                  if( community.Size() > 0 ) {
                  const std::set<long>& nodes = community.Nodes();
                  for( auto it = nodes.begin(); it != nodes.end(); ++it ) {
                  stream << graph.ReMap(*it) << " ";
                  } 
                  stream << std::endl;
                  }
                  }
                  */
    }

    double TestRemove( const Cover& cover, const long nodeId, const long community, double alpha, double overlapp) {
        return 0.0; 
    }

    double TestInsert( const Cover& cover, const long nodeId, const long comunity, double alpha, double overlapp) {
        return 0.0;
    }

    static void ComputeMembershipStats( Cover& cover, double alpha ) {
        const Graph& graph = *cover.m_Graph;
        for( long i = 0; i < graph.GetNumNodes(); ++i ) {
            long r = 0;
            for( long c : cover.m_NodeMemberships[i]) {
                r+=cover.m_Communities[c].size()-1;
            }
            long din = 0;
            long dinPrima = 0;
            const long* adjacencies = graph.GetNeighbors(i);
            long degree = graph.GetDegree(i);
            for( long j = 0; j < degree; ++j ) {
                long neighbor = adjacencies[j];
                std::set<long> intersection;
                std::set_intersection(cover.m_NodeMemberships[i].begin(),
                        cover.m_NodeMemberships[i].end(),
                        cover.m_NodeMemberships[neighbor].begin(),
                        cover.m_NodeMemberships[neighbor].end(),
                        std::inserter(intersection, intersection.begin()));
                long intersectionSize = intersection.size();
                din += static_cast<long>(intersectionSize > 0);
                dinPrima += intersectionSize;
                cover.m_Weights[graph.GetEdgeIndex(i,j)] = intersectionSize;
            }

            double denominator = graph.GetDegree(i) + alpha*( r - (dinPrima)) ;
            cover.m_MembershipStats[i].m_Score = denominator > 0 ? (din) / denominator : 0;
            denominator = graph.GetDegree(i) + alpha*( r+1 - (dinPrima+1)) ;
            cover.m_MembershipStats[i].m_ConnectedAdd = denominator > 0 ? (din+1) / denominator : 0;
            cover.m_MembershipStats[i].m_ConnectedAdd -= cover.m_MembershipStats[i].m_Score;
            denominator = graph.GetDegree(i) + alpha*( r+1 - dinPrima) ;
            cover.m_MembershipStats[i].m_NotConnectedAdd = denominator > 0 ? (din) / denominator : 0;
            cover.m_MembershipStats[i].m_NotConnectedAdd -= cover.m_MembershipStats[i].m_Score;
            denominator = graph.GetDegree(i) + alpha*( r-1 - (dinPrima-1)) ;
            cover.m_MembershipStats[i].m_ConnectedRemove = denominator > 0 ? (din-1) / denominator : 0;
            cover.m_MembershipStats[i].m_ConnectedRemove -= cover.m_MembershipStats[i].m_Score;
            denominator = graph.GetDegree(i) + alpha*( r-1 - dinPrima) ;
            cover.m_MembershipStats[i].m_NotConnectedRemove = denominator > 0 ? (din) / denominator : 0;
            cover.m_MembershipStats[i].m_NotConnectedRemove -= cover.m_MembershipStats[i].m_Score;
            cover.m_MembershipStats[i].m_R = r;
            cover.m_MembershipStats[i].m_Din = din;
            cover.m_MembershipStats[i].m_DinPrima = dinPrima;
        } 
    }

    /** @brief Compares two node clustering.
     *  @param a The first tuple to sort.
     *  @param b The second tuple to sort.   
     *  @return -1 if a goes before b. 1 if a goes after b. 0 if a and b are equal.*/
    static bool CompareClusterings(const std::tuple<long, double,long >& a, const std::tuple<long, double, long> b) {
        return std::tie( std::get<1>(a), std::get<2>(a) ) > std::tie( std::get<1>(b), std::get<2>(b) ) ;
    }

    void Initialize( Cover& cover, double alpha ) {
        // Computing clustering coefficients
        const Graph& graph = *cover.m_Graph;
        std::vector< std::tuple<  long, double,  long > > nC( graph.GetNumNodes() ); 
        for(  long i = 0; i < graph.GetNumNodes(); ++i ) { 
            nC.push_back( std::tuple<  long, double,  long >( i,
                        ClusteringCoefficient( graph, i ),
                        graph.GetDegree(i)) 
                    );  
        }   
        std::sort(nC.begin(), nC.end(), CompareClusterings );                     //Sorts the nC
        std::vector<bool> visited( graph.GetNumNodes(), false );                //Vector to track those visited members.

        // Creating initial communities
        for( auto it = nC.begin(); it != nC.end(); ++it ) {    
            long nodeId = std::get<0>(*it);
            if(!visited[nodeId]) {
                visited[nodeId] = true;
                cover.m_Communities[nodeId].insert(nodeId);
                cover.m_NodeMemberships[nodeId].insert(nodeId);
                const  long* adjacencies = graph.GetNeighbors(nodeId);
                long degree = graph.GetDegree(nodeId);
                for( long j = 0; j < degree; ++j ) { 
                    long neighbor = adjacencies[j];
                    if(!visited[neighbor]) {
                        cover.m_Communities[nodeId].insert(neighbor);
                        cover.m_NodeMemberships[neighbor].insert(nodeId);
                        visited[neighbor] = true;
                    }
                }   
            }   
        }   
        ComputeMembershipStats(cover, alpha);
    }

    enum MovementType {
        E_REMOVE,
        E_INSERT,
        E_TRANSFER,
    };

    struct Movement {
        MovementType m_Type;
        long m_NodeId;
        long m_CommunityId1;
        long m_CommunityId2;
        double m_Improvement;
    };

    static std::vector<Movement> PerformMovements( Cover& cover, const long nodeId, double alpha, double overlapp, double bestScore ) {
        std::vector<Movement> movements;
        const long* adjacencies = cover.m_Graph->GetNeighbors(nodeId);
        for( long c : cover.m_NodeMemberships[nodeId] ) {
            double increment = 0.0;
            long dinLost = 0;
            long dinPrimaLost = 0;
            for( long n : cover.m_Communities[c] ) {
                long pos = cover.m_Graph->HasNeighbor(nodeId,n);
                if( pos != -1 ) {
                    long edgeIndex = cover.m_Graph->GetEdgeIndex(nodeId, pos);
                    if( cover.m_Weights[edgeIndex] == 1 ) {
                        increment+=cover.m_MembershipStats[n].m_ConnectedRemove;
                        dinLost++;
                    }
                    dinPrimaLost++;
                } else {
                    increment+=cover.m_MembershipStats[n].m_NotConnectedRemove;
                }
            }
            double denominator = cover.m_Graph->GetDegree(nodeId) + alpha*(cover.m_MembershipStats[nodeId].m_R-cover.m_Communities[c].size() - (cover.m_MembershipStats[nodeId].m_DinPrima - dinPrimaLost));
            double newScore = denominator > 0 ?  (cover.m_MembershipStats[nodeId].m_Din - dinLost) / denominator : 0.0;
            increment += newScore - cover.m_MembershipStats[nodeId].m_Score;
            if( increment > 0.0 ) {
                Movement movement;
                movement.m_Type = E_REMOVE;
                movement.m_NodeId = nodeId;
                movement.m_CommunityId1 = c;
                movements.push_back(movement);
            }
        }

       return movements; 
    }

    void RefineCommunities( Cover& cover, double alpha, double overlap ) {
        long bestCoverSize = 0;
        Cover* bestCover = Create(cover.m_Graph);
        Copy(*bestCover, cover);
        double bestScore = Score(*bestCover, alpha, overlap );
        double currentScore = bestScore;
        std::cout << "Initial Score: " << bestScore << std::endl;
        int lookahead = LOOKAHEAD;
        while(lookahead > 0) {
            lookahead--;
            std::cout << "Starting Iteration" << std::endl;
            std::vector<Movement> movements;
            for( long i = 0; i < cover.m_Graph->GetNumNodes(); ++i ) {
                std::vector<Movement> aux = PerformMovements(cover,i, alpha, overlap, currentScore);
                movements.insert(movements.end(), aux.begin(), aux.end());
            }
            for( long i = 0; i < movements.size(); ++i ) {
                Movement& movement = movements[i];
                switch(movement.m_Type) {
                    case E_REMOVE:
                        cover.m_Communities[movement.m_CommunityId1].erase(movement.m_NodeId);
                        cover.m_NodeMemberships[movement.m_NodeId].erase(movement.m_CommunityId1);
                        /*if(cover.NumCommunities(movement.m_NodeId) == 0) {
                            Cover::Community* community = new Cover::Community(cover,cover.NumCommunities());
                            cover.m_Communities.push_back(community);
                            community->Add(movement.m_NodeId);
                        }
                        */
                        break;
                    case E_INSERT:
                        cover.m_Communities[movement.m_CommunityId1].insert(movement.m_NodeId);
                        cover.m_NodeMemberships[movement.m_NodeId].insert(movement.m_CommunityId1);
                        break;
                    /*case E_TRANSFER:
                        community1 = &(cover.GetCommunity(movement.m_CommunityId1));
                        community2 = &(cover.GetCommunity(movement.m_CommunityId2));
                        community1->Remove(movement.m_NodeId);
                        community2->Add(movement.m_NodeId);
                        break;
                        */
                };
            }
            std::cout << "Number of movements performed: " << movements.size() << std::endl;
            ComputeMembershipStats(cover,alpha);
            currentScore = Score(cover, alpha, overlap );
            std::cout << "New Score " << currentScore << std::endl;
            double diff = currentScore - bestScore;
            if( diff > 0.001) {
                std::cout << "Cover Improved" << std::endl;
                Destroy(bestCover);
                Copy(*bestCover,cover);
                bestScore = currentScore;
                lookahead = LOOKAHEAD;
            } 
        }        
        Copy(cover,*bestCover);
    }
}
