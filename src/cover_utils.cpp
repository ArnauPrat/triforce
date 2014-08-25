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
#include "graph_utils.h"
#include <algorithm>
#include <cassert>
#include <iterator>
#include <cmath>
#include <list>
#include <set>

#define LOOKAHEAD 5 

namespace triforce {

    static long MaxOverlap( long degree, double overlap ) {
       return static_cast<long>(std::floor(overlap*static_cast<double>(degree)) + 1);
    }

    double Score( const Cover& cover, double alpha, double overlap ) {
        double score = 0.0;
        for( long i = 0; i < cover.m_Graph->GetNumNodes(); ++i ) {
            score+=Score(cover,i,alpha,overlap); 
            if(cover.m_NodeMemberships[i].size() > MaxOverlap(cover.m_Graph->GetDegree(i),overlap)  ) {
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
        const Graph& graph = *cover.m_Graph;
        for( long i = 0; i < cover.m_Communities.size(); ++i ) {
            const std::set<long>& nodes = *cover.m_Communities[i];
            if( nodes.size() > 0 ) {
                for( long n : nodes ) {
                    stream << graph.ReMap(n) << " ";
                } 
                stream << std::endl;
            }
        }
    }

    static void ComputeMembershipStats( Cover& cover, double alpha ) {
        const Graph& graph = *cover.m_Graph;
        for( long i = 0; i < graph.GetNumNodes(); ++i ) {
            long r = 0;
            for( long c : cover.m_NodeMemberships[i]) {
                r+=static_cast<long>(cover.m_Communities[c]->size()-1);
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
                long intersectionSize = static_cast<long>(intersection.size());
                din += static_cast<long>(intersectionSize > 0);
                dinPrima += intersectionSize;
                cover.m_Weights[graph.GetEdgeIndex(i,j)] = intersectionSize;
            }

            double denominator = graph.GetDegree(i) + alpha*( r - (dinPrima)) ;
            cover.m_MembershipStats[i].m_Score = denominator > 0.0 ? din / denominator : 0.0;

            denominator = graph.GetDegree(i) + alpha*( r+1 - (dinPrima+1)) ;
            cover.m_MembershipStats[i].m_ConnectedAdd = denominator > 0 ? (din+1) / denominator : 0.0;
            cover.m_MembershipStats[i].m_ConnectedAdd -= cover.m_MembershipStats[i].m_Score;

            denominator = graph.GetDegree(i) + alpha*( r+1 - dinPrima) ;
            cover.m_MembershipStats[i].m_NotConnectedAdd = denominator > 0 ? (din) / denominator : 0.0;
            cover.m_MembershipStats[i].m_NotConnectedAdd -= cover.m_MembershipStats[i].m_Score;

            denominator = graph.GetDegree(i) + alpha*( r-1 - (dinPrima-1)) ;
            cover.m_MembershipStats[i].m_ConnectedRemove = denominator > 0 ? (din-1) / denominator : 0.0;
            cover.m_MembershipStats[i].m_ConnectedRemove -= cover.m_MembershipStats[i].m_Score;

            denominator = graph.GetDegree(i) + alpha*( r-1 - dinPrima) ;
            cover.m_MembershipStats[i].m_NotConnectedRemove = denominator > 0 ? (din) / denominator : 0.0;
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
                cover.m_Communities.push_back(new std::set<long>());
                long communityId = static_cast<long>(cover.m_Communities.size() - 1);
                cover.m_Communities[communityId]->insert(nodeId);
                cover.m_NodeMemberships[nodeId].insert(communityId);
                const  long* adjacencies = graph.GetNeighbors(nodeId);
                long degree = graph.GetDegree(nodeId);
                for( long j = 0; j < degree; ++j ) { 
                    long neighbor = adjacencies[j];
                    if(!visited[neighbor]) {
                        cover.m_Communities[communityId]->insert(neighbor);
                        cover.m_NodeMemberships[neighbor].insert(communityId);
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

    static bool CompareMovement( const Movement& a, const Movement& b ) {
        return a.m_Improvement < b.m_Improvement;
    }

    double TestRemove( const Cover& cover, const long nodeId, const long community, double alpha, double overlapp) {
        double increment = 0.0;
        long dinLost = 0;
        long dinPrimaLost = 0;
        for (long n : *cover.m_Communities[community]) {
            long pos = cover.m_Graph->HasNeighbor(nodeId, n);
            if (pos != -1) {
                long edgeIndex = cover.m_Graph->GetEdgeIndex(nodeId, pos);
                if (cover.m_Weights[edgeIndex] == 1) {
                    increment += cover.m_MembershipStats[n].m_ConnectedRemove;
                    dinLost++;
                }
                dinPrimaLost++;
            }
            else {
                increment += cover.m_MembershipStats[n].m_NotConnectedRemove;
            }
        }
        double denominator = cover.m_Graph->GetDegree(nodeId) + alpha*(cover.m_MembershipStats[nodeId].m_R - cover.m_Communities[community]->size() - (cover.m_MembershipStats[nodeId].m_DinPrima - dinPrimaLost));
        double newScore = denominator > 0 ? (cover.m_MembershipStats[nodeId].m_Din - dinLost) / denominator : 0.0;
        increment += newScore - cover.m_MembershipStats[nodeId].m_Score;
        return increment;
    }

    double TestInsert( const Cover& cover, const long nodeId, const long community, double alpha, double overlapp) {
            double increment = 0.0;
            long dinAdded = 0;
            long dinPrimaAdded = 0;
            for( long n : *cover.m_Communities[community] ) {
                long pos = cover.m_Graph->HasNeighbor(nodeId,n);
                if( pos != -1 ) {
                    long edgeIndex = cover.m_Graph->GetEdgeIndex(nodeId, pos);
                    if( cover.m_Weights[edgeIndex] == 0 ) {
                        increment+=cover.m_MembershipStats[n].m_ConnectedAdd;
                        dinAdded++;
                    }
                    dinPrimaAdded++;
                } else {
                    increment+=cover.m_MembershipStats[n].m_NotConnectedAdd;
                }
            }
            double denominator = cover.m_Graph->GetDegree(nodeId) + alpha*(cover.m_MembershipStats[nodeId].m_R+cover.m_Communities[community]->size() - (cover.m_MembershipStats[nodeId].m_DinPrima + dinPrimaAdded));
            double newScore = denominator > 0 ?  (cover.m_MembershipStats[nodeId].m_Din + dinAdded) / denominator : 0.0;
            increment += newScore - cover.m_MembershipStats[nodeId].m_Score;
        return increment;
    }


    static void PerformMovements( Cover& cover, const long nodeId, double alpha, double overlap, double bestScore, std::vector<Movement>& retRemoves, std::vector<Movement>& retInserts, std::vector<Movement>& retTransfers ) {
        std::vector<Movement> removes;
        const long* adjacencies = cover.m_Graph->GetNeighbors(nodeId);
        for( long c : cover.m_NodeMemberships[nodeId] ) {
            double improvement = TestRemove(cover, nodeId, c, alpha, overlap);
            if( removes.size() == 0 ) {
                Movement movement;
                movement.m_Type = E_REMOVE;
                movement.m_NodeId = nodeId;
                movement.m_CommunityId1 = c;
                movement.m_Improvement = improvement;
                removes.push_back(movement);
            } else if( removes[0].m_Improvement < improvement ) {
                removes[0].m_NodeId = nodeId;
                removes[0].m_CommunityId1 = c;
                removes[0].m_Improvement = improvement;
            }
        }

        std::vector<Movement> inserts;
        std::set<long> candidates;
        long degree = cover.m_Graph->GetDegree(nodeId);
        for( long i = 0; i < degree; ++i ) {
            long neighbor = adjacencies[i];
            for( long c : cover.m_NodeMemberships[neighbor] ) {
                candidates.insert(c);
            }
        }

        for( long c : candidates ) {
            double improvement = TestInsert(cover, nodeId, c, alpha, overlap);
            if( improvement > 0.0 ) {
                if( inserts.size() == 0 ) {
                    Movement movement;
                    movement.m_Type = E_INSERT;
                    movement.m_NodeId = nodeId;
                    movement.m_CommunityId1 = c;
                    movement.m_Improvement = improvement;
                    inserts.push_back(movement);
                } else if( inserts[0].m_Improvement < improvement ) {
                    inserts[0].m_NodeId = nodeId;
                    inserts[0].m_CommunityId1 = c;
                    inserts[0].m_Improvement = improvement;
                }
            }
        }

        std::sort(removes.begin(), removes.end(), CompareMovement);
        std::sort(inserts.begin(), inserts.end(), CompareMovement);
        long rindex = 0, iindex = 0;
        long slots = MaxOverlap(degree,overlap) - static_cast<long>(cover.m_NodeMemberships[nodeId].size());
        while( ( (rindex < removes.size()) && (removes[rindex].m_Improvement != 0.0)) || (iindex < inserts.size() && (inserts[iindex].m_Improvement != 0.0))) {
            bool movementFound = false;
            if( rindex < removes.size() && removes[rindex].m_Improvement > 0.0 ) {
                retRemoves.push_back(removes[rindex]);
                movementFound = true;
                rindex++;
                slots++;
            }
            if( iindex < inserts.size() ) {
                if( slots > 0 && inserts[iindex].m_Improvement > 0.0 ) {
                    retInserts.push_back(inserts[iindex]);
                    movementFound = true;
                    iindex++;
                    slots--;
                } else if ( (rindex < removes.size()) && (inserts[iindex].m_Improvement > 0.0) && (removes[rindex].m_Improvement <= 0.0) &&
                            (inserts[iindex].m_Improvement + removes[rindex].m_Improvement) > 0.0 ) {
                    Movement movement;
                    movement.m_Type = E_TRANSFER;
                    movement.m_NodeId = nodeId;
                    movement.m_CommunityId1 = removes[rindex].m_CommunityId1;
                    movement.m_CommunityId2 = inserts[iindex].m_CommunityId1;
                    movement.m_Improvement = inserts[iindex].m_Improvement + removes[rindex].m_Improvement;
                    retTransfers.push_back(movement);
                    movementFound = true;
                    iindex++;
                    rindex++;
                }
            }
            if(!movementFound) return;
        } 
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
            std::vector<Movement> removes;
            std::vector<Movement> inserts;
            std::vector<Movement> transfers;
            for( long i = 0; i < cover.m_Graph->GetNumNodes(); ++i ) {
                PerformMovements(cover,i, alpha, overlap, currentScore, removes, inserts, transfers);
            }
            for( Movement& movement : removes ) {
                cover.m_Communities[movement.m_CommunityId1]->erase(movement.m_NodeId);
                cover.m_NodeMemberships[movement.m_NodeId].erase(movement.m_CommunityId1);
            }

            for( Movement& movement : inserts ) {
                cover.m_Communities[movement.m_CommunityId1]->insert(movement.m_NodeId);
                cover.m_NodeMemberships[movement.m_NodeId].insert(movement.m_CommunityId1);
            }

            for( Movement& movement : transfers ) {
                cover.m_Communities[movement.m_CommunityId1]->erase(movement.m_NodeId);
                cover.m_NodeMemberships[movement.m_NodeId].erase(movement.m_CommunityId1);
                cover.m_Communities[movement.m_CommunityId2]->insert(movement.m_NodeId);
                cover.m_NodeMemberships[movement.m_NodeId].insert(movement.m_CommunityId2);
            }

            std::list<long> slots;
            for( int i = 0; i < cover.m_Communities.size(); ++i ) {
                if( cover.m_Communities[i]->size() == 0 ) {
                    slots.push_back(i);
                }
            }

            for( int i = 0; i < cover.m_Graph->GetNumNodes(); ++i ) {
                if( cover.m_NodeMemberships[i].size() == 0 ) {
                    if( slots.size() > 1 ) {
                        long community = slots.front();
                        slots.pop_front();
                        cover.m_Communities[community]->insert(i);
                        cover.m_NodeMemberships[i].insert(community);
                    } else {
                        cover.m_Communities.push_back(new std::set<long>());
                        long community = static_cast<long>(cover.m_Communities.size() - 1);
                        cover.m_Communities[community]->insert(i);
                        cover.m_NodeMemberships[i].insert(community);
                    }
                }
            }

            std::cout << "Number of movements performed: " << inserts.size()+removes.size()+transfers.size() << std::endl;
            ComputeMembershipStats(cover,alpha);
            currentScore = Score(cover, alpha, overlap );
            std::cout << "New Score " << currentScore << std::endl;
            double diff = currentScore - bestScore;
            if( diff > 0.001) {
                std::cout << "Cover Improved" << std::endl;
                Destroy(bestCover);
                bestCover = Create(cover.m_Graph);
                Copy(*bestCover,cover);
                bestScore = currentScore;
                lookahead = LOOKAHEAD;
            } 
        }        
        Copy(cover,*bestCover);
    }
}
