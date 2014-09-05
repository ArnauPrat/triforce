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
#include <cstring>
#include <cmath>
#include <list>
#include <set>

#define LOOKAHEAD 5 

namespace triforce {

   struct CommunityRelation {
       long     m_CommunityId1;
       long     m_CommunityId2;

       bool operator() (const CommunityRelation& a, const CommunityRelation& b) const
       {
           if ( a.m_CommunityId1 < b.m_CommunityId1 ) return true;
           if ( b.m_CommunityId1 < a.m_CommunityId1 ) return false;
           return a.m_CommunityId2 < b.m_CommunityId2;
       }
   };

   struct CommunityRelationStats {
       long     m_Connected1;
       long     m_Connected2;
       long     m_DoutConnected1;
       long     m_DoutConnected2;
       long     m_DinConnected1;
       long     m_DinConnected2;
       long     m_Edges;
   };

   std::map<CommunityRelation, CommunityRelationStats, CommunityRelation>  communityRelations;

    static long MaxOverlap( long degree, double overlap ) {
       return static_cast<long>(std::floor(overlap*static_cast<double>(degree)) + 1);
    }

    double Score( const Cover& cover, double alpha, double overlap ) {
        double score = 0.0;
        for( long i = 0; i < cover.m_Graph->GetNumNodes(); ++i ) {
            score+=Score(cover,i,alpha,overlap); 
            if(static_cast<long>(cover.m_NodeMemberships[i].size()) > MaxOverlap(cover.m_Graph->GetDegree(i),overlap)  ) {
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
        for( unsigned long i = 0; i < cover.m_Communities.size(); ++i ) {
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
        if (cover.m_CommunityStats) delete[] cover.m_CommunityStats;
        cover.m_CommunityStats = new CommunityStats[cover.m_Communities.size()];
        std::memset(cover.m_CommunityStats, 0, sizeof(CommunityStats)*cover.m_Communities.size());
        const Graph& graph = *cover.m_Graph;
        for( long i = 0; i < graph.GetNumNodes(); ++i ) {
            long r = 0;
            for( long c : cover.m_NodeMemberships[i]) {
                r+=static_cast<long>(cover.m_Communities[c]->size()-1);
            }

            std::map<long, long> candidates;
            std::map<long, long> inDegrees;
            long din = 0;
            long dinPrima = 0;
            const long* adjacencies = graph.GetNeighbors(i);
            long degree = graph.GetDegree(i);
            for (long j = 0; j < degree; ++j) {
                long neighbor = adjacencies[j];
                std::set<long> intersection;
                std::set_intersection(cover.m_NodeMemberships[i].begin(),
                        cover.m_NodeMemberships[i].end(),
                        cover.m_NodeMemberships[neighbor].begin(),
                        cover.m_NodeMemberships[neighbor].end(),
                        std::inserter(intersection, intersection.begin()));
                long intersectionSize = static_cast<long>(intersection.size());

                if (intersectionSize > 0) {
                    din++;
                    for (long c : intersection) {
                        auto p = inDegrees.insert(std::pair<long,long>(c,0));
                        cover.m_CommunityStats[c].m_InDegree++;
                        (*p.first).second++;
                    }

                }

                if( intersectionSize == 0 ) {
                    for (long c : cover.m_NodeMemberships[neighbor]) {
                        auto p = candidates.insert(std::pair<long,long>(c,0));
                        (*p.first).second++;
                    }
                }

                for (long c : cover.m_NodeMemberships[i]) {
                    if( intersection.find(c) == intersection.end()) {
                        cover.m_CommunityStats[c].m_OutDegree++;
                    }
                }
                dinPrima += intersectionSize;
                cover.m_Weights[graph.GetEdgeIndex(i,j)] = intersectionSize;
            }

            for( auto c : candidates ) {
                for( long r : cover.m_NodeMemberships[i] ) {
                    CommunityRelation rel;
                    rel.m_CommunityId1 = std::min(c.first, r);
                    rel.m_CommunityId2 = std::max(c.first, r);
                    CommunityRelationStats stats;
                    stats.m_Connected1 = 0;
                    stats.m_Connected2 = 0;
                    stats.m_DoutConnected1 = 0;
                    stats.m_DoutConnected2 = 0;
                    stats.m_DinConnected1 = 0;
                    stats.m_DinConnected2 = 0;
                    stats.m_Edges = 0;
                    auto p = communityRelations.insert(std::pair<CommunityRelation,CommunityRelationStats>(rel,stats));
                    auto ref = p.first;
                    (*ref).second.m_Edges+=c.second;
                    if( r == rel.m_CommunityId1 ) {
                        (*ref).second.m_Connected1++;
                        long inDegree = 0;
                        auto p2 = inDegrees.find(r);
                        if( p2 != inDegrees.end() ) {
                            inDegree = (*p2).second;
                        }
                        (*ref).second.m_DinConnected1+= inDegree;
                        (*ref).second.m_DoutConnected1+=cover.m_Graph->GetDegree(i) - inDegree;
                    } else {
                        (*ref).second.m_Connected2++;
                        long inDegree = 0;
                        auto p2 = inDegrees.find(r);
                        if( p2 != inDegrees.end() ) {
                            inDegree = (*p2).second;
                        }
                        (*ref).second.m_DinConnected2+= inDegree;
                        (*ref).second.m_DoutConnected2+=cover.m_Graph->GetDegree(i) - inDegree;
                    }

                    assert(r!=c.first);
                }
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
        std::map<long, long> candidates;
        long degree = cover.m_Graph->GetDegree(nodeId);
        for( long i = 0; i < degree; ++i ) {
            long neighbor = adjacencies[i];
            for( long c : cover.m_NodeMemberships[neighbor] ) {
                if( cover.m_NodeMemberships[nodeId].find(c) == cover.m_NodeMemberships[nodeId].end() ) {
                    auto p = candidates.insert(std::pair<long,long>(c,0));
                    (*p.first).second++;
                }
            }
        }

        for( auto c : candidates ) {
            double improvement = TestInsert(cover, nodeId, c.first, alpha, overlap);
            if( improvement > 0.0 ) {
                if( inserts.size() == 0 ) {
                    Movement movement;
                    movement.m_Type = E_INSERT;
                    movement.m_NodeId = nodeId;
                    movement.m_CommunityId1 = c.first;
                    movement.m_Improvement = improvement;
                    inserts.push_back(movement);
                } else if( inserts[0].m_Improvement < improvement ) {
                    inserts[0].m_NodeId = nodeId;
                    inserts[0].m_CommunityId1 = c.first;
                    inserts[0].m_Improvement = improvement;
                }
            }
        }

        std::sort(removes.begin(), removes.end(), CompareMovement);
        std::sort(inserts.begin(), inserts.end(), CompareMovement);
        long rindex = 0, iindex = 0;
        long slots = MaxOverlap(degree,overlap) - static_cast<long>(cover.m_NodeMemberships[nodeId].size());
        while( ( (rindex < static_cast<long>(removes.size())) && (removes[rindex].m_Improvement != 0.0)) || (iindex < static_cast<long>(inserts.size()) && (inserts[iindex].m_Improvement != 0.0))) {
            bool movementFound = false;
            if( rindex < static_cast<long>(removes.size()) && removes[rindex].m_Improvement > 0.0 ) {
                retRemoves.push_back(removes[rindex]);
                movementFound = true;
                rindex++;
                slots++;
            }
            if( iindex < static_cast<long>(inserts.size()) ) {
                if( slots > 0 && inserts[iindex].m_Improvement > 0.0 ) {
                    retInserts.push_back(inserts[iindex]);
                    movementFound = true;
                    iindex++;
                    slots--;
                } else if ( (rindex < static_cast<long>(removes.size())) && (inserts[iindex].m_Improvement > 0.0) && (removes[rindex].m_Improvement <= 0.0) &&
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


    double TestMerge( const Cover& cover, const long communityId1, const long communityId2, double alpha, double overlap ) {
        const Graph& graph = *cover.m_Graph;
        const std::set<long>& community1 = *cover.m_Communities[communityId1];
        const std::set<long>& community2 = *cover.m_Communities[communityId2];
        double improvement = 0.0;
        for( long nodeId1 : community1 ) {
            double din = 0;
            double dinPrima = 0;
            for( long nodeId2 : community2 ) {
                long pos = graph.HasNeighbor( nodeId1, nodeId2 );
                if( pos != -1 ) {
                    long index = graph.GetEdgeIndex( nodeId1, pos );
                    if( cover.m_Weights[index] == 0 ) {
                        din++;
                    }
                    dinPrima++;
                }
            }
            double denominator = cover.m_Graph->GetDegree(nodeId1) + alpha*(cover.m_MembershipStats[nodeId1].m_R + cover.m_Communities[communityId2]->size() - (cover.m_MembershipStats[nodeId1].m_DinPrima + dinPrima));
            double newScore = denominator > 0 ? (cover.m_MembershipStats[nodeId1].m_Din + din) / denominator : 0.0;
            improvement += newScore - cover.m_MembershipStats[nodeId1].m_Score;
        }

        for( long nodeId1 : community2 ) {
            double din = 0;
            double dinPrima = 0;
            for( long nodeId2 : community1 ) {
                long pos = graph.HasNeighbor( nodeId1, nodeId2 );
                if( pos != -1 ) {
                    long index = graph.GetEdgeIndex( nodeId1, pos );
                    if( cover.m_Weights[index] == 0 ) {
                        din++;
                    }
                    dinPrima++;
                }
            }
            double denominator = cover.m_Graph->GetDegree(nodeId1) + alpha*(cover.m_MembershipStats[nodeId1].m_R + cover.m_Communities[communityId1]->size() - (cover.m_MembershipStats[nodeId1].m_DinPrima + dinPrima));
            double newScore = denominator > 0 ? (cover.m_MembershipStats[nodeId1].m_Din + din) / denominator : 0.0;
            improvement += newScore - cover.m_MembershipStats[nodeId1].m_Score;
        }
        return improvement;
    }

    struct CommunityMerge {
        long m_CommunityId1;
        long m_CommunityId2;
        long m_Improvement;
    };

    static bool CompareCommunityMerge( const CommunityMerge& a, const CommunityMerge& b ) {
        return a.m_Improvement < b.m_Improvement;
    }


   double TryMerge(const Cover& cover, const CommunityRelation& rel, const CommunityRelationStats& stats, double alpha) {
       double r1 = cover.m_Communities[rel.m_CommunityId1]->size();
       double r2 = cover.m_Communities[rel.m_CommunityId2]->size();
       double numConnected1 = stats.m_Connected1;
       double numConnected2 = stats.m_Connected2;
       double numNotConnected1 = r1 - numConnected1;
       double numNotConnected2 = r2 - numConnected2;
       double dinConnected1 = stats.m_DinConnected1 / numConnected1;
       double dinConnected2 = stats.m_DinConnected2 / numConnected2;
       double dinNotConnected1 = numNotConnected1 > 0 ? (cover.m_CommunityStats[rel.m_CommunityId1].m_InDegree - stats.m_DinConnected1) / numNotConnected1 : 0; 
       double dinNotConnected2 = numNotConnected2 > 0 ? (cover.m_CommunityStats[rel.m_CommunityId2].m_InDegree - stats.m_DinConnected2) / numNotConnected2 : 0; 
       double d1 = stats.m_Edges/(numConnected1*2.0f);
       double d2 = stats.m_Edges/(numConnected2*2.0f);
       double doutConnected1 = numConnected1 > 0 ? stats.m_DoutConnected1/numConnected1 : 0; 
       double doutConnected2 = numConnected2 > 0 ? stats.m_DoutConnected2/numConnected2 : 0; 
       double doutNotConnected1 = numNotConnected1 > 0 ? (cover.m_CommunityStats[rel.m_CommunityId1].m_OutDegree - stats.m_DoutConnected1)/numNotConnected1 : 0; 
       double doutNotConnected2 = numNotConnected2 > 0 ? (cover.m_CommunityStats[rel.m_CommunityId2].m_OutDegree - stats.m_DoutConnected1)/numNotConnected2 : 0; 

       double before = numNotConnected1*( dinNotConnected1 / (dinNotConnected1+doutNotConnected1+alpha*(r1-1-dinNotConnected1))) + 
                       numNotConnected2*( dinNotConnected2 / (dinNotConnected2+doutNotConnected2+alpha*(r2-1-dinNotConnected2))) +
                       numConnected1*( dinConnected1 / (dinConnected1+doutConnected1+alpha*(r1-1-dinConnected1))) + 
                       numConnected2*( dinConnected2 / (dinConnected2+doutConnected2+alpha*(r2-1-dinConnected2)));

       double after =  numNotConnected1*( dinNotConnected1 / (dinNotConnected1+doutNotConnected1+alpha*(r1+r2-1-dinNotConnected1))) + 
                       numNotConnected2*( dinNotConnected2 / (dinNotConnected2+doutNotConnected2+alpha*(r1+r2-1-dinNotConnected2))) +
                       numConnected1*( (dinConnected1+d1) / (dinConnected1+doutConnected1+alpha*(r1+r2-1-dinConnected1-d1))) + 
                       numConnected2*( (dinConnected2+d2) / (dinConnected2+doutConnected2+alpha*(r1+r2-1-dinConnected2-d2)));

       double res = after - before;
       if( res > 0.0 ) {
           std::cout << "r1: " << r1 << std::endl;
           std::cout << "r2: " << r2 << std::endl;
           std::cout << "numConnected1: " << numConnected1 << std::endl;
           std::cout << "numConnected2: " << numConnected2 << std::endl;
           std::cout << "numNotConnected1: " << numNotConnected1 << std::endl;
           std::cout << "numNotConnected2: " << numNotConnected2 << std::endl;
           std::cout << "dinConnected1: " << dinConnected1 << std::endl;
           std::cout << "dinConnected2: " << dinConnected2 << std::endl;
           std::cout << "dinNotConnected1: " << dinNotConnected1 << std::endl;
           std::cout << "dinNotConnected2: " << dinNotConnected2 << std::endl;
           std::cout << "d1: " << d1 << std::endl;
           std::cout << "d2: " << d2 << std::endl;
           std::cout << "doutConnected1: " << doutConnected1 << std::endl;
           std::cout << "doutConnected2: " << doutConnected2 << std::endl;
           std::cout << "doutNotConnected1: " << doutNotConnected1 << std::endl;
           std::cout << "doutNotConnected2: " << doutNotConnected2 << std::endl;

           std::cout << "Community1:";
           for( long c : *cover.m_Communities[rel.m_CommunityId1] ) {
               std::cout << " " << cover.m_Graph->ReMap(c);
           }
           std::cout << std::endl;

           std::cout << "Community2:";
           for( long c : *cover.m_Communities[rel.m_CommunityId2] ) {
               std::cout << " " << cover.m_Graph->ReMap(c);
           }
           std::cout << std::endl;
       }
       return  res;
   }

    static void MergeCommunities( Cover& cover, double alpha, double overlapp ) {
        long filtered = 0;
        long total = 0;
        long improved = 0;
        std::vector<CommunityMerge> merges;
        for( auto rel : communityRelations ) {
            long c1 = rel.first.m_CommunityId1;
            long c2 = rel.first.m_CommunityId2;
            if( cover.m_Communities[c1]->size() != 0 && cover.m_Communities[c2]->size() != 0 ) {
                total++;
                double improvement = TryMerge(cover, rel.first, rel.second, alpha) ;
                if( improvement > 0.0 ) {
                    improvement = TestMerge(cover, c1, c2, alpha, overlapp);
                    if (improvement > 0.0) {
                        std::cout << "IMPROVED" << std::endl;
                        CommunityMerge merge;
                        merge.m_CommunityId1 = c1;
                        merge.m_CommunityId2 = c2;
                        merge.m_Improvement = improvement;
                        merges.push_back(merge);
                        improved++;
                    }
                }
                else {
                    filtered++;
                }
            }
        }
        std::cout << "Number of merges filtered by heuristic: " << filtered << " out of " << total << std::endl;
        std::cout << "Heuristic precision: " << improved / static_cast<double>(total - filtered) << std::endl;
        std::cout << "Number of possible merges: " << merges.size() << std::endl;
        std::sort(merges.begin(), merges.end(), CompareCommunityMerge );
        std::set<long> touched;
        for( CommunityMerge m : merges ) {
            if( touched.find(m.m_CommunityId1) == touched.end() &&
                    touched.find(m.m_CommunityId2) == touched.end() ) {
                std::set<long>& community1 = *cover.m_Communities[m.m_CommunityId1];
                std::set<long>& community2 = *cover.m_Communities[m.m_CommunityId2];
                for( long n : community2 ) {
                    cover.m_NodeMemberships[n].erase(m.m_CommunityId2);
                    cover.m_NodeMemberships[n].insert(m.m_CommunityId1);
                }
                community1.insert(community2.begin(), community2.end());
                community2.clear();
                touched.insert(m.m_CommunityId1);
                touched.insert(m.m_CommunityId2);
            }
        }
    }

    void RefineCommunities( Cover& cover, double alpha, double overlap ) {
        Cover* bestCover = Create(cover.m_Graph);
        Copy(*bestCover, cover);
        double bestScore = Score(*bestCover, alpha, overlap );
        double currentScore = bestScore;
        std::cout << "Initial Score: " << bestScore << std::endl;
        int lookahead = LOOKAHEAD;
        while(lookahead > 0) {
            lookahead--;
            communityRelations.clear();
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
            for( int i = 0; i < static_cast<int>(cover.m_Communities.size()); ++i ) {
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
            if( diff <= 0.001) {
                std::cout << "Movements did not improve, trying merges " << std::endl;
                MergeCommunities(cover, alpha, overlap);
                ComputeMembershipStats(cover,alpha);
                currentScore = Score(cover, alpha, overlap );
                std::cout << "New Score after merging " << currentScore << std::endl;
                diff = currentScore - bestScore;
            } 

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
