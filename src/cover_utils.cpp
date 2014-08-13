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
#include <cfloat>

#define LOOKAHEAD 5 

namespace triforce {

    bool ValidOverlapp( const Cover& cover,  long nodeId, double overlapp ) {
        /*const std::set< long>& communities = cover.NodeCommunities( nodeId );        
        for( auto it1 = communities.begin(); it1 != communities.end(); ++it1 ) {
            for( auto it2 = communities.begin(); it2 != communities.end(); ++it2 ) {
                if( *it1 != *it2 && Similar(cover, cover.GetCommunity(*it1), cover.GetCommunity(*it2), overlapp)) {
                    return false; 
                }
            }
        }
        return true;
        */
        const Graph& graph = cover.GetGraph();
        long degree = graph.GetDegree(nodeId);
        return ((cover.NumNodeCommunities(nodeId)-1) / (double)degree) <= overlapp;
    }

    bool CanOverlapMore( const Cover& cover, const long nodeId, double overlapp ) {
        const Graph& graph = cover.GetGraph();
        long degree = graph.GetDegree(nodeId);
        return ((cover.NumNodeCommunities(nodeId)) / (double)degree) <= overlapp;
    }

    bool Similar( const Cover& cover, const Cover::Community& communityA, const Cover::Community& communityB, double overlapp ) {
        const std::set< long>& nodesA = communityA.Nodes();
        const std::set< long>& nodesB = communityB.Nodes();
        auto it1 = nodesA.begin();
        auto it2 = nodesB.begin();
         long common = 0;
        while( (it1 != nodesA.end()) && (it2 != nodesB.end()) ) { 
            if( *it1 < *it2 ) { 
               ++it1; 
            }   
            else if( *it1 > *it2 ) { 
                ++it2;
            } else {
                common++;
                ++it1;
                ++it2;
            }   
        }   
        //long common = cover.m_CommunityMatrix.Get(communityA.Id(), communityB.Id());
        return (common / (double)( std::min(communityA.Size(), communityB.Size()))) > overlapp;
    }

    double Score( const Cover& cover,  long nodeId, double alpha, double overlapp, bool& validOverlapp ) {
        const Graph& graph = cover.GetGraph();
        if( ValidOverlapp( cover, nodeId, overlapp ) ) {
            std::set< long> unionSet = CommunityUnion( cover, nodeId );
            const  long* adjacencies = graph.GetNeighbors( nodeId );
             long degree = graph.GetDegree( nodeId );
             long kin = 0;
             long kout = 0;
            for(  long i = 0; i < degree; ++i ) {
                if( unionSet.find(adjacencies[i]) != unionSet.end() ) {
                    kin++;
                } else {
                    kout++;
                }
            }
            long denom = kout + unionSet.size() - 1;
            validOverlapp = true;
            if( denom == 0 ) return 0.0;
            return kin /(double)( kout + kin + alpha*(unionSet.size()-kin - 1));
        }
        validOverlapp = false;
        return 0.0;
    }

    double Score( const Cover& cover, double alpha, double overlapp ) {
        const Graph& graph = cover.GetGraph();
        long numNodes = graph.GetNumNodes();
        double sum = 0;
        for(  long i = 0; i < numNodes; ++i ) {
            bool validOverlapp;
            double aux = Score( cover, i, alpha, overlapp, validOverlapp );
            if(!validOverlapp) return 0.0;
            sum += aux;
        }
        return sum / (double)numNodes;
    }

    std::set< long> CommunityUnion( const Cover& cover, 
                                            const  long nodeId ) {
        std::set< long> unionSet;
        const std::set< long>& communities = cover.NodeCommunities( nodeId );        
        for( auto it1 = communities.begin(); it1 != communities.end(); ++it1 ) {
            const std::set< long>& communityNodes = cover.GetCommunity(*it1).Nodes();
            for( auto it2 = communityNodes.begin(); it2 != communityNodes.end(); ++it2) {
                unionSet.insert(*it2);
            }
        }
        return unionSet;
    }

    void Print( const Cover& cover, std::ostream& stream ) {
        const Graph& graph = cover.GetGraph();
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
    }

    void PrintZero( const Cover& cover, std::ostream& stream, double alpha, double overlapp ) {
        const Graph& graph = cover.GetGraph();
         long numNodes = graph.GetNumNodes();
        double sum = 0;
        for(  long i = 0; i < numNodes; ++i ) {
            bool validOverlapp;
            double aux = Score( cover, i, alpha, overlapp, validOverlapp );
            if( aux == 0.0 ) stream << graph.ReMap(i) << " has score zero" << std::endl;;
        }
    }

    static double CheckForIncrement( long r, long din, long dout, long cout, double pin, double alpha ) {
//        std::cout << "Size: " << r << " din:  "<< din  << " dout: " << dout << " cout: " << cout << " pin " << pin << " alpha " << alpha << std::endl;
        double avgKout;
        double avgKin;
        if (r > 0) {
            avgKout = (cout - din) / (double) r;
            avgKin = (r-1)*pin; 
        } else {
            avgKout = 0.0;
            avgKin = 0.0;
        }
        double A = 0.0;
        double denom1 = avgKout + dout + avgKin + alpha*(r-1-avgKin);
        double denom2 = avgKout + 1 + avgKin + alpha*(r-1-avgKin);
        if (denom1 != 0.0 ) {
            A = (avgKin+1) / denom1;
        }
        if(denom2 != 0.0) {
            A -= (avgKin) / denom2;
        }
        double B = 0.0;
        denom1 = avgKout + avgKin + alpha*(r-avgKin);
        denom2 = avgKout + avgKin + alpha*(r-avgKin -1);
        if (denom1 != 0.0) {
            B = avgKin / denom1;
        }
        if (denom2 != 0.0) {
            B -= avgKin / denom2;
        }

        double C = 0.0;
        denom1 = dout + din + alpha*(r-1);;
        if (denom1 != 0.0 ) {
            C = din/denom1;
        }
        double res = (A*din + (r - din)*B + C);
        return res;
    }

    double TestRemove( const long nodeId, const Cover::Community& community, double alpha, double overlapp, bool& validOverlapp) {
        const Cover& cover = community.GetCover();
        const Graph& graph = community.GetGraph(); 
        Cover::Community& comAux = const_cast<Cover::Community&>(community);
/*        const long* adjacencies = graph.GetNeighbors(nodeId);
        long degree = graph.GetDegree(nodeId);
        long din = 0;
        long dout = 0;
        std::set<long> inNeighbors;
        for( long i = 0; i < degree; ++i ) {
            long neighbor = adjacencies[i];
            if( community.Contains(neighbor) ) {
                din++;
                inNeighbors.insert(neighbor);
            } else {
                dout++;
            }
        }
        */
        double before = 0.0;
        double after = 0.0;
        long r = community.Size();
        bool tempValidOverlapp;
        bool violatesOverlapp = false;
        std::set<long> nodes = community.Nodes();
        for( auto it = nodes.begin(); it != nodes.end(); ++it ) {
            before += Score( cover, *it, alpha, overlapp, tempValidOverlapp );
            violatesOverlapp = violatesOverlapp || !tempValidOverlapp;
        }
        if( violatesOverlapp ) before = 0.0;

        comAux.Remove(nodeId);
        nodes = community.Nodes();
        violatesOverlapp = false;
        for( auto it = nodes.begin(); it != nodes.end(); ++it ) {
            after += Score( cover, *it, alpha, overlapp, tempValidOverlapp );
            violatesOverlapp = violatesOverlapp || !tempValidOverlapp;
        }
        after += Score( cover, nodeId, alpha, overlapp, tempValidOverlapp );
        violatesOverlapp = violatesOverlapp || !tempValidOverlapp;
        if( violatesOverlapp ) {
            validOverlapp = false;
        } else {
            validOverlapp = true;
        }
        comAux.Add(nodeId);
        /*double pin = ((community.GetInternalEdges()-din)*2)/(double)((community.Size()-1)*(community.Size()-2));
        double increment = CheckForIncrement( community.Size() - 1, din, dout, community.GetExternalEdges()-dout+din, pin, alpha);
        */
        //return  (after - before)/graph.GetNumNodes();
        return  (after - before);
    }

    double TestInsert( const long nodeId, const Cover::Community& community, double alpha, double overlapp, bool& validOverlapp ) {
        const Cover& cover = community.GetCover();
        const Graph& graph = community.GetGraph(); 
        Cover::Community& comAux = const_cast<Cover::Community&>(community);
        /*const long* adjacencies = graph.GetNeighbors(nodeId);
        long degree = graph.GetDegree(nodeId);
        long din = 0;
        long dout = 0;
        for( long i = 0; i < degree; ++i ) {
            long neighbor = adjacencies[i];
            if( community.Contains(neighbor) ) {
                din++;
            } else {
                dout++;
            }
        }
        double increment = CheckForIncrement( community.Size(), din, dout, community.GetInternalEdges()*2/(double)((community.Size())*(community.Size() -1)),community.GetExternalEdges(), alpha);
        */

        double before = 0.0;
        double after = 0.0;
        bool tempValidOverlapp;
        bool violatesOverlapp = false;
        std::set<long> nodes = comAux.Nodes();
        for( auto it = nodes.begin(); it != nodes.end(); ++it ) {
            before += Score( cover, *it, alpha, overlapp, tempValidOverlapp );
            violatesOverlapp = violatesOverlapp || !tempValidOverlapp;
        }
        before += Score( cover, nodeId, alpha, overlapp, tempValidOverlapp );
        violatesOverlapp = violatesOverlapp || !tempValidOverlapp;
        if( violatesOverlapp ) before = 0.0;
        
        violatesOverlapp = false;
        comAux.Add(nodeId);
        nodes = comAux.Nodes();
        for( auto it = nodes.begin(); it != nodes.end(); ++it ) {
            after += Score( cover, *it, alpha, overlapp, tempValidOverlapp );
            violatesOverlapp = violatesOverlapp || !tempValidOverlapp;
        }
        comAux.Remove(nodeId);
        if( violatesOverlapp ) {
            validOverlapp = false;
        } else {
            validOverlapp = true;
        }

        //return (after/graph.GetNumNodes()) - (before/graph.GetNumNodes());
        return (after) - (before);
    }

    bool ViolatesOverlappInsert( const long nodeId, const Cover::Community& community) {
/*        const Cover& cover = community.GetCover(); 
        const std::set<long>& nodeCommunities = cover.NodeCommunities( nodeId ); 
        for( auto it = nodeCommunities.begin(); it != nodeCommunities.end(); ++it ) {
            const Cover::Community& communityB = cover.GetCommunity(*it);
            long common = cover.m_CommunityMatrix.Get( community.Id(), communityB.Id() );
            if( (common + 1)/ ((double) community.Size() + communityB.Size() + 1) > OVERLAPPING_THRESHOLD &&
                (common)/ ((double) community.Size() + communityB.Size()) <= OVERLAPPING_THRESHOLD ) {
                return true;
            }
        }
        */
        return false;
    }

    bool ViolatesOverlappRemove( const long nodeId, const Cover::Community& community) {
        /*const Cover& cover = community.GetCover(); 
        const std::set<long>& nodeCommunities = cover.NodeCommunities( nodeId ); 
        for( auto it = nodeCommunities.begin(); it != nodeCommunities.end(); ++it ) {
            const Cover::Community& communityB = cover.GetCommunity(*it);
            long common = cover.m_CommunityMatrix.Get( community.Id(), communityB.Id() );
            if( (common - 1)/ ((double) community.Size() + communityB.Size() - 1) <= OVERLAPPING_THRESHOLD &&
                (common)/ ((double) community.Size() + communityB.Size()) > OVERLAPPING_THRESHOLD ) {
                return false;
            }
        }
        */
        return true;
    }

    /** @brief Compares two node clustering.
     *  @param a The first tuple to sort.
     *  @param b The second tuple to sort.   
     *  @return -1 if a goes before b. 1 if a goes after b. 0 if a and b are equal.*/
    static bool CompareClusterings(const std::tuple<long, double,long >& a, const std::tuple<long, double, long> b) {
        //return std::tie( std::get<1>(a), std::get<2>(a) ) < std::tie( std::get<1>(b), std::get<2>(b) ) ;
        //return std::tie( std::get<2>(a), std::get<1>(a) ) < std::tie( std::get<2>(b), std::get<1>(b) ) ;
        return std::tie( std::get<1>(a), std::get<2>(a) ) > std::tie( std::get<1>(b), std::get<2>(b) ) ;
    }

    void Initialize( Cover& cover, double alpha, double overlapp ) {

        const Graph& graph = cover.GetGraph();
        std::vector< std::tuple<  long, double,  long > > nC( graph.GetNumNodes() ); 
        for(  long i = 0; i < graph.GetNumNodes(); ++i ) { 
            nC.push_back( std::tuple<  long, double,  long >( i,
                        ClusteringCoefficient( graph, i ),
                        graph.GetDegree(i)) 
                    );  
        }   
        std::sort(nC.begin(), nC.end(), CompareClusterings );                     //Sorts the nC
        std::vector<bool> visited( graph.GetNumNodes(), false );                //Vector to track those visited members.
        long nextLabel = 0;
        for( auto it = nC.begin(); it != nC.end(); ++it ) {    
            long nodeId = std::get<0>(*it);
            if(!visited[nodeId]) {
                visited[nodeId] = true;
                Cover::Community* community = new Cover::Community( cover, nextLabel );
                cover.m_Communities.push_back(community);
                community->Add(nodeId);
                const  long* adjacencies = graph.GetNeighbors(nodeId);
                long degree = graph.GetDegree(nodeId);
                    for( long j = 0; j < degree; ++j ) { 
                        long neighbor = adjacencies[j];
                        if(!visited[neighbor]) {
                            community->Add(neighbor);
                            visited[neighbor] = true;
                        }
                    }   
                nextLabel++;
            }   
        }   
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

  bool PerformBestMovement( Cover& cover, const long nodeId, double alpha, double overlapp, double bestScore, Movement& movement ) {
        bool movementFound = false;
        double bestRemoveImprovement = -FLT_MAX;
        long bestRemoveCommunity = -1;
        bool bestRemoveFound = false;
        bool bestRemoveOverlapp = false;
        std::set<long> nodeCommunities = cover.NodeCommunities( nodeId );
        for(auto it = nodeCommunities.begin(); it != nodeCommunities.end(); ++it) {
            Cover::Community& community = cover.GetCommunity(*it);
            bool auxValidOverlapp;
            double aux = TestRemove( nodeId, community, alpha, overlapp, auxValidOverlapp);
            if( aux >= bestRemoveImprovement  ) {
                bestRemoveImprovement = aux; 
                bestRemoveCommunity = community.Id();
                bestRemoveOverlapp = auxValidOverlapp;
                if( bestRemoveImprovement > 0.0 ) bestRemoveFound = true;
            }
        }

        const long* adjacencies = cover.m_Graph.GetNeighbors(nodeId); 
        long degree = cover.m_Graph.GetDegree(nodeId);
        std::set<long> possibleCommunities;
        for( long i = 0; i < degree; ++i ) {
            long neighbor = adjacencies[i];
            const std::set<long>& neighborCommunities = cover.NodeCommunities( neighbor );
            for( auto it = neighborCommunities.begin(); it != neighborCommunities.end(); ++it ) {
                const Cover::Community& community = cover.GetCommunity(*it);
                if( nodeCommunities.find(*it) == nodeCommunities.end())  {
                    possibleCommunities.insert(*it);
                }
            }
        }

        double bestInsertImprovement = -FLT_MAX;
        long bestInsertCommunity = -1;
        bool bestInsertFound = false;
        bool bestInsertOverlapp = false;
        for( auto it = possibleCommunities.begin(); it != possibleCommunities.end(); ++it ) {
            Cover::Community& community = cover.GetCommunity(*it);
            bool auxOverlapp;
            double aux = TestInsert( nodeId, community, alpha, overlapp, auxOverlapp );
            if( aux >= bestInsertImprovement ) {
                bestInsertImprovement = aux; 
                bestInsertCommunity = community.Id();
                bestInsertOverlapp = auxOverlapp;
                if( bestInsertImprovement > 0.0 ) bestInsertFound = true;
            }
        }

        double bestImprovement = 0.0;
        if( bestRemoveFound && bestRemoveOverlapp ) {
            bestImprovement = bestRemoveImprovement;
            movement.m_Type = E_REMOVE;
            movement.m_NodeId = nodeId;
            movement.m_CommunityId1 = bestRemoveCommunity;
            movement.m_Improvement = bestImprovement;
            movementFound = true;
        }
        if( bestInsertFound && bestInsertOverlapp && (bestInsertImprovement > bestRemoveImprovement) ) {
            bestImprovement = bestInsertImprovement;
            movement.m_Type = E_INSERT;
            movement.m_NodeId = nodeId;
            movement.m_CommunityId1 = bestInsertCommunity;
            movement.m_Improvement = bestImprovement;
            movementFound = true;
        }
        double bestRemoveAndInsert = (bestInsertImprovement + bestRemoveImprovement);
        if( (bestRemoveAndInsert > 0.0) &&  (bestRemoveAndInsert > bestImprovement) && ValidOverlapp(cover, nodeId, overlapp)) {
            bestImprovement = bestRemoveAndInsert;
            movement.m_Type = E_TRANSFER;
            movement.m_NodeId = nodeId;
            movement.m_CommunityId1 = bestRemoveCommunity;
            movement.m_CommunityId2 = bestInsertCommunity;
            movement.m_Improvement = bestImprovement;
            movementFound = true;
        } 
        //std::cout << cover.m_Graph.ReMap(nodeId) << " "<< bestRemoveAndInsert << " " << bestInsertImprovement << " " << bestRemoveImprovement << " " << bestImprovement << " " << movement.m_Type << " " << movementFound << std::endl;

//        nodeCommunities = NodeCommunities( nodeId );
//        for(auto it = nodeCommunities.begin(); it != nodeCommunities.end(); ++it) {
//            Cover::Community& community1 = this->GetCommunity(*it);
//            for( auto it2 = possibleCommunities.begin(); it2 != possibleCommunities.end(); ++it2 ) {
//                Cover::Community& community2 = this->GetCommunity(*it2);
//                if( (*it != *it2) && community1.Contains(nodeId) && !community2.Contains(nodeId) ) {
//                /*    community1.Remove(nodeId);
//                    community2.Add(nodeId);
//                    double newScore = Score( *this, alpha, overlapp );
//                    community1.Add(nodeId);
//                    community2.Remove(nodeId);
//                    double scoreDiff = newScore - bestScore ;
//                    if( scoreDiff > bestImprovement) {
//                        bestImprovement = newScore - bestScore; 
//                        movement.m_Type = E_TRANSFER;
//                        movement.m_NodeId = nodeId;
//                        movement.m_CommunityId1 = community1.Id();
//                        movement.m_CommunityId2 = community2.Id();
//                        movementFound = true;
//                    }
//                    */
//                    double aux = TestInsert( nodeId, community2, alpha, overlapp ) + TestRemove( nodeId, community1, alpha, overlapp );
//                    if( aux > bestImprovement ) {
//                        bestImprovement = aux; 
//                        movement.m_Type = E_TRANSFER;
//                        movement.m_NodeId = nodeId;
//                        movement.m_CommunityId1 = community1.Id();
//                        movement.m_CommunityId2 = community2.Id();
//                        movementFound = true;
//                    }
//                }
//            }
//        }
        return movementFound;
    }

  static void Merge( Cover& cover, long id1, long id2, double alpha, double overlapp ) {
      Cover::Community& community1 = cover.GetCommunity( id1 );
      Cover::Community& community2 = cover.GetCommunity( id2 );
      
      std::set<long> temp1 = community1.Nodes();;
      std::set<long> temp2 = community2.Nodes();;
      std::set<long> unionSet = community1.Nodes();
      std::set<long> intersectionSet = community2.Nodes();;
      for( long n : temp1 ) {
          unionSet.insert(n);
      }

      for( long n : temp2 ) {
          unionSet.insert(n);
      }
      
      double before = 0.0;
      for( long n : unionSet ) {
          if( temp1.find(n) != temp1.end() && temp2.find(n) != temp2.end() ) {
              intersectionSet.insert(n);
          }
          bool valid;
          before += Score(cover, n, alpha, overlapp, valid);
          if(!valid) {
              before=0.0;
              break;
          }
      }

      for( long n : temp2 ) {
          community2.Remove(n);
          community1.Add(n);
      }

      double after = 0.0;
      for( long n : unionSet ) {
          bool valid;
          after += Score(cover, n, alpha, overlapp, valid);
          if(!valid) {
              after=0.0;
              break;
          }
      }

      if( after < before ) {
          for( long n : temp2 ) {
              community1.Remove(n);
              community2.Add(n);
          }
          for( long n : intersectionSet ) {
              community1.Add(n);
          }

      }

  }

  static void MergeCommunities( Cover& cover, double alpha, double overlapp ) {
      const Graph& graph = cover.GetGraph();
      std::set<std::tuple<long,long> > mergeTries;
      long numNodes = graph.GetNumNodes();
      for( long i = 0; i < numNodes; ++i ) {
          const std::set<long> communities = cover.NodeCommunities( i ); 
          const long* adjacencies = graph.GetNeighbors(i);
          long degree = graph.GetDegree(i);
          for( long j = 0; j < degree; ++j ) {
              long neighbor = adjacencies[j];
              if( i < neighbor ) {
                  const std::set<long> nCommunities = cover.NodeCommunities( neighbor );
                  for( auto it = communities.begin(); it != communities.end(); ++it ) {
                      for( auto it2 = nCommunities.begin(); it2 != nCommunities.end(); ++it2 ) {
                          mergeTries.insert(std::tuple<long,long>(*it, *it2));
                      }
                  }
              }
          }
      }
      for( std::tuple<long,long> t : mergeTries ) {
          Merge( cover, std::get<0>(t),std::get<1>(t), alpha, overlapp );
      }
  }

  void RefineCommunities( Cover& cover, double alpha, double overlapp ) {
      long bestCoverSize = 0;
      long* bestCover = cover.Serialize( bestCoverSize );
      double bestScore = Score(cover, alpha, overlapp );
      double currentScore = bestScore;
      std::cout << "Initial Score: " << bestScore << std::endl;
      int lookahead = LOOKAHEAD;
      while(lookahead > 0) {
          lookahead--;
          std::cout << "Starting Iteration" << std::endl;
          //Compute new cover.
          std::vector<Movement> movements;
          for( long i = 0; i < cover.m_Graph.GetNumNodes(); ++i ) {
              Movement movement;
              if(PerformBestMovement(cover,i, alpha, overlapp, currentScore, movement)) {
                  movements.push_back(movement);
              }
          }
          std::cout << "Number of movements performed " << movements.size() << std::endl;
          Cover::Community* community1; 
          Cover::Community* community2;
          for( long i = 0; i < movements.size(); ++i ) {
              Movement& movement = movements[i];
              switch(movement.m_Type) {
                  case E_REMOVE:
                      community1 = &(cover.GetCommunity(movement.m_CommunityId1));
                      community1->Remove(movement.m_NodeId);
                      if(cover.NumCommunities(movement.m_NodeId) == 0) {
                          Cover::Community* community = new Cover::Community(cover,cover.NumCommunities()); 
                          cover.m_Communities.push_back(community);
                          community->Add(movement.m_NodeId);
                      }
                      break;
                  case E_INSERT:
                      community1 = &(cover.GetCommunity(movement.m_CommunityId1));
                      community1->Add(movement.m_NodeId);
                      break;
                  case E_TRANSFER:
                      community1 = &(cover.GetCommunity(movement.m_CommunityId1));
                      community2 = &(cover.GetCommunity(movement.m_CommunityId2));
                      community1->Remove(movement.m_NodeId);
                      community2->Add(movement.m_NodeId);
                      break;
              };
          }
          currentScore = Score(cover, alpha, overlapp );
          std::cout << "New Score " << currentScore << std::endl;
          double diff = currentScore - bestScore;
          if( diff > 0.0 ) {
              std::cout << "Cover Improved" << std::endl;
              delete [] bestCover;
              bestCover = cover.Serialize( bestCoverSize );
              bestScore = currentScore;
              lookahead = LOOKAHEAD;
              if( diff <= 0.001 ) {
                  break;
              }
          } /*else {
              // Merge communities
              cover.Deserialize(bestCover, bestCoverSize);
              MergeCommunities( cover, alpha, overlapp);
              currentScore = Score(cover, alpha, overlapp );
              std::cout << "New Score after merging communities " << currentScore << std::endl;
              double diff = currentScore - bestScore;
              if( diff > 0.0 ) {
                  std::cout << "Cover Improved" << std::endl;
                  delete [] bestCover;
                  bestCover = cover.Serialize( bestCoverSize );
                  bestScore = currentScore;
                  lookahead = LOOKAHEAD;
              }               
              
          }
          */
      }        
      cover.Deserialize(bestCover, bestCoverSize);
      delete [] bestCover;
  }
}
