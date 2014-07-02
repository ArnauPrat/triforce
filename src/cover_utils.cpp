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
#include <set>
#include <algorithm>


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
            if( denom == 0 ) return 0;
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

    double TestRemove( const long nodeId, const Cover::Community& community, double alpha, double overlapp ) {
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
        bool validOverlapp;
        std::set<long> nodes = community.Nodes();
        for( auto it = nodes.begin(); it != nodes.end(); ++it ) {
            before += Score( cover, *it, alpha, overlapp, validOverlapp );
        }
        comAux.Remove(nodeId);
        for( auto it = nodes.begin(); it != nodes.end(); ++it ) {
            after += Score( cover, *it, alpha, overlapp, validOverlapp );
        }
        after += Score( cover, nodeId, alpha, overlapp, validOverlapp );
        comAux.Add(nodeId);
        double increment = after - before;

        /*double pin = ((community.GetInternalEdges()-din)*2)/(double)((community.Size()-1)*(community.Size()-2));
        double increment = CheckForIncrement( community.Size() - 1, din, dout, community.GetExternalEdges()-dout+din, pin, alpha);
        */
        return  increment;
    }

    double TestInsert( const long nodeId, const Cover::Community& community, double alpha, double overlapp ) {
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
        long r = community.Size();
        bool validOverlapp;
        std::set<long> nodes = community.Nodes();
        for( auto it = nodes.begin(); it != nodes.end(); ++it ) {
            before += Score( cover, *it, alpha, overlapp, validOverlapp );
        }
        before += Score( cover, nodeId, alpha, overlapp, validOverlapp );
        comAux.Add(nodeId);
        for( auto it = nodes.begin(); it != nodes.end(); ++it ) {
            after += Score( cover, *it, alpha, overlapp, validOverlapp );
        }
        comAux.Remove(nodeId);
        double increment = after - before;
        return increment;
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
}
