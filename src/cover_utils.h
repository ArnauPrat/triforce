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


#ifndef COVER_UTILS_H 
#define COVER_UTILS_H

#include "cover.h"
#include <ostream>

namespace triforce {

    /** @brief Computes the score of a cover.
     *  @param cover the cover to compute the score from.
     *  @return The score of the cover.*/
    double Score( const Cover& cover, double alpha, double overlapp );

    /** @brief Computes the score of a node in a cover.
     *  @param cover The cover to compute the score from.
     *  @param nodeId The id of the node to compute the score.
     *  @return The score of the node in the cover.*/
    double Score( const Cover& cover,  long nodeId, double alpha, double overlapp, bool& validOverlapp );

    /** @brief Computes if two communities are significantly similar.
     *  @param communityA The first community.  
     *  @param communityB The second community.
     *  @return true if both communities are significantly similar. false otherwise.*/
    bool Similar( const Cover& cover, const Cover::Community& communityA, const Cover::Community& communityB, double overlapp );

    /** @brief Computes a the union of all the communities a node belongs.
     *  @param cover The cover.
     *  @param nodeId The node to compute the union from.
     *  @return The set of nodes in the union.*/
    std::set< long> CommunityUnion( const Cover& cover, 
                                            const  long nodeId );
    /** @brief Prints a cover to a stream.
     *  @param cover The cover to print.
     *  @param stream The stream to print the over.*/
    void Print( const Cover& cover, std::ostream& stream );

    /** @brief Countes the number of nodes .
     *  @param cover The cover to print.
     *  @param stream The stream to print the over.
     *  @param alpha The alpha parameter controllign the cohesion of the communities*/
    void PrintZero( const Cover& cover, std::ostream& stream, double alpha, double overlapp );

    /** @brief Tests if a node should be removed from a community.
     *  @param nodeID The id of the node.
     *  @param community The community to test.
     *  @return The improvement*/
    double TestRemove( const long nodeId, const Cover::Community& comunity, double alpha, double overlapp, bool& validOverlapp );

    /** @brief Tests if a node should be inserted into a community.
     *  @param nodeID The id of the node.
     *  @param community The community to test.
     *  @return The improvement.*/
    double TestInsert( const long nodeId, const Cover::Community& comunity, double alpha, double overlapp, bool& validOverlapp );

    /** @brief Tests if a node has a valid overlapp.
     *  @param nodeId The id of the node.
     *  @return true if the has a valid overlap.*/
    bool ValidOverlapp( const Cover& cover,  long nodeId, double overlapp );

    /** @brief Tests if a node can overlap to a new community.
     *  @param nodeID The id of the node.
     *  @return true if the node can overlap.*/
    bool CanOverlapMore( const Cover& cover, const long nodeId, double overlapp );

    /** @brief Tests if inserting a node into a community violates the overlapp.
     *  @param nodeId The id of the node.
     *  @param community The community to test.
     *  @return true if inserting the node violates overlapp*/
    bool ViolatesOverlappInsert( const long nodeId, const Cover::Community& community);

    /** @brief Tests if removing a node into a community violates the overlapp.
     *  @param nodeId The id of the node.
     *  @param community The community to test.
     *  @return true if removing the node violates overlapp*/
    bool ViolatesOverlappRemove( const long nodeId, const Cover::Community& community);

    /** @brief Initializes a cover with an initial partition
     *  @param alpha The alpha parameter to use.
     *  @param overlapp The maximum alowed overlapp.*/
    void Initialize( Cover& cover, double alpha, double overlapp );

    /** @brief Refines the communities in the cover.
     *  @param alpha The alpha value controling the cohesivness of the communities
     *  @param overlapp The level of overlapp allowed.*/
    void RefineCommunities( Cover& cover, double alpha , double overlapp );

    bool PerformBestMovement( Cover& cover, const long nodeId, double alpha, double overlapp, double bestScore, Movement& movement );

}

#endif
