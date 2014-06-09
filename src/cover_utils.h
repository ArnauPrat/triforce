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

namespace triforce {


    /** @brief Computes the score of a cover.
     *  @param cover the cover to compute the score from.
     *  @return The score of the cover.*/
    double Score( const Cover& cover );

    /** @brief Computes the score of a node in a cover.
     *  @param cover The cover to compute the score from.
     *  @param nodeId The id of the node to compute the score.
     *  @return The score of the node in the cover.*/
    double Score( const Cover& cover,  long nodeId );

    /** @brief Computes if two communities are significantly similar.
     *  @param communityA The first community.  
     *  @param communityB The second community.
     *  @return true if both communities are significantly similar. false otherwise.*/
    bool Similar( const Cover::Community& communityA, const Cover::Community& communityB );

    /** @brief Computes a the union of all the communities a node belongs.
     *  @param cover The cover.
     *  @param nodeId The node to compute the union from.
     *  @return The set of nodes in the union.*/
    std::set< long> CommunityUnion( const Cover& cover, 
                                            const  long nodeId );
}

#endif
