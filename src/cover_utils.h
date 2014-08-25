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

    void Initialize( Cover& cover, double alpha );

    /** @brief Computes the score of a cover.
     *  @param cover the cover to compute the score from.
     *  @return The score of the cover.*/
    double Score( const Cover& cover, double alpha, double overlapp );

    /** @brief Computes the score of a node in a cover.
     *  @param cover The cover to compute the score from.
     *  @param nodeId The id of the node to compute the score.
     *  @param alpha The alpha coefficient.
     *  @param overlap The allowed overlap.
     *  @return The score of the node in the cover.*/
    double Score( const Cover& cover,  long nodeId, double alpha, double overlap);


    /** @brief Prints a cover to a stream.
     *  @param cover The cover to print.
     *  @param stream The stream to print the over.*/
    void Print( const Cover& cover, std::ostream& stream );

    /** @brief Tests if a node should be removed from a community.
     *  @param cover The cover.
     *  @param nodeId The id of the node.
     *  @param community The community to test.
     *  @param alpha The alpha coefficient.
     *  @param overlap The allowed overlap.
     *  @return The improvement*/
    double TestRemove( const Cover& cover, const long nodeId, const long community, double alpha, double overlapp);

    /** @brief Tests if a node should be removed from a community.
     *  @param cover The cover.
     *  @param nodeId The id of the node.
     *  @param community The community to test.
     *  @param alpha The alpha coefficient.
     *  @param overlap The allowed overlap.
     *  @return The improvement*/
    double TestInsert( const Cover& cover, const long nodeId, const long comunity, double alpha, double overlapp);

    /** @brief Refines the communities in the cover.
     *  @param alpha The alpha value controling the cohesivness of the communities
     *  @param overlapp The level of overlapp allowed.*/
    void RefineCommunities( Cover& cover, double alpha , double overlapp );
}

#endif
