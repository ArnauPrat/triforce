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

int main( int argc, char** argv ) {
    triforce::Graph graph;
    graph.Load(argv[1], 1);
    triforce::Cover cover(graph);
    std::cout << "Score: " << triforce::Score(cover) << std::endl;
    std::cout << "Number of communities: " << cover.NumCommunities() << std::endl;
    return 0;
}
