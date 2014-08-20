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
#include <iostream>
#include <fstream>
#include <string.h>

#define CHECK_ARGUMENT_STRING(index, option,variable,setVariable) \
	if( strcmp(argv[index],option) == 0 ){ \
			setVariable = true; \
			if( (index+1) < argc ) { \
				variable = argv[index+1]; \
			} else { \
				printf( "Invalid options.\n" ); \
				return 1;\
			}\
		}

#define CHECK_ARGUMENT_FLOAT(index, option,variable,setVariable) \
	if( strcmp(argv[index],option) == 0 ){ \
			setVariable = true; \
			if( (index+1) < argc ) { \
				variable = atof(argv[index+1]); \
			} else { \
				printf( "Invalid options.\n" ); \
				return 1;\
			}\
		}

#define CHECK_ARGUMENT_INT(index, option,variable,setVariable) \
	if( strcmp(argv[index],option) == 0 ){ \
			setVariable = true; \
			if( (index+1) < argc ) { \
				variable = atoi(argv[index+1]); \
			} else { \
				printf( "Invalid options.\n" ); \
				return 1;\
			}\
		}

#define CHECK_FLAG(index, option,setVariable) \
	if( strcmp(argv[index],option) == 0 ){ \
			setVariable = true; \
		}

void PrintUsage() {
    std::cout << "Wrong arguments" << std::endl;
}

int main( int argc, char** argv ) {

    bool graphFileNameSet = false;
    bool outputFileNameSet = false;
    bool alphaSet = false;
    bool overlappSet = false;
    char * graphFileName = NULL;
    char * outputFileName = NULL;
    double alpha = 1.0;
    double overlapp = 0.5;

    for (uint32_t i = 1; i < argc; i++) {
        CHECK_ARGUMENT_STRING(i, "-f", graphFileName, graphFileNameSet)
        CHECK_ARGUMENT_STRING(i, "-o", outputFileName, outputFileNameSet)
        CHECK_ARGUMENT_FLOAT(i, "-a", alpha, alphaSet)
        CHECK_ARGUMENT_FLOAT(i, "-s", overlapp, overlappSet)
    }

    if( !graphFileNameSet || !outputFileNameSet ) {
        PrintUsage();
        exit(1);
    }

    triforce::Graph graph;
    graph.Load(graphFileName, 1);
    triforce::Cover* cover =  Create(&graph);
    triforce::Initialize(*cover,alpha);
    triforce::RefineCommunities(*cover, alpha, overlapp);
    std::ofstream outputFile;
    outputFile.open(outputFileName);
    Print(*cover, outputFile);
    outputFile.close();
    std::cout << "Score: " << triforce::Score( *cover, alpha, overlapp ) << std::endl;
    Destroy(cover);
    return 0;
}
