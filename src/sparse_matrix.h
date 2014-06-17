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


#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include <map>
namespace triforce {

    template <class T>
    class SparseMatrix {
        public:
            SparseMatrix() {
            }
            ~SparseMatrix() {
            }

            /** @brief Sets the values at a given position.
             *  @param i The ith row of the matrix. 
             *  @param j The jth column of the matrix.
             *  @param val The value to set. */ 
            void Set( const long i, const long j, const T& val ) {
                Position p;
                if( i < j ) {
                    p.m_I = i;
                    p.m_J = j;
                } else {
                    p.m_J = i;
                    p.m_I = j;
                }
                m_Matrix[p] = val;
            }

            /** @brief Gets the values of a given position.
             *  @param i The ith row of the matrix. 
             *  @param j The jth column of the matrix. 
             *  @return The current value.*/ 
            T Get( const long i, const long j ) {
                Position p;
                if( i < j ) {
                    p.m_I = i;
                    p.m_J = j;
                } else {
                    p.m_J = i;
                    p.m_I = j;
                }
                return m_Matrix[p];
            }

            /** @brief Increments the value at a given position.
             *  @param i The ith row of the matrix. 
             *  @param j The jth column of the matrix. */ 
            void Inc( const long i, const long j) {
                Position p;
                if( i < j ) {
                    p.m_I = i;
                    p.m_J = j;
                } else {
                    p.m_J = i;
                    p.m_I = j;
                }
                auto it = m_Matrix.find(p);
                if( it != m_Matrix.end() ) {
                    (*it).second++;
                } else {
                    m_Matrix[p]=0;
                }
            }

            /** @brief Decrements the value at a given position.
             *  @param i The ith row of the matrix. 
             *  @param j The jth column of the matrix. */ 
            void Dec( const long i, const long j) {
                Position p;
                if( i < j ) {
                    p.m_I = i;
                    p.m_J = j;
                } else {
                    p.m_J = i;
                    p.m_I = j;
                }
                m_Matrix[p]--;
            }


            void Clear() {
                m_Matrix.clear();
            }


        private:
            struct Position {
                long m_I;   /**< @brief The row of the position.*/
                long m_J;   /**< @brief The column of the position.*/

                bool operator<( const Position& b ) const {
                    if( m_I < b.m_I ) return true;
                    if( m_I > b.m_I ) return false;
                    if( m_J < b.m_J ) return true;
                    return false;
                }
            };

            std::map<Position, T>   m_Matrix;   /**< @brief A map storing the matrix.*/
    };
}

#endif

