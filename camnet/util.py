# Copyright (c) 2013 Riccardo Lucchese, riccardo.lucchese at gmail.com
#
# This software is provided 'as-is', without any express or implied
# warranty. In no event will the authors be held liable for any damages
# arising from the use of this software.
#
# Permission is granted to anyone to use this software for any purpose,
# including commercial applications, and to alter it and redistribute it
# freely, subject to the following restrictions:
#
#    1. The origin of this software must not be misrepresented; you must not
#    claim that you wrote the original software. If you use this software
#    in a product, an acknowledgment in the product documentation would be
#    appreciated but is not required.
#
#    2. Altered source versions must be plainly marked as such, and must not be
#    misrepresented as being the original software.
#
#    3. This notice may not be removed or altered from any source
#    distribution.

import numpy
import networkx as nx


def numpy_matrix_grow(mat, nr_rows=1, nr_cols=1):
    # grow the matrix 'mat' by appending rows and columns of zeros
    assert isinstance(mat, numpy.ndarray)
    assert isinstance(nr_rows, int)
    assert isinstance(nr_cols, int)
    assert nr_rows >= 0
    assert nr_cols >= 0

    # extend the transition matrix with nr_rows rows and nr_cols columns
    rows, cols = mat.shape
    newmat = numpy.zeros((rows+nr_rows,cols+nr_cols))
    newmat[:rows,:cols] = mat
    return newmat


def numpy_matrix_row_str(mat):
    _str = "["
    for e in mat:
        _str += "%.2f " % e
                
    return _str[:-1] + ']'

