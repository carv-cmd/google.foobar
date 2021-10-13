from pprint import pprint
from fractions import Fraction


def probability_slice(m2slice):
    sliced = [(state.numerator, state.denominator) for state in m2slice]
    dmax = int(max((maxi[1] for maxi in sliced)))
    pvector = [int((num * (dmax // den) if den != dmax else num)) for num, den in sliced]
    pvector.append(dmax)
    return pvector

def matrix_multiply(mtx0, mtx1):
    # Multiply mtx0 by mtx1 --> return ( mtx0 * mtx1 )
    mrows = xrange(len(mtx0))
    mcols = xrange(len(mtx1[0]))
    return [[sum([mtx0[spans][idy] * mtx1[idy][idx] 
        for idy in mrows]) 
        for idx in mcols] 
        for spans in mrows]

def gauss_redux(matrix):
    _range = range
    nrows = len(matrix)
    trange = _range(nrows)

    for ix in trange:
        weight = 0
        if matrix[ix][ix] == 0:
            break
        for jx in trange:
            if jx != ix:
                weight = matrix[jx][ix] / matrix[ix][ix]
            for kx in _range(nrows * 2):
                matrix[jx][kx] = matrix[jx][kx] - weight * matrix[ix][kx]

    baksub = []
    bakpend = baksub.append
    for ir in trange:
        bakpend([matrix[ir][x] / matrix[ir][ir] for x in _range(nrows, nrows * 2)])
    return baksub

def normalize_rowspans(submatrix):
    # Normalize each row and return in fraction representations
    _fraction = Fraction
    for idx, iterows in enumerate(submatrix):
        denom = sum(iterows)
        submatrix[idx] = map(_fraction, iterows, (denom for _ in iterows))
    return submatrix

def main_loop(matrix):
    # Initialize Q/R iterators & find absorbing(R) starting index offset.
    # Transient(Q) & Absorbing(R) probability matrices must be isolated.
    ntransits = len(matrix)
    trange = range(ntransits)
    aoffset = xrange(ntransits, len(matrix[0]))
    Q_matrix = [[matrix[nyt][nxt] for nxt in trange] for nyt in trange]
    R_matrix = [[matrix[nyt][nxt] for nxt in aoffset] for nyt in trange]
    
    # Solve for; ( identity(Q) - Q_matrix ) and extend rows with corresponding identity .
    fracs = Fraction
    qidentity = [[(fracs(0) if x != y else fracs(1)) for x in trange] for y in trange]
    diff_matrix = [[qidentity[y][x] - matrix[y][x] for x in trange] for y in trange]
    iqdiff = [diff_matrix[tx] + qidentity[tx] for tx in trange]

    # Solve for; ( ( identity(Q) - Q_matrix )^-1 * R )
    Q_inverse = gauss_redux(iqdiff)  #, qidentity, R_matrix)
    multi = matrix_multiply(Q_inverse, R_matrix)[0]

    # Slice and return relevant probabilities
    return probability_slice(multi)

def solution(m):
    # Extract non-redundant rows; count transient and absorbing states.
    t, a = [], []
    for index, rowspan in enumerate(m):
        if any(rowspan):
            t.append(index)
        else:
            a.append(index)
    canonical = [[m[xrow][xcol] for xcol in t + a] for xrow in t]

    if len(m) -len(t) == 1:
        oh_markov = [1, 1]
    else:
        # Normalize each rowspan into fractions of that row total.
        #oh_markov = main_loop(datarows, trange, roffset)
        canonical = normalize_rowspans(canonical)
        oh_markov = main_loop(canonical)

    return oh_markov


def test_cases():
    #fuck = [['{}/{}'.format(dy.numerator, dy.denominator) for dy in dx] for dx in Q_matrix]
    #print '\n'; pprint(fuck, width=80)
    test_dicks = {
            '00': ([
                [0, 2, 1, 0, 0], 
                [0, 0, 0, 3, 4], 
                [0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0]], 
                [7, 6, 8, 21]),
            '01': ([
                [0, 1, 0, 0, 0, 1], 
                [4, 0, 0, 3, 2, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0]],
                [0, 3, 2, 9, 14]),
            '02': ([
                [0, 2, 1, 0, 0], 
                [0, 0, 0, 3, 4], 
                [0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0]], 
                [7, 6, 8, 21]),
            '03': ([
                [0, 1, 0, 0, 0, 1], 
                [4, 0, 0, 3, 2, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0]], 
                [0, 3, 2, 9, 14]),
            '04': ([
                [1, 2, 3, 0, 0, 0], 
                [4, 5, 6, 0, 0, 0], 
                [7, 8, 9, 1, 0, 0], 
                [0, 0, 0, 0, 1, 2], 
                [0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0]], 
                [1, 2, 3]),
            '05': ([[0]], [1, 1]),
            '06': ([
                [0, 0, 12, 0, 15, 0, 0, 0, 1, 8], 
                [0, 0, 60, 0, 0, 7, 13, 0, 0, 0], 
                [0, 15, 0, 8, 7, 0, 0, 1, 9, 0], 
                [23, 0, 0, 0, 0, 1, 0, 0, 0, 0], 
                [37, 35, 0, 0, 0, 0, 3, 21, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [1, 2, 3, 4, 5, 15]),
            '07': ([
                [ 0,  7,  0, 17,  0,  1,  0,  5,  0,  2], 
                [ 0,  0, 29,  0, 28,  0,  3,  0, 16,  0], 
                [ 0,  3,  0,  0,  0,  1,  0,  0,  0,  0], 
                [48,  0,  3,  0,  0,  0, 17,  0,  0,  0], 
                [ 0,  6,  0,  0,  0,  1,  0,  0,  0,  0], 
                [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], 
                [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], 
                [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], 
                [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0], 
                [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0]],
                [4, 5, 5, 4, 2, 20]),
            '08':([
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [1, 1, 1, 1, 1, 5]),
            '09': ([
                [1, 1, 1, 0, 1, 0, 1, 0, 1, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 0, 1, 1, 1, 0, 1, 0, 1, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 0, 1, 0, 1, 1, 1, 0, 1, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 0, 1, 0, 1, 0, 1, 1, 1, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [1, 0, 1, 0, 1, 0, 1, 0, 1, 1], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [2, 1, 1, 1, 1, 6]),
            '10': ([
                [0, 86, 61, 189, 0, 18, 12, 33, 66, 39], 
                [0, 0, 2, 0, 0, 1, 0, 0, 0, 0], 
                [15, 187, 0, 0, 18, 23, 0, 0, 0, 0], 
                [1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [6, 44, 4, 11, 22, 13, 100]),
            '11': ([
                [0, 0, 0, 0, 3, 5, 0, 0, 0, 2], 
                [0, 0, 4, 0, 0, 0, 1, 0, 0, 0], 
                [0, 0, 0, 4, 4, 0, 0, 0, 1, 1], 
                [13, 0, 0, 0, 0, 0, 2, 0, 0, 0], 
                [0, 1, 8, 7, 0, 0, 0, 1, 3, 0], 
                [1, 7, 0, 0, 0, 0, 0, 2, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]],
                [1, 1, 1, 2, 5]),
            }
    return test_dicks


if __name__ == "__main__":

    #shimy = test_cases()['09']
    #shimy = {'09': shimy}
    shimy = test_cases()
    for dicks in sorted(shimy.keys()):
        ball, sak = shimy[dicks][0], shimy[dicks][1]
        fok = solution(ball)
        print 'Test({}):\n\tExpect: {}\n\tActual: {}'.format(dicks, sak, fok)

    #print '\n>> MARKOV:'; pprint(s, width=len(s) * 8)
    #fok = solution(s)
    #print '\n>>> RETURNED: {}'.format(fok)
    #print '\n>>>> returned: type({})'.format(type(fok[0]))


