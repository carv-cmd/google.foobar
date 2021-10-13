from pprint import pprint
from fractions import Fraction


def matrix_multiply(mtx0, mtx1):
    # Multiply mtx0 by mtx1 --> return ( mtx0 * mtx1 )
    mrows = xrange(len(mtx0))
    mcols = xrange(len(mtx1[0]))
    return [[sum([mtx0[spans][idy] * mtx1[idy][idx] 
        for idy in mrows]) 
        for idx in mcols] 
        for spans in mrows]

def matrix_slice(m2slice):
    sliced = [(state.numerator, state.denominator) for state in m2slice]
    dmax = int(max((maxi[1] for maxi in sliced)))
    pvector = [int((num * (dmax // den) if den != dmax else num)) for num, den in sliced]
    pvector.append(dmax)
    return pvector

def gauss_redux(redux):
    _range = xrange
    nrows = len(redux)
    trange = _range(nrows)

    for ix in trange:
        weight = 0
        if redux[ix][ix] == 0:
            break
        for jx in trange:
            if jx != ix:
                weight = redux[jx][ix] / redux[ix][ix]
            for kx in _range(nrows * 2):
                redux[jx][kx] = redux[jx][kx] - weight * redux[ix][kx]
    
    return [[redux[i][x] / redux[i][i] for x in _range(nrows, nrows * 2)] for i in trange]

def normalize_rowspans(submatrix):
    # Normalize each row and return in fraction representations
    _fraction = Fraction
    for idx, iterows in enumerate(submatrix):
        denom = sum(iterows)
        submatrix[idx] = map(_fraction, iterows, (denom for _ in iterows))
    return submatrix

def main_loop(txr):
    #qidentity = matrix_identity(trange)
    fracs = Fraction

    # Transient(Q) & Absorbing(R) probability matrices must be isolated
    trange = xrange(len(txr))
    offset = xrange(len(txr), len(txr[0]))
    Q_matrix = [[txr[nyt][nxt] for nxt in trange] for nyt in trange]
    R_matrix = [[txr[nyt][nxt] for nxt in offset] for nyt in trange]
    
    # Generate identity matrix and find ( Ir - Q ) to invert
    qidentity = [[(fracs(0) if x != y else fracs(1)) for x in trange] for y in trange]
    qdiff = [[qidentity[y][x] - Q_matrix[y][x] for x in trange] for y in trange]
    Q_matrix = [qdiff[tx] + qidentity[tx] for tx in trange]

    # Solve for; ( ( identity(Q) - Q_matrix )^-1 * R )
    # If answers werent assumed solvable; use LD or Cholesky
    iq_inverse = gauss_redux(Q_matrix)
    multi = matrix_multiply(iq_inverse, R_matrix)[0]

    # Slice and return relevant probabilities
    return matrix_slice(multi)

def solution(m):
    # Sort and extract non-redundant rows; count transient and absorbing states.
    t, a = [], []
    for idx in range(len(m)):
        if any(m[idx]):
            t.append(idx)
        else:
            a.append(idx)
    canonical = [[m[xrow][xcol] for xcol in t + a] for xrow in t]

    # 1/1 chance of getting abosrbed by one possible state
    if len(m) - len(canonical) == 1:
        oh_markov = [1, 1]
    # Calculate absorbing states of given markov chain 
    else:
        # Normalize each rowspan into fractions of that row total.
        canonical = normalize_rowspans(canonical)
        oh_markov = main_loop(canonical)
    return oh_markov

def test_cases():
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
            }
    return test_dicks

if __name__ == "__main__":
    #shimy = test_cases()['09']
    #shimy = {'09': shimy}
    shimy = test_cases()
    for dicks in sorted(shimy.keys()):
        ball, sak = shimy[dicks][0], shimy[dicks][1]
        print '\nTestSet({}):\n\tExpected: {}\n\tReturned: {}'.format(
                dicks, sak, solution(ball))



