from fractions import Fraction
from pprint import pprint


def matrix_inverse(mtx):
    # Return inverse of matrix --> matrix^-1
    a, b, c, d = (mtx[0][0], mtx[0][1], mtx[1][0], mtx[1][1])
    scalar = Fraction(1, 1) / ((a * d) - (b * c))
    inverter = [[d, -1*b], [-1*c, a]]
    return [[scalar * inverter[iy][ix] for ix in range(2)] for iy in range(2)]

def matrix_multiply(mtx0, mtx1):
    # Multiply mtx0 by mtx1 --> return ( mtx0 * mtx1 )
    mrun_0 = xrange(len(mtx0))
    mrun_1 = xrange(len(mtx1[0]))
    return [[sum([mtx0[spans][idy] * mtx1[idy][idx] 
        for idy in mrun_0]) for idx in mrun_1] for spans in mrun_0]

def probability_vector(n_matrix, absorber):
    # Multiply N-matrix by absorbing matrix, yank & return e[0] simplified
    prob_matrix = matrix_multiply(n_matrix, absorber)
    #print '\nbakVector:\n'; pprint(prob_matrix, width=len(prob_matrix[0]) * 20); print

    stable = [(state.numerator, state.denominator) for state in prob_matrix[0]]
    dmax = max((maxi[1] for maxi in stable))
    pvector = [(num * (dmax // den) if den != dmax else num) for num, den in stable]
    pvector.append(dmax)
    return pvector

def matrix_identity(idlen):
    # Generate identity matrices
    fracs = Fraction
    return [[(fracs(0) if x != y else fracs(1)) for x in idlen] for y in idlen]

def gauss_elim(ivm):
    _range = range
    nrows = len(ivm)
    weight = 0

    for ix in _range(nrows):
        if ivm[ix][ix] == 0:
            break
        for jx in _range(nrows):
            if jx != ix:
                weight = ivm[jx][ix] / ivm[ix][ix]
            for kx in _range(nrows * 2):
                ivm[jx][kx] = ivm[jx][kx] - weight * ivm[ix][kx]
    
    return [[ivm[i][x] / ivm[i][i] for x in _range(nrows, nrows * 2)] for i in _range(nrows)]

def gauss_redux(matrix, trange):
    qidentity = matrix_identity(trange)
    diff_matrix = [[qidentity[y][x] - matrix[y][x] for x in trange] for y in trange]
    #print '\nInverting:\n'; pprint(diff_matrix, width=(len(diff_matrix[0]) * 18)); print
    
    inverter = [diff_matrix[tx] + qidentity[tx] for tx in trange]
    resi = gauss_elim(inverter)

    #print '\nInverted:\n'; pprint(resi, width=len(resi[0]) * 18); print
    return resi

def normalize_rowspans(submatrix):
    # Normalize each row and return in fraction representations
    _fraction = Fraction
    for idx, iterows in enumerate(submatrix):
        denom = sum(iterows)
        submatrix[idx] = map(_fraction, iterows, (denom for _ in iterows))
    return submatrix

def solution(m):
    # Extract non-redundant rows; count transient and absorbing states.
    datarows = filter(any, m)
    ntransits = len(datarows)

    # Initialize Q/R iterators
    trange = xrange(ntransits)
    roffset = xrange(ntransits, ntransits + (len(m) - ntransits))

    # Normalize each rowspan into fractions of that row total.
    txr = normalize_rowspans(datarows)

    # Transient(Q) & Absorbing(R) probability matrices must be isolated
    # Canonical form given; therefore 'datarows' split into respective Q/R matrices
    Q_matrix = [[txr[nyt][nxt] for nxt in trange] for nyt in trange]
    R_matrix = [[txr[nyt][nxt] for nxt in roffset] for nyt in trange]

    # Solves for; ( identity(Q) - Q_matrix )^-1
    iq_inverse = gauss_redux(Q_matrix, trange)

    # Matrix multiplication finds stable matrix, stable[0] returned
    return probability_vector(iq_inverse, R_matrix)

if __name__ == "__main__":
    # [7, 6, 8, 21]
    sol01 = [[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], 
	    [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
    # [0, 3, 2, 9, 14]
    sol02 = [[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], 
	    [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], 
	    [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
    # []
    sol03 = [
            [0, 2, 1, 1, 0, 0, 0], 
            [0, 0, 0, 3, 4, 0, 2], 
	    [2, 1, 0, 1, 0, 1, 1], 
            [0, 0, 1, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0], 
            [0, 0, 0, 0, 0, 0, 0]
            ]

    sol04 = [
            [0, 2, 1, 4, 0, 0, 0, 0, 0, 1],
            [1, 0, 2, 1, 0, 3, 0, 0, 0, 0],
            [1, 0, 1, 1, 0, 0, 0, 2, 0, 0],
            [0, 3, 0, 0, 4, 0, 2, 0, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            #[2, 0, 3, 0, 0, 1, 1, 0, 0, 2],
            #[3, 0, 0, 2, 2, 4, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            ]

    for foo in [sol01, sol02, sol03, sol04]:
        print '\nProbability Vector: {}\n'.format(solution(foo))


