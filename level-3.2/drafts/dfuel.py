from fractions import Fraction
from pprint import pprint


def scalar_multiply(scalar, matrix):
    # Multiply a scalar across all values of 'matrix'
    mlen = len(matrix)
    return [[scalar * matrix[iy][ix] for ix in range(mlen)] for iy in range(mlen)]

def matrix_multiply(mtx0, mtx1):
    # Multiply mtx0 by mtx1 --> return ( mtx0 * mtx1 )
    mx0 = len(mtx0)
    mx1 = len(mtx1[0])
    return [[sum([mtx0[px][iy] * mtx1[iy][ix] for iy in range(mx0)]) 
        for ix in range(mx1)] for px in range(mx0)]

def matrix_subtract(mtx0, mtx1):
    # Subtract mtx1 from mtx0 --> return ( mtx0 - mtx1 )
    mlen = len(mtx0)
    return [[mtx1[y][x] - mtx0[y][x] for x in range(mlen)] for y in range(mlen)]

def matrix_inverse(mtx):
    # Return inverse of matrix --> matrix^-1
    a, b, c, d = (mtx[0][0], mtx[0][1], mtx[1][0], mtx[1][1])
    foos = Fraction(1, 1) / ((a * d) - (b * c))
    return scalar_multiply(foos, [[d, -1*b], [-1*c, a]])

def matrix_identity(alen):
    # Generate identity matrices
    fracs = Fraction
    return [[(fracs(0) if x != y else fracs(1)) for x in range(alen)] for y in range(alen)]

def extract_data(n_matrix, absorber):
    # Multiply fund_inv by absorbing matrix, yank & return e[0] normalized
    prob_matrix = matrix_multiply(n_matrix, absorber)[0]
    dmax = max((i.denominator for i in prob_matrix))
    pvector = [(frak.numerator * (dmax // frak.denominator) 
        if frak.numerator != dmax else frak.numerator) 
        for frak in prob_matrix]
    pvector.append(dmax)
    return pvector

def fundamental_inverse(numt, tmatrix):
    # Find  ( I - M_transient ) ^ -1
    m_matrix = matrix_subtract(tmatrix, matrix_identity(numt))
    return matrix_inverse(m_matrix)

def isolate_qr(tnum, txr):
    # From the ( Q x R ) matrix, isolate Q(transient) and R(absorbing) separately
    anum = len(txr[0]) - tnum
    transit = [[txr[nyt][nxt] for nxt in range(tnum)] for nyt in range(tnum)]
    absorb = [[txr[nya][nxa] for nxa in range(tnum, tnum + anum)] for nya in range(tnum)]
    return transit, absorb

def normalize_rowspans(submatrix):
    # Normalize each row and return in fraction representations
    _sum, _map, _fraction = sum, map, Fraction
    for idx, iterows in enumerate(submatrix):
        denom = _sum(iterows)
        submatrix[idx] = _map(_fraction, iterows, (denom for _ in iterows))
    return submatrix

def solution(m):
    # Extract non-redundant rows; count transient and absorbing states.
    datarows = filter(any, m)
    ntransits = len(datarows)
    nabsorbs = len(m) - ntransits

    # Normalize each rowspan into fractions of that row total.
    # 'fractions' module used through final vector extraction.
    txr = normalize_rowspans(datarows)

    # Transient(Q) & Absorbing(R) probability matrices must be isolated
    # Canonical form given, 'datarows' split into respective matrices
    transient = [[txr[nyt][nxt] for nxt in range(ntransits)] 
            for nyt in range(ntransits)]

    absorbing = [[txr[nya][nxa] for nxa in range(ntransits, ntransits + nabsorbs)] 
            for nya in range(ntransits)]

    # 
    pbmatrix = fundamental_inverse(num_transient, transient)

    # 
    return extract_data(pbmatrix, absorbing)

def matrix_states(matrix):
    # Identify transient and absorbing states in given matrix
    usingrows = filter(any, matrix)
    num_transients = len(usingrows)
    num_absorbing = len(matrix) - num_transients
    usingrows = normalize_rowspans(usingrows)
    return num_absorbing, num_transients, usingrows

def main_loop(tests):
    for testing in tests:
        print '\nTestSet:'
        pprint(testing, width=40)
        nabsorb, ntransit, using = matrix_states(testing)
        transient, absorbing = isolate_qr(nabsorb, ntransit, using)
        fundmatrix = fundamental_inverse(ntransit, transient)
        print "\nAbsorbing Probabilities:\n\t{}".format(
                extract_data(fundmatrix, absorbing))


if __name__ == "__main__":
    # [7, 6, 8, 21]
    sol01 = [[0, 2, 1, 0, 0], [0, 0, 0, 3, 4], 
	    [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
    # [0, 3, 2, 9, 14]
    sol02 = [[0, 1, 0, 0, 0, 1], [4, 0, 0, 3, 2, 0], 
	    [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], 
	    [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]]
    #main_loop()
    print solution(sol02)


