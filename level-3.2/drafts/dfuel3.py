from fractions import Fraction
from pprint import pprint

def gauss_redux(matrix, trange):
    # Gauss-Jordan elimination
    _range = range
    nrows = len(matrix)
    qidentity = matrix_identity(trange)
    diff_matrix = [[qidentity[y][x] - matrix[y][x] for x in trange] for y in trange]
    inv = [diff_matrix[tx] + qidentity[tx] for tx in trange]
    weight = 0
    for ix in _range(nrows):
        if inv[ix][ix] == 0:
            break
        for jx in _range(nrows):
            if jx != ix:
                weight = inv[jx][ix] / inv[ix][ix]
            for kx in _range(nrows * 2):
                inv[jx][kx] = inv[jx][kx] - weight * inv[ix][kx]
        weight = 0
    return [[inv[i][x] / inv[i][i] for x in _range(nrows, nrows * 2)] for i in _range(nrows)]

def matrix_multiply(mtx0, mtx1):
    # Multiply mtx0 by mtx1 --> return ( mtx0 * mtx1 )
    mrun_0 = xrange(len(mtx0))
    mrun_1 = xrange(len(mtx1[0]))
    return [[sum([mtx0[spans][idy] * mtx1[idy][idx] 
        for idy in mrun_0]) for idx in mrun_1] for spans in mrun_0]

def probability_vector(n_matrix, absorber):
    # Multiply N-matrix by absorbing matrix, yank & return e[0] simplified
    prob_matrix = matrix_multiply(n_matrix, absorber)[0]
    stable = [(state.numerator, state.denominator) for state in prob_matrix]
    dmax = max((maxi[1] for maxi in stable))
    pvector = [(num * (dmax // den) if den != dmax else num) for num, den in stable]
    pvector.append(dmax)
    return pvector

def matrix_identity(idlen):
    # Generate identity matrices
    fracs = Fraction
    return [[(fracs(0) if x != y else fracs(1)) for x in idlen] for y in idlen]

def normalize_rowspans(submatrix):
    # Normalize each row and return in fraction representations
    _fraction = Fraction
    for idx, iterows in enumerate(submatrix):
        denom = sum(iterows)
        submatrix[idx] = map(_fraction, iterows, (denom for _ in iterows))
    return submatrix

def lud_multiply(lower, upper):
    # Multiply mtx0 by mtx1 --> return ( mtx0 * mtx1 )
    #return [
    #[sum([mtx0[spans][idy] * mtx1[idy][idx] 
    #       for idy in mrun_0]) 
    #           for idx in mrun_1] 
    #               for spans in mrun_0]
    xrun = xrange(len(lower[0]))
    for 
    foos = [[sum((lower[slot][iy] * upper[iy][ix] 
        for iy in xrun)) 
        for ix in xrun] 
        for slot in xrun]
    return foos

def ludcmp(m, trange):
    # EQ(L * U = A)
    #   Decompose matrix into upper and lower triangluar matrices
    _range = range
    fracs = Fraction
    nrows = len(m)
    qidentity = matrix_identity(trange)

    diff = [[qidentity[y][x] - m[y][x] for x in trange] for y in trange]
    lower = [[(diff[x][y] if y < x + 1 else fracs(0)) for y in trange] for x in trange]
    upper = [[(diff[x][y] if y >= x else fracs(0)) for y in trange] for x in trange]
    multi = lud_multiply(lower, upper)

    displays =  (('Upper', upper), ('Initial', diff), ('Lower', lower), ('multi', multi))
    for t, itr in displays:
        fiz = [['{}/{}'.format(i.numerator, i.denominator) for i in item] for item in itr]
        print '\n{} Matrix'.format(t); pprint(fiz, width=len(fiz[0] * 15))
    #print '\nUpper * Lower'; pprint(multi, width=120)

    # EQ( 2.3.8 );
    #   [ i < j ]: (a[i][1]*B[1][j]) + (a[i][2]*B[2][j]) +,...,+ ( a[i][i] * B[i][j] ) = a[i][j]
    # EQ( 2.3.9 );
    #   [ i = j ]: (a[i][1]*B[1][j]) + (a[i][2]*B[2][j]) +,...,+ ( a[i][i] * B[j][j] ) = a[i][j]
    # EQ( 2.3.10 );
    #   [ i > j ]: (a[i][1]*B[1][j]) + (a[i][2]*B[2][j]) +,...,+ ( a[i][j] * B[j][j] ) = a[i][j]
    # EQ( 2.3.11 );
    #   a[i][i] === 1   i = 1, ..., N
    #inv = [diff_matrix[tx] + qidentity[tx] for tx in trange]
    #print matrix

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
    #iq_inverse = gauss_redux(Q_matrix, trange)
    iq_inverse = ludcmp(Q_matrix, trange)
    if True: return None

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
    # LU-decomp test
    sol03 = [
            [3, 0, 5, 2, 4],
            [3, 0, 2, 4, 0],
            [2, 3, 0, 2, 7],
            [0, 2, 2, 0, 4],
            [0, 3, 5, 7, 0],
            ]
    sol04 = [
            [3, 0, 5],
            [3, 0, 2],
            [2, 3, 0],
            ]
    ludcmp(sol04, xrange(len(sol04)))

    #print solution(sol01)
    #print solution(sol02)


