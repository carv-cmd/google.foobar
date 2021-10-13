from argparse import ArgumentParser
from sys import stdout

from collections import defaultdict
from operator import xor, mod
from array import array
from math import log as log2
from math import pow as expo

###########################################
def parse_args():
    _parser = ArgumentParser(description='Parse cryptos')
    _parser.add_argument('-f')
    _parser.add_argument('-s')
    _parser.add_argument('-o')
    _parsed = _parser.parse_args()
    return _parsed.f, _parsed.s, _parsed.o

###########################################
def solar_queue01(start, offset):
    stdw, stdf = stdout.write, stdout.flush
    rcvector = []
    for _col in range(start, start + offset):
        for _row in range((start + offset) - _col):
            stdw('\033[2K>>> Col {} Row: {}\r'.format(_col, _row))
            rcvector.append(_row * offset + _col)
            stdf()
    print '\033[2K'
    return reduce(lambda x, y: x ^ y, rcvector)

###########################################
def solar_queue02(start, length):
    stdw, stdf = stdout.write, stdout.flush
    _xrange, _reduce, _xor = xrange, reduce, xor
    xors = 0
    for _row in range(length):
        stdout.write('\r\033[2K>>> RowSpan: {} / {}'.format(_row, length - 1))
        xors ^= _reduce(_xor, _xrange(start, start + length - _row))
        start += length
        stdout.flush()
    print '\033[2K'
    return xors

###########################################
def solar_queue03(start, length):
    stdw, stdf = stdout.write, stdout.flush
    _xrange, _reduce, _enum, _xor = xrange, reduce, enumerate, xor
    foobaz = array('L')
    foobaz.extend([_reduce(_xor, _xrange(start + (ix * length), start + (ix * length) + length - ix)) 
            for ix, ic in _enum(_xrange(start, start + length))])
    print 'Performing Final XOR Iteration'
    redux = reduce(_xor, foobaz)
    print '\033[2K'
    return redux

###########################################
def solar_queue04(start, length):
    stdw, stdf = stdout.write, stdout.flush
    _xrange, _reduce, _enum, _xor = xrange, reduce, enumerate, xor
    _array = array('L')
    _arrayappend = _array.append

    for _ex, _ix in _enum(_xrange(start, start + length)):
        stdw('\r\033[2K>>> ColSpan: {} / {}'.format(_ex, length - 1))
        _arrayappend(_reduce(_xor, _xrange(_ix, _ix + (length * (length - _ex)), length)))
        stdf()

    print '\n>>> Performing Final XOR Iteration\n'
    print '\033[2K'
    return reduce(_xor, _array)

###########################################
def solar_queue05(start, length):
    _xrange = xrange
    _reduce = reduce
    _enum = enumerate
    _xor = xor
    return reduce(_xor, [_reduce(_xor, _xrange(_ix, _ix + (length * (length - _ex)), length)) 
        for _ex, _ix in _enum(_xrange(start, start + length))])

###########################################
def aux_xor(sxs, exe):
    # aux_xor operations
    # if head odd; decrement count and pop head
    # if tail even; decrement count and pop tail
    # if remaining list divisible by 2, but not 4; add odd parity

    pass

def solar_queue06(start, length):
    #_xrange _reduce, _maxi, _xor = reduce, enumerate, max, xor
    chksum = 0
    for _row in xrange(length):
        idx01 = start + _row * length
        idx02 = length - _row
        print '\nStart/Length({}/{}): row: {}'.format(start, length, _row)
        print idx01, idx02
        # [_reduce(_xor, _xrange(_ix, _ix + (length * (length - _ex)), length)) 
        pass
    
    #listing = array('L')
    #appender = listing.append
    #for i in xrange(start, start + length):
    #    appender()
    #return reduce(_xor, [_reduce(_xor, _xrange(_ix, _ix + (length * (length - _ex)), length)) 
    #    for _ex, _ix in _enum(_xrange(start, start + length))])

###########################################
def _test_cases(func, strt, leng):
    # -> start[ID-Range][0]  -->  [ 0 <= r <= 2000000000 ]
    # -> length[ilength][1]  -->  [ 1 < x <= 44270 ]
    # Verify-Even: [ (e|o),(e|e),(e|1),(Z|e),(e|M) || (o|e),(o|o),(Z|o),(o|1),(o|M) ]
    gkcopy = globals().copy()
    gcopy = list(filter(lambda f: func in f, gkcopy.keys()))[-1]
    funkshin = gkcopy[gcopy]
    print '\nExecuting Function: {}\n'.format(gcopy)
    if int(leng) != 0:
        tests = {'UserIn': (int(strt), int(leng))}
    else: 
        tests = {'Zero/Even': (0, 4), 'Zero/Odd': (0, 3),
                'Even/One': (6, 1), 'Odd/One': (17, 1),
                'Even/Even': (26, 4), 'Even/Odd': (4, 7),
                'Odd/Odd': (3, 3), 'Odd/Even': (17, 4),
                'Zero/One:': (0, 1), 'One/One': (100, 1)}
    for key in sorted(tests.keys()):
        _st, _en = tests[key]
        squeue =  funkshin(_st, _en)
        print '>>> Start/Length: [ {} / {} ] --> ( {} : {} )'.format(
                _st, _en, key, squeue)


if __name__ == "__main__":
    _test_cases(*parse_args())

















































