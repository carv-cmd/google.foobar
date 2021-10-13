from argparse import ArgumentParser
from sys import stdout
#from dis import dis


###########################################
def parse_args():
    _parser = ArgumentParser(description='Parse cryptos')
    _parser.add_argument('-f')
    _parser.add_argument('-s')
    _parser.add_argument('-o')
    _parsed = _parser.parse_args()
    return _parsed.f, _parsed.s, _parsed.o


###########################################
def solar_queue00(start, length):
    _xrange, _reduce, _enum, _xor = xrange, reduce, enumerate, xor
    _array = array('L')
    _arrayappend = _array.append
    for _ex, _ix in _enum(_xrange(start, start + length)):
        _arrayappend(
                _reduce(_xor, 
                    _xrange(_ix, _ix + (length * (length - _ex)), length)))
    return reduce(_xor, _array)


###########################################
def xor_rowspan(start, length, _rxor=0):
    if not length:
        _rxor = 0
    elif length == 1:
        _rxor = start
    elif length == 2:
        _rxor = start ^ start + 1
    else:
        if start & 1:
            _rxor ^= start
            length -= 1
            start += 1
        if length & 1:
            length -= 1
            _rxor ^= start + (length)
        _rxor ^= (1 if length % 4 else 0)
    return _rxor

def solar_queue01(start, length):
    _xorbits = xor
    _rowxor = xor_rowspan
    _chksum = (_rowxor(start + _row * length, length - _row) 
            for _row in xrange(length))
    return reduce(_xorbits, _chksum)


###########################################
def solar_queue02(start, length):
    chksum = 0
    for irow in xrange(length):
        istart = start + irow * length
        ilength = length - irow
        print '\nProcessing: {}'.format(range(istart, istart + ilength))

        if ilength < 3:
            chksum ^= istart
            print '\tilength({}) < 3: chksum={}'.format(ilength, chksum)

            if ilength & 2: 
                chksum ^= istart + 1
                print '\tilength({}) == 2: chksum={}'.format(ilength, chksum)

        else:
            print '\n\tL > 2: (ilength: {}, istart: {}, chksum: {})\n'.format(
                    ilength, istart, chksum)

            if istart & 1:
                ilength -= 1; istart += 1
                chksum ^= istart
                print '\t\tiS: odd: (ilength: {}, istart: {}, chksum: {})'.format(
                        ilength, istart, chksum)

            if ilength & 1:
                ilength -= 1
                chksum ^= istart + ilength
                print '\t\tiL: odd: (ilength: {}, istart: {}, chksum: {})'.format(
                        ilength, istart, chksum)

            if ilength % 4:
                chksum ^= 1
                print '\t\tlen % 4: (parity:{} : chksum:{})'.format(
                        ilength % 4, chksum)

    return chksum


###########################################
def solar_queue03(start, length):
    chksum = 0
    for irow in xrange(length):
        istart = start + irow * length
        ilength = length - irow

        _astart = [istart, 1, istart + 1, 0]
        _alength = [ilength, 1, ilength + 1, 0]

        smod = _astart[istart % 4]
        lmod = _alength[ilength % 4]

        print smod, lmod

        chksum ^= lmod ^ smod

    return chksum

###########################################
def _test_cases(func, start, length):
    gkcopy = globals().copy()
    gcopy = list(filter(lambda f: func in f, gkcopy.keys()))[-1]
    print '\nExecuting Function: {} w/ ( {} | {} )'.format(gcopy, start, length)

    funkshin = gkcopy[gcopy]
    squeue =  funkshin(int(start), int(length))
    _col = '\033[38;5;201;48;5;154m'
    _msg = """\r
    \r{}>>>--Start( {} )--Length( {} )-->>>\033[00m
    \n{}>>>--[ {} ]-->>>\033[00m
     """.format(_col, start, length, _col, squeue)

    stdout.write(_msg)
    stdout.flush()


if __name__ == "__main__":
    _test_cases(*parse_args())






















