#_CLEAR = "\033[00m"
#_COLOR = "\r\033[2K" + "\033[38;5;11m"
#_MSG = _COLOR + '>>> Loading: '
#_printer = _MSG + '{} / {}'.format(_idx, len(_xs) - 1) + _CLEAR
#stdout.write(_printer)
#stdout.flush()

from argparse import ArgumentParser
import random as rnd
from datetime import datetime

from collections import defaultdict

#####################################################################
def parse_args():
    _parser = ArgumentParser()
    _parser.add_argument('-i')
    _parser.add_argument('-m')
    _parsed = _parser.parse_args()
    return int(_parsed.i), str(_parsed.m)
 

#####################################################################
def gen_list():
    with open('array-file.txt', 'r+') as fuckhead:
        _ttitti = [rnd.randint(-10, 10) for ipx in range(int(n))]
        for i in _ttitti:
            fuckhead.write('%s\n' % i)


#####################################################################
def quick_solar00(xs, _parity_chk=-1000):
    _neg = []
    _pos = []
    print 'Preprocess: 00'
    for _elem in xs:
        if _elem > 0:
            _pos.append(_elem)
        elif _elem < 0:
            _neg.append(_elem)
            if _elem > _parity_chk:
                _parity_chk = _elem

    print 'Parity Check: 00'
    if len(_neg) % 2:
        _neg.remove(_parity_chk)

    try:
        print 'Reduce: 00'
        _redux = reduce(lambda x, y: x * y, _neg + _pos)
    except TypeError:
        _redux = 0
    return str(_redux)

def quick_solar01(axs): 
    if len(axs) <= 1:
        _redox = axs.pop()
    else:
        _redox = quick_solar00(axs)
    return str(_redox)


#####################################################################
def quick_solar02(_xsarray):
    _zref = len(_xsarray)
    _absmin = -1001
    _product = 1
    if _zref == 1:
        _product = _xsarray[-1]
    else:
        print 'Preprocess: 02'
        for _elem in _xsarray:
            if not _elem:
                _zref -= 1
                continue
            if _absmin < _elem < 0: 
                _absmin = _elem
            _product *= _elem

        print 'Parity Check: 02'
        if _zref == 0:
            _product = 0
        elif _product < 0:
            _product //= int(_absmin)
    return _product


#####################################################################
def quick_solar03(_xsarray):
    if len(_xsarray) == 1:
        _product = str(_xsarray[-1])
    else:
        _dikt = defaultdict(int)
        for _elem in _xsarray:
            if _elem not in _dikt:
                _dikt[_elem] = _elem
            else:
                _dikt[_elem] *= _elem
        try:
            _filth = filter(lambda v: v != 0, _dikt.values())
            assert any(_filth)
            _product = reduce(lambda x, y: x * y, _filth)
            if _product < 0:
                _product //= max(filter(lambda m: m < 0, _dikt.keys()))
        except AssertionError:
            _product = 0
    return _product


#####################################################################
def quick_solar04(_xsarray):
    if len(_xsarray) == 1:
        _product = str(_xsarray[-1])

    else:
        _dikt = defaultdict(list)
        for _elem in _xsarray:
            if _elem not in _dikt:
                _dikt[_elem] = [_elem]
            else:
                _dikt[_elem].append(_elem)

        _keys = _dikt.keys()
        _filth = filter(lambda v: v != 0, _keys)
        _genxy = [fixy ** len(_dikt[fixy]) for fixy in _filth]
        _product = reduce(lambda x, y: x * y, _genxy)

        if _product < 0:
            if len(_genxy) == 1:
                _product = 0
            else:
                _product //= max(filter(lambda m: m < 0, _keys))

    return _product


#####################################################################
def get_time():
    _ts = datetime.now()
    return _ts.strftime('%H:%M:%S.%f')


def _test_cases(n=0, _mode=None):
    _iC='\033[38;5;{}m{}\033[00m'
    _exes = ([2, -3, 1, 0, -5], 
            [2, 0, 2, 2, 0], 
            [-2, -3, 4, -5], 
            [-2,], 
            [0, 0, 0, 0, 0, -3, 0, 0, 0, 0, 0, 0])

    if n != 0:
        _exes = []
        with open('array-file.txt', 'r+') as filehand:
            for ipx, ttitti in enumerate(filehand):
                if ipx > n:
                    break
                _exes.append(int(ttitti[:-1]))
        _exes = [_exes]

    _cols = 14
    print '\n>> Populated Array:\t\033[38;5;10m{}\033[00m'.format(get_time())
    for subset in _exes:
        if _mode in '1':
            print _iC.format(_cols, '\n>>> Executing: [ 1 ] ')
            _SOLAR = quick_solar01(subset)
        elif _mode in '2':
            print _iC.format(_cols, '\n>>> Executing: [ m: 2 | n : ' + str(n) + ' ]')
            _SOLAR = quick_solar02(subset)
        elif _mode in '3':
            print _iC.format(_cols, '\n>>> Executing: [ m: 3 | n : ' + str(n) + ' ]')
            _SOLAR = quick_solar03(subset)
        else:
            print _iC.format(_cols, '\n>>> Executing: [ m: 4 | n : ' + str(n) + ' ]')
            _SOLAR = quick_solar04(subset)
        if n != 0:
            print '\n<<< ohfuk.Coo.com @:\t\033[38;5;10m{}\033[00m'.format(get_time())
        else:
            print _iC.format(13, '\n>>> ARRAY-CONF >>>\t{}'.format(subset))
            print _iC.format(10, '\n>>>> ZERO-DAYS >>>\t[ {} ]\n'.format(_SOLAR))


#####################################################################
if __name__ == "__main__":
    print '\n> Start SolarWinds:\t\033[38;5;10m{}\033[00m'.format(get_time())
    _test_cases(*parse_args())
















































