from argparse import ArgumentParser
from random import random
from math import log, log1p
import re

def parse_args():
    _parser = ArgumentParser(description='Parse cryptos')
    _parser.add_argument('-i')
    _parsed = _parser.parse_args()
    return int(_parsed.i)


def solution1(_total_lambs):
    def multiplexer(_total, _payout=1, _prev=None):
        if _total < _payout:
            return 0
        elif _prev is not None:
            return 1 + multiplexer(_total - _payout, _payout + _prev, _payout)
        else:
            return 1 + multiplexer(_total - _payout, _payout * 2, None)
        return 1 + multiplexer(*_temp)
    stingy = multiplexer(_total_lambs, _prev=0)
    generous = multiplexer(_total_lambs, _prev=None)
    return stingy - generous


def solution2(_total_lambs):
    def multiplexer(_total, _payout=1, _prev=None):
        if _total < _payout:
            return 0
        try:
            _tmp = (_payout + _prev, _payout)
        except TypeError:
            _tmp = (_payout + _payout, None)
        finally:
            return 1 + multiplexer(_total - _payout, *_tmp)
    stingy = multiplexer(_total_lambs, _prev=0)
    generous = multiplexer(_total_lambs, _prev=None)
    return stingy - generous


def solution3(_total_lambs):
    def stingy(_total, _payout=1, _prev=0):
        if _total <= _payout:
            return 0
        else:
            return 1 + stingy(_total - _prev, _payout + _prev, _payout)
    generous = int(log(_total_lambs, 2) // 1.0000001)
    return stingy(_total_lambs) - generous


if __name__ == "__main__":

    _ARGS = parse_args()
    _regx = re.compile('solution[0-9]')
    _globals = globals().copy()
    _gcopy = list(filter(lambda igx: re.match(_regx, igx), _globals.keys()))

    for _func in sorted(_gcopy):
        print '\nUsing( {} ): Solving[ {} ] -> [ {} ]'.format(
                _func,
                _ARGS,
                _globals[_func](_ARGS),
                )


























