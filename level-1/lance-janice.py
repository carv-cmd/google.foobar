from argparse import ArgumentParser
from string import ascii_lowercase


def parse_args():
    _parser = ArgumentParser(description='Parse cryptos')
    _parser.add_argument('-i')
    return _parser.parse_args()


def solution(_input):
    print '\n> Encoded: {}'.format(str(_input.i))
    _alpha = ascii_lowercase
    _revr = list(_alpha)
    _revr.reverse()
    _codex = {_alpha[pi]: _revr[pi] for pi in xrange(len(_alpha))}

    _decoded = ''
    for _char in list(_input.i):
        if _char.islower():
            _decoded += _codex[_char]
        else:
            _decoded += _char

    return _decoded


_solute = solution(parse_args())
print '\n>>> Decoded: {}'.format(_solute)





