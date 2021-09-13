#!/usr/bin/env python
"""Tools for data handling."""
#
#
#
#  Created by Hartmut Lentz on 13.05.13.
#


def cdf2histogram(c_in):
    """Read cdf as list and returns histogram."""
    h = []
    h.append(c_in[0])
    for i in range(1, len(c_in)):
        h.append(c_in[i] - c_in[i-1])
    return h


def dict2file(d, nameoffile='dict.txt', sorted=True):
    """
    Write dictionary (or list or tuple) to a textfile.

    Sorted by keys, if sorted=True.
    """
    def list2dict(li):
        x = {}
        for (i, el) in enumerate(li):
            x[i] = el
        return x

    if not isinstance(d, dict):
        d = list2dict(d)

    dk = list(d.keys())
    if sorted:
        dk.sort()

    # if d={ 1: [a,b,c,...], 2:[d,e,f,...],... }
    s = list(d.values())[0]
    if isinstance(s, dict) or isinstance(s, list) or isinstance(s, tuple):
        laenge = len(d.values()[0])
        with open(nameoffile, "w+") as g:
            # g = file(nameoffile, 'w+')
            for k in dk:
                wstring = ''
                for i in range(laenge):
                    wstring += '\t' + str(d[k][i])
                g.writelines((str(k) + wstring + '\n'))
            g.close
        return

    g = open(nameoffile, 'w+')
    for k in dk:
        g.writelines((str(k), '\t', str(d[k]), '\n'))
    g.close

    return


if __name__ == "__main__":
    pass
