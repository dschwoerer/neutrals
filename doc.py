#!/usr/bin/env python3

from glob import glob

src = glob("*.cxx")

hdr = glob("*.hxx")
# opts=[]
defaults = {}


models = []

var = {}

for f in src:
    with open(f) as s:
        for l in s:
            if "::" in l:
                cls = l.split("::")[0].split()[-1]
                if cls not in defaults:
                    defaults[cls] = {}
            if "OPTION" in l:
                #print(l, end='')
                t = l.split(',')
                defaults[cls][t[1]] = t[2].split(')')[0]
# done

# print(defaults)

parents = {}
bases = []
coms = {}

for f in hdr:
    cc = None
    com = ''
    with open(f) as s:
        for l in s:
            if "class" in l:
                if ";" in l:
                    continue
                if '///' in l:
                    continue
                if ":" in l:
                    # is derived
                    l2 = l.split(':')
                    # print(f,l2)
                    cls = l2[0].strip()[6:].strip()
                    base = l2[1].strip()[6:-1].strip()
                    parents[cls] = base
                    # print(cls,base)
                else:
                    cls = l.strip()[6:-1].strip()
                    bases.append(cls)
                    # print(cls)
                    # pass
                # print('var',var)
                # print('cls',cls)
                var[cls] = {}
                cc = cls
                com = ''
            elif cc:
                if '///<' in l:
                    ls = l.split('///<')
                    com = ls[1]
                    vari = ls[0]
                    # print(vari,com)
                    var[cc][vari] = com
                elif '///' in l:
                    com += l.split('///')[1]
                elif com:
                    vari = l
                    # print(vari,com)
                    var[cc][vari] = com
                    com = ''

# print('var',var)


def print_class(cls=None, lvl=0, mode=None):
    if cls is None:
        for c in bases:
            print_class(c, mode=mode)
    else:
        print('    ' * lvl, cls)
        if mode == 'members':
            if cls in var:
                for v in var[cls]:
                    print('    ' * lvl, '*', v, end='')
                    for l in var[cls][v].split('\n'):
                        print('    ' * lvl, '   ', l)
        if mode == 'defaults':
            if cls in defaults:
                for v in defaults[cls]:
                    print('    ' * lvl, '*', "%s [%s]" % (v, defaults[cls][v]))
                    if cls in var:
                        for v2 in var[cls]:
                            if v in v2:
                                for l in var[cls][v2].split('\n'):
                                    if l:
                                        print('    ' * lvl, '   ', l)
        for k in parents:
            if parents[k] == cls:
                print_class(k, lvl + 1)

if __name__ == "__main__":
    import sys
    if sys.argv[1] == 'def':
        print_class('Neutrals', mode='defaults')
    elif sys.argv[1] == 'mem':
        print_class('Neutrals', mode='members')
    elif sys.argv[1] == 'tree':
        print_class()
    else:
        sys.exit(1)
else:
    print(__name__)


# for
# print(src)

# for f in
