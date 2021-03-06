#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 dal22 <dal22@fenrir.ess.washington.edu>
#
# Distributed under terms of the MIT license.

"""

"""
import sys

convn = sys.argv[1]
conv = convn
if conv.lower() == 'ftw':
    ACC = 0.1714
    VFUNC = 'VSSA2000M017My'
elif conv.lower() in ['ftwuc', 'uc']:
    convn = 'ftw_uc'
    ACC = 0.1
    VFUNC = 'VSSA2000M01My'
elif conv.lower() == 'hiacc':
    ACC = 0.6
    VFUNC = 'VSSA2000M06My'
else:
    raise ValueError('Unknown conv')

rc = sys.argv[2]

oop = sys.argv[3]
oop_line = '  OOPlane23 Strain Name = String "OOP23"\n  OOPlane13 Strain Name = String "OOP13"'
DUDY_LINE = '  dudy = Real 0.0'

bm = sys.argv[4]
if float(bm) == 0.0:
    bm_val = '1.0e-12'
else:
    bm_val = bm
    if bm_val == '1.0e-3':
        VFUNC = VFUNC + '001bm'
    elif bm_val == '1.0e-2':
        VFUNC = VFUNC + '01bm'
    elif bm_val == '2.0e-2':
        VFUNC = VFUNC + '02bm'
    else:
        raise ValueError('Unknown bm')
ACC = float(ACC) + float(bm)

subs = {'conv': conv, 'rc':rc, 'oop': oop, 'bm': bm, 'ACC':ACC, 'VFUNC':VFUNC, 
        'oop_lines': oop_line, 'bm_val':bm_val, 'DUDY_LINE':DUDY_LINE, 'convn': convn}

template_fn = 'template_fromoop.sif'
with open(template_fn, 'r') as fin:
    with open('a_{conv}_rc{rc}_oop{oop}_bm{bm}.sif'.format(conv=conv, rc=rc, oop=oop, bm=bm), 'w') as fout:
        for line in fin:
            for sub in subs:
                if '{%s}' % sub in line:
                    line = line.format(**subs)
                    break
            fout.write(line)

template_fn = 'template_oop2deform.sif'
with open(template_fn, 'r') as fin:
    with open('d_{conv}_rc{rc}_oop{oop}_bm{bm}.sif'.format(conv=conv, rc=rc, oop=oop, bm=bm), 'w') as fout:
        for line in fin:
            for sub in subs:
                if '{%s}' % sub in line:
                    line = line.format(**subs)
                    break
            fout.write(line)
