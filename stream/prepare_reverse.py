#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 dal22 <dal22@fenrir.ess.washington.edu>
#
# Distributed under terms of the MIT license.

"""
Make a file to return to sliding steady state
"""

import sys

name = sys.argv[1]

template_fn = 'template_todeform.sif'

bm = float(name.split('_')[3][2:])

ACC = 0.1 + bm

subs = {'fn': name, 'rc': name.split('_')[1][2:], 'bm_val': name.split('_')[3][2:], 'ACC': str(ACC)}
print(subs)

with open(template_fn, 'r') as fin:
    with open('d_' + name + '.sif', 'w') as fout:
        for line in fin:
            for sub in subs:
                if '{%s}' % sub in line:
                    line = line.format(**subs)
                    break
            fout.write(line)
