#! /usr/bin/env python
# coding: utf-8

import sys
from modeltools.lib import fastvtulib
from scipy.interpolate import interp1d
import numpy as np
from modeltools.lib import gmshlib

if len(sys.argv) >= 2:
    rescale = float(sys.argv[1])
else:
    rescale = 10.

if len(sys.argv) >= 3:
    if sys.argv[2] == 'stream':
        in_files = ['stream/deform_merged0076.vtu']
        out_files = ['relaxed_stream_tall.geo']
        reses = [500]
    elif sys.argv[2] == 'relaxed_stream_tall':
        in_files = ['relaxed_stream_tall/deform0151.vtu']
        out_files = ['relaxed_stream_tall_fabric.geo']
        reses = [250]
    elif sys.argv[2] == 'rstf':
        in_files = ['relaxed_stream_tall_fabric/deform_diffuse_rr0545.vtu']
        out_files = ['rstf.geo']
        reses = [250]
    elif sys.argv[2] == 'rstf2':
        in_files = ['rstf/deform2e-3_0102.vtu']
        out_files = ['rstf.geo']
        reses = [250]
    elif sys.argv[2] == 'stream_ftw':
        in_files = ['rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_1002.vtu']
        out_files = ['stream_ftw.geo']
        reses = [250]
    else:
        raise ValueError('Which one')
else:
    in_files = ['stream/deform_merged0076.vtu', 'relaxed_stream_tall/deform0151.vtu']
    out_files = ['relaxed_stream_tall.geo', 'relaxed_stream_tall_fabric.geo']
    reses = [500, 250]

for in_file, out_file, res in zip(in_files, out_files, reses):
    a = fastvtulib.get_structured_vtu(in_file)
    tcoordsx = a.rawdata_dict['coordsX'][a.rawdata_dict['zs top'] > 0]
    tcoordsy = a.rawdata_dict['coordsY'][a.rawdata_dict['zs top'] > 0]
    interper = interp1d(tcoordsx, tcoordsy)

    x = 100000
    # Double the res here so we don't have artificial fining
    xv = np.arange(x, -x - res * 4., -res * 4.)
    yv = interper(xv)
    yv *= rescale

    top = np.vstack((xv, yv)).T
    gmshlib.gmsh_outline(out_file, [[(-x, 0)], [(x, 0)], top[:-1, :], [top[-1, :]]], res, cuts=None, cuts_lc=None, points=None, points_lc=None, spline_or_line='Line')
