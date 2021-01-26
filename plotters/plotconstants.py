#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2017 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""
I want it to be easier to choose constants
"""

import matplotlib as mpl

mpl_ver = mpl.__version__.split('.')
if int(mpl_ver[0]) < 2 and int(mpl_ver[1]) < 9:
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
else:
    colors = ['C{:d}'.format(d) for d in range(10)]

colors = ['#e41a1c',
          '#377eb8',
          '#4daf4a',
          '#984ea3',
          '#ff7f00',
          '#ffff33',
          '#a65628',
          '#f781bf']

colors = ['#1b9e77', '#d95f02', '#e7298a', '#e6ab02','#a6761d','#666666', '#7570b3', '#66a61e']

colors = ['maroon', 'goldenrod', 'mediumseagreen', 'royalblue'] + colors

class Constants():
    def __init__(self):
        self.density = None
        self.fontsize = None
        self.medsize = None
        self.bigfontsize = None
        self.fontname = None
        self.full_font_black = None
        self.folder = None
        self.std_leg_width = None
        self.filetype = None
        self.lettering = None
        self.cbar_height = None


class AnimationConstants(Constants):
    def __init__(self):
        self.density = 220
        self.fontsize = 16
        self.medfontsize = 20
        self.bigfontsize = 32
        self.fontname = 'Computer Modern Sans Serif'
        self.fontfamily = 'sans-serif'
        self.full_font_black = '{:d}p,{:s},black'.format(self.fontsize, self.fontname)
        self.folder = 'forward'
        self.std_leg_width = 1.5
        self.filetype = 'mp4'
        self.lettering = ['a', 'b', 'c', 'd', 'e', 'f']
        self.cbar_height = 0.2


class PrintedConstants(Constants):
    def __init__(self):
        self.density = 300
        self.fontsize = 8
        self.medfontsize = 12
        self.bigfontsize = 16
        self.fontname = 'Times-Roman'
        self.fontfamily = 'serif'
        self.full_font_black = '{:d}p,{:s},black'.format(self.fontsize, self.fontname)
        self.folder = 'pub'
        self.std_leg_width = 1.5
        self.filetype = 'png'
        self.lettering = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
        self.cbar_height = 0.1


class ScreenConstants(Constants):
    def __init__(self):
        self.density = 450
        self.fontsize = 12
        self.medfontsize = 14
        self.bigfontsize = 20
        self.fontname = 'Computer Modern Sans Serif'
        self.fontfamily = 'sans-serif'
        self.full_font_black = '{:d}p,{:s},black'.format(self.fontsize, self.fontname)
        self.folder = 'slideshow'
        self.std_leg_width = 1.5
        self.filetype = 'pdf'
        self.lettering = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
        self.cbar_height = 0.1


def load_constants(cname, **kwargs):
    print('loading constants')
    if cname in ['serif', 'print', 'pub']:
        constants = PrintedConstants()
    elif cname in ['animation', 'movie']:
        constants = AnimationConstants()
    elif cname in ['sans', 'sans-serif', 'screen']:
        constants = ScreenConstants()
    else:
        raise ValueError('bad cname')

    for key, value in kwargs.items():
        setattr(constants, key, value)

    if False:
        try:
            mpl.use('macosx')
        except:
            if constants.filetype == 'pdf':
                mpl.use('PDF')
            elif constants.filetype == 'png':
                mpl.use('agg')
            else:
                mpl.use('TKagg')

    mpl.rc('text', usetex=True)
    mpl.rc('font', **{'family': constants.fontfamily, constants.fontfamily: constants.fontname})
    if constants.fontfamily == 'sans-serif':
        mpl.rc('text.latex', preamble=r'\usepackage{color} \usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault} \usepackage{cmbright}')
    else:
        mpl.rc('text.latex', preamble=r'\usepackage{color}')
    try:
        mpl.rcParams['patch.force_edgecolor'] = True
    except KeyError:
        # This is going to happen with some old versions of mpl
        pass

    return constants
