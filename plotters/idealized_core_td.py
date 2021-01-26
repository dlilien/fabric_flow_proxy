#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 dal22 <dal22@fenrir.ess.washington.edu>
#
# Distributed under terms of the MIT license.

"""
Plots with time on the x axis, something on the y.
"""
from matplotlib import rc
rc('font',**{'family':'sans-serif', 'sans-serif':['Linux Biolinum O']})

import os.path
from glob import glob
import argparse
import pickle

from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt
from modeltools.lib import fastvtulib
import numpy as np
import matplotlib.gridspec as gridspec
from lib import anisolib, fabricplotlib

debug = False

fs = 9
bfs = 12
colors = ['C0', 'C1', 'C2']
oop_colors = ['C3', 'C4', 'C5']
name_fmt = 's_{:s}_rc{:s}_oop{:s}_bm{:s}'
b_name_fmt = 'd_{:s}_rc{:s}_oop{:s}_bm{:s}'
ob_name_fmt = 'a_{:s}_rc{:s}_oop{:s}_bm{:s}'
name_dict = {'ftw': 'Converge\nupstream',
             'uc': 'Converge\neverywhere',
             'hiacc': 'High\naccumulation'}
short_name_dict = {'ftw': 'Conv. upstream',
                   'uc': 'Conv. everywhere',
                   'hiacc': 'High acc.',
                   'nuc': 'NUC',
                   'nuclc': 'NUCLC',
                   'ftwlc': 'FTW LC',
                   'midacc': 'Mid Acc'}


def get_vars_or_cache(x, y, fmt, fn, folder='../stream/rstf'):
    """Check modification times and load cache if available.

    Exists because loading is really slow since it is many, many Gb.
    If no .result, we try just the cache. Unhandled exception if neither
    exists.
    """
    cache_fn = os.path.join(folder, fmt + '.pickle')
    result_fn = os.path.join(folder, fmt + '.result')
    if (not os.path.exists(result_fn)) or (os.path.exists(cache_fn) and (
            os.path.getmtime(cache_fn) > os.path.getmtime(result_fn))):
        print('Using cached value for {:s}'.format(fmt))
        dat = pickle.load(open(cache_fn, 'rb'))
    else:
        td_vtu = fastvtulib.get_td_vtu(fn)
        dat = td_vtu.get_pts_2d(['fabric 1', 'fabric 2', 'fabric 3',
                                 'fabric 4', 'fabric 5',
                                 'eigenv 1', 'eigenv 2', 'eigenv 3',
                                 'eigenv 4'], x, y)
        pickle.dump(dat, open(cache_fn, 'wb'))
    return dat


def oop(poster=False):
    """Plot with out of page components."""
    fig = plt.figure(constrained_layout=True, figsize=(7.8, 5.5))
    gs = gridspec.GridSpec(nrows=4, ncols=10, hspace=0.5, wspace=0.0,
                           left=0.075, bottom=0.07, top=0.91, right=0.95,
                           width_ratios=(1, 0.25, 0.15, 0.05, 0.03, 0.17, 0.35, 0.05, 0.03, 0.12),
                           height_ratios=(0.5, 0.5, 0.9, 0.2))
    ax1 = fig.add_subplot(gs[0:2, 0:3])
    ax2 = fig.add_subplot(gs[2, 0])

    axb1 = fig.add_subplot(gs[0, 6], projection=fabricplotlib.PRJ)
    axb2 = fig.add_subplot(gs[1, 6], projection=fabricplotlib.PRJ)
    axb_cb = fig.add_subplot(gs[0:2, 8])

    ax8 = fig.add_subplot(gs[2, 2:10], sharex=ax2)
    ax_leg = fig.add_subplot(gs[3, :])
    ax_leg.axis('off')
    dist = 75.0
    x, y = np.ones(3) * 1000.0 * dist, np.array([1000., 1500.])

    ftws = ['ftw', 'uc', 'hiacc']
    oops = ['1.0e-3', '1.0e-2', '1.0e-1']
    bms = ['0.0', '1.0e-2', '2.0e-2']
    rcs = ['2.0e-3', '2.0e-4', '2.0e-2']

    oop_names = [name_fmt.format(ftws[0], rcs[0], oop, bms[0]) for oop in oops]
    bck_names = [ob_name_fmt.format(ftws[0], rcs[0], oop, bms[0]) for oop in oops]

    fmts = oop_names + bck_names
    pretty_fmts = {}
    pretty_fmts_oop = {}
    for fmt in fmts:
        fm_spl = fmt.split('_')
        pretty_fmts[fmt] = r'%s $\lambda$=%1.4g $\dot b$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[2][2:]), float(fm_spl[4][2:]))
        if fmt[0] == 's':
            pretty_fmts_oop[fmt] = r'$\dot\varepsilon_{xy}^{(max)}$=%1.4g Fwd' % (float(fm_spl[3][3:]))
        else:
            pretty_fmts_oop[fmt] = r'$\dot\varepsilon_{xy}^{(max)}$=%1.4g Bck' % (float(fm_spl[3][3:]))


    if not debug:
        files = [glob('../stream/stream_ftw/' + fmt + '_????.vtu') for fmt in fmts]
        inds = [np.argsort(np.array([float(fn[-8:-4]) for fn in filel])) for filel in files]
        fns = [[file_list[ind] for ind in inds_i] for inds_i, file_list in zip(inds, files)]

        times = np.hstack(([0.0, 0.1], np.arange(10.0, 100010.0, 10.0)))

        a_s = {}
        for i, fmt in enumerate(fmts):
            a_s[fmt] = get_vars_or_cache(x, y, fmt, fns[i], folder='../stream/stream_ftw')
        timess = {fmt: times[:min(a_s[fmt]['fabric 1'].shape[1], len(times))] for fmt, fn in zip(fmts, fns)}
        taus = {name: val / (dist * 1000.0 / 50.0) for name, val in timess.items()}

    else:
        taus = {name: np.linspace(0, 3, 100) for name in fmts}
        vtu = {name: np.ones((3, 100)) for name in ['fabric 1', 'fabric 2', 'fabric 3', 'eigenv 1', 'eigenv 2', 'eigenv 3', 'eigenv 4']}
        a_s = {name: vtu for name in fmts}

    def do_plot(fmt, color, lw, ls):
        if fmt[0] == 's':
            label = pretty_fmts_oop[fmt]
        else:
            label = None
        a_s[fmt]['eigenv 3'][1, :][a_s[fmt]['eigenv 3'][1, :] > 1.0] = 1.0
        ax2.plot(timess[fmt] / 1000.0, a_s[fmt]['eigenv 3'][1, :], color=color, linewidth=lw, linestyle=ls, label=label)

        if fmt[0] == 'a':
            label = pretty_fmts_oop[fmt]
        else:
            label = None

        ax8.plot(timess[fmt][1:] / 1000.0, fabricplotlib.fabric_to_hor_rot(a_s[fmt]['fabric 1'][1, 1:],
                                                             a_s[fmt]['fabric 2'][1, 1:],
                                                             a_s[fmt]['fabric 5'][1, 1:]),
                 color=color, linewidth=lw, linestyle=ls, label=pretty_fmts_oop[fmt])

    for fmt, color in zip(oop_names, oop_colors):
        do_plot(fmt, color, 1, 'solid')

    for fmt, color in zip(bck_names, oop_colors):
        do_plot(fmt, color, 1, 'dashed')

    ax1.set_ylabel('Depth (m)', fontsize=fs)
    ax1.set_xlabel('Distance (km)', fontsize=fs)
    ax2.set_ylabel(r'$a^{(2)}_{1}$', fontsize=fs)
    ax8.set_ylabel(r'$\theta$', fontsize=fs)
    ax2.set_xlabel(r'Time (kyr)', fontsize=fs)
    ax8.set_xlabel(r'Time (kyr)', fontsize=fs)
    ax2.set_ylim(0.66666, 1.0)
    ax2.set_yticks([0.66666, 5. / 6., 1.])
    ax2.set_yticklabels([r'$\frac{2}{3}$', r'$\frac{5}{6}$', '1'])
    ax2.set_xlim(0, 3)
    ax8.set_ylim(0, 45)
    h, l = ax8.get_legend_handles_labels()
    ax_leg.legend([h[0], h[3], h[1], h[4], h[2], h[5]], [l[0], l[3], l[1], l[4], l[2], l[5]], loc='upper left', frameon=False, ncol=3, fontsize=fs)
    # ax2.legend(loc='upper right', frameon=False, fontsize=fs)
    # ax5.legend(loc='upper right', frameon=False, fontsize=fs)

    vtus = [fastvtulib.get_structured_vtu(fns[1][-1])]
    tris = [Triangulation(np.array([rc[0] + 100000. for rc in vtu.raw_coords[:, 0]]) / 1000., np.array([rc[1] for rc in vtu.raw_coords[:, 0]]), vtu.simptt) for vtu in vtus]

    a12_axes = [ax1]
    ax_c_a12 = fig.add_subplot(gs[0:2, 4])

    for axa, tri, vtu in zip(a12_axes, tris, vtus):
        axa.set_xlim(0, 175)
        axa.set_xticks([0., 50., 100., 150.])
        axa.set_xticklabels(['' for tick in axa.get_xticklabels()])

        cm3 = axa.tricontourf(tri, vtu.rawdata_dict['eigenv 3'], cmap='summer', levels=np.linspace(0.3333, 1, 101), extend='neither')
        for c in cm3.collections:
            c.set_edgecolor("face")
        axa.set_xticklabels(['0', '50', '100', '150'])
        fabricplotlib.quiver(axa, vtu, scale=25, width=0.003)

    a12_axes[0].scatter(-1000, -1000, marker=r'$\uparrow$', label='Single max. in x-z', color='k')
    a12_axes[0].legend(loc='lower left', bbox_to_anchor=(0.1, 1.0), ncol=2, fontsize=fs, framealpha=1.0)
    a12_axes[0].set_xlim(0, 175)
    a12_axes[0].set_ylim(0, 2200)

    cbr = plt.colorbar(cm3, cax=ax_c_a12, orientation='vertical', ticks=(1. / 3., 2. / 3., 1.))
    cbr.set_label(label=r'$a^{(2)}_{1}$', size=fs)
    cbr.ax.set_yticklabels([r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)

    # Cartoons
    x, y = np.array([40000, 40000]), np.array([1100, 150])
    ax1.text(x[0] / 1000. + 100, y[0], 'b', fontsize=fs, ha='center', va='center', bbox=dict(boxstyle='square,pad=0.1', facecolor='white', alpha=0.75))
    ax1.text(x[1] / 1000. + 100, y[1], 'c', fontsize=fs, ha='center', va='center', bbox=dict(boxstyle='square,pad=0.1', facecolor='white', alpha=0.75))
    fab_at_pts = vtu.get_pts_2d(anisolib.fabs, x, y)
    a2 = anisolib.fabric_dict_to_a2(fab_at_pts)

    for letter, ax in zip('bcde', [axb1, axb2]):
        ax.text(0.00, 0.9, letter, transform=ax.transAxes, fontsize=bfs)

    fabricplotlib.a2plot(a2[:, :, 0], ax=axb1, cbr=False, show=False, levels=13)
    cm = fabricplotlib.a2plot(a2[:, :, 1], ax=axb2, cbr=False, show=False, levels=13)
    cbr = plt.colorbar(cm, cax=axb_cb, orientation='vertical')
    cbr.set_label(label=r'ODF($\theta,\phi$)', size=fs)
    # cbr.ax.set_xticklabels(['0', '1', '2'])
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)

    for letter, ax in zip('adefghijklmnop', (ax1, ax2, ax8)):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis='both', which='major', labelsize=fs)

    fig.savefig('../plots/idealized_core_oop_td.png', dpi=300)
    fig.savefig('../plots/poster_idealized_core_oop_td.png', dpi=300, transparent=True)


def slow():
    x, y = np.array([75000., 75000., 75000., 50000., 50000., 50000., 25000., 25000., 25000.]), np.array([1500., 1000., 500., 1500., 1000., 500., 1500., 1000., 500.])

    ftws = ['ftw', 'uc', 'hiacc', 'nuc', 'nuclc', 'ftwlc', 'midacc']
    oops = ['0.0', '1.0e-3', '1.0e-2', '1.0e-1']
    bms = ['0.0', '1.0e-2', '2.0e-2']
    rcs = ['2.0e-3', '2.0e-4', '2.0e-2']

    ftw_names = [name_fmt.format(ftw, rcs[0], oops[0], bms[0]) for ftw in ftws]
    bm_names = [name_fmt.format(ftws[0], rcs[0], oops[0], bm) for bm in bms[1:]]
    rc_names = [name_fmt.format(ftws[0], rc, oops[0], bms[0]) for rc in rcs[1:]]
    oop_names = [name_fmt.format(ftws[0], rcs[0], oop, bms[0]) for oop in oops[1:]]

    fmts = ftw_names + bm_names + rc_names + oop_names
    pretty_fmts = {}
    pretty_fmts_oop = {}
    for fmt in fmts:
        fm_spl = fmt.split('_')
        pretty_fmts[fmt] = r'%s $\lambda$=%1.4g $\dot b$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[2][2:]), float(fm_spl[4][2:]))
        pretty_fmts_oop[fmt] = r'%s $\dot\varepsilon_{xy}^{(2)}$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[3][3:]))

    linewidths = [2, 1, 0.5]
    linestyles = ['solid', 'dashed', 'dotted']

    if not debug:
        files = [glob('../stream/rstf/' + fmt + '_????.vtu') for fmt in fmts]
        inds = [np.argsort(np.array([float(fn[-8:-4]) for fn in filel])) for filel in files]
        fns = [[file_list[ind] for ind in inds_i] for inds_i, file_list in zip(inds, files)]

        times = np.hstack(([0.0, 0.1], np.arange(10.0, 100010.0, 10.0)))

        a_s = {}
        for i, fmt in enumerate(fmts):
            a_s[fmt] = get_vars_or_cache(x, y, fmt, fns[i], folder='../stream/rstf')
        timess = {fmt: times[:min(a_s[fmt]['fabric 1'].shape[1], len(times))] for fmt, fn in zip(fmts, fns)}
        taus = {name: val / (75000.0 / 50.0) for name, val in timess.items()}
        taus25 = {name: val / (25000. / 50.0) for name, val in timess.items()}
        taus50 = {name: val / (50000. / 50.0) for name, val in timess.items()}

    else:
        timess = {name: np.linspace(0, 4500, 100) for name in fmts}
        taus = {name: np.linspace(0, 3, 100) for name in fmts}
        vtu = {name: np.ones((9, 100)) for name in ['fabric 1', 'fabric 2', 'fabric 3', 'eigenv 1', 'eigenv 2', 'eigenv 3', 'eigenv 4']}
        a_s = {name: vtu for name in fmts}

    def make_fig(name, dist_num=0, taus=taus):
        fig = plt.figure(constrained_layout=True, figsize=(8, 5))
        gs = gridspec.GridSpec(nrows=5, ncols=3, hspace=0.14, wspace=0.12, left=0.075, bottom=0.0, top=0.965, right=0.98, height_ratios=(1, 1, 1, 0.06, 0.67))
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
        ax3 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
        ax4 = fig.add_subplot(gs[1, 0], sharex=ax1)
        ax5 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax4)
        ax6 = fig.add_subplot(gs[1, 2], sharex=ax1, sharey=ax4)
        ax7 = fig.add_subplot(gs[2, 0], sharex=ax1)
        ax8 = fig.add_subplot(gs[2, 1], sharex=ax1, sharey=ax7)
        ax9 = fig.add_subplot(gs[2, 2], sharex=ax1, sharey=ax7)
        ax_leg = fig.add_subplot(gs[4, :])
        ax_leg.axis('off')

        def do_plot(fmt, color, lw, ls, dist_num=0, taus=taus, ax1=ax1, ax2=ax2, ax3=ax3, ax4=ax4, ax5=ax5, ax6=ax6, ax7=ax7, ax8=ax8, ax9=ax9, fmts=pretty_fmts):
            ax1.plot(taus[fmt], a_s[fmt]['fabric 1'][0 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax2.plot(taus[fmt], a_s[fmt]['fabric 1'][1 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax3.plot(taus[fmt], a_s[fmt]['fabric 1'][2 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

            ax4.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][0 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][0 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax5.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][1 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][1 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax6.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][2 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][2 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

            ax7.plot(taus[fmt], a_s[fmt]['fabric 2'][0 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax8.plot(taus[fmt], a_s[fmt]['fabric 2'][1 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax9.plot(taus[fmt], a_s[fmt]['fabric 2'][2 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

            # ax7.plot(taus[fmt], np.abs(np.log10(np.abs(a_s[fmt]['eigenv 3'][0 + 3 * dist_num, :] / a_s[fmt]['eigenv 2'][0, :])) / np.log10(np.abs(a_s[fmt]['eigenv 2'][0, :] / a_s[fmt]['eigenv 1'][0, :]))), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            # ax8.plot(taus[fmt], np.abs(np.log10(np.abs(a_s[fmt]['eigenv 3'][1 + 3 * dist_num, :] / a_s[fmt]['eigenv 2'][1, :])) / np.log10(np.abs(a_s[fmt]['eigenv 2'][1, :] / a_s[fmt]['eigenv 1'][1, :]))), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            # ax9.plot(taus[fmt], np.abs(np.log10(np.abs(a_s[fmt]['eigenv 3'][2 + 3 * dist_num, :] / a_s[fmt]['eigenv 2'][2, :])) / np.log10(np.abs(a_s[fmt]['eigenv 2'][2, :] / a_s[fmt]['eigenv 1'][2, :]))), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

        for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9):
            ax.axhline(1. / 3., color='k', linestyle='dotted')

        # for ax in (ax7, ax8, ax9):
        #     ax.axhline(1., color='k', linestyle='dotted')

        for fmt, color in zip(ftw_names, ['C0', 'C1', 'C2', 'C6', 'C7', 'C8', 'C9', 'g']):
            lw = linewidths[0]
            ls = linestyles[0]
            do_plot(fmt, color, lw, ls, taus=taus, dist_num=dist_num)

        ax1.set_ylabel(r'$a^{(2)}_{xx}$', fontsize=fs)
        ax4.set_ylabel(r'$a^{(2)}_{yy}$', fontsize=fs)
        ax7.set_ylabel(r'$a^{(2)}_{zz}$', fontsize=fs)
        # ax7.set_ylabel(r'Woodcock $k$', fontsize=fs)
        ax1.set_title('Depth = {:d} m'.format(int(2000.0 - y[0])), fontsize=fs)
        ax2.set_title('Depth = {:d} m'.format(int(2000.0 - y[1])), fontsize=fs)
        ax3.set_title('Depth = {:d} m'.format(int(2000.0 - y[2])), fontsize=fs)
        ax1.set_ylim(0, 1)
        if np.max(taus[fmts[0]]) < 10:
            ax1.set_xlim(0, 3)
            ax8.set_xlabel(r'Time / $\tau_s$', fontsize=fs)
        else:
            ax1.set_xlim(0, 10000)
            ax8.set_xlabel('Time (yrs)', fontsize=fs)
        ax4.set_ylim(0, 1)
        ax7.set_ylim(0, 1)
        # ax7.set_ylim(0, 2)

        handles, labels = ax7.get_legend_handles_labels()
        ax_leg.legend(handles, labels, loc='upper left', frameon=False, ncol=3, fontsize=fs)

        for ax in (ax1, ax2, ax3, ax4, ax5, ax6):
            ax.xaxis.set_tick_params(labelbottom=False)

        for ax in (ax2, ax3, ax5, ax6, ax8, ax9):
            ax.yaxis.set_tick_params(labelleft=False)

        for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)):
            ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
            ax.tick_params(axis='both', which='major', labelsize=fs)

        fig.tight_layout(pad=0.1)
        fig.savefig(name, dpi=300)
    make_fig('../plots/idealized_core_75_time_slow.pdf', dist_num=0, taus=timess)
    make_fig('../plots/idealized_core_25_time_slow.pdf', dist_num=2, taus=timess)
    make_fig('../plots/idealized_core_50_time_slow.pdf', dist_num=1, taus=timess)


def forward(poster=False):
    x, y = np.array([75000., 75000., 75000., 50000., 50000., 50000., 25000., 25000., 25000.]), np.array([1500., 1000., 500., 1500., 1000., 500., 1500., 1000., 500.])

    ftws = ['ftw', 'uc', 'hiacc']
    oops = ['0.0', '1.0e-3', '1.0e-2', '1.0e-1']
    bms = ['0.0', '1.0e-2', '2.0e-2']
    rcs = ['2.0e-3', '2.0e-4', '2.0e-2']

    ftw_names = [name_fmt.format(ftw, rcs[0], oops[0], bms[0]) for ftw in ftws]
    # bm_names = [name_fmt.format(ftws[0], rcs[0], oops[0], bm) for bm in bms[1:]]
    # rc_names = [name_fmt.format(ftws[0], rc, oops[0], bms[0]) for rc in rcs[1:]]
    oop_names = [name_fmt.format(ftws[0], rcs[0], oop, bms[0]) for oop in oops[1:]]

    if poster:
        fmts=ftw_names
    else:
        fmts = ftw_names + oop_names
    pretty_fmts = {}
    pretty_fmts_oop = {}
    for fmt in fmts:
        fm_spl = fmt.split('_')
        pretty_fmts[fmt] = r'%s $\lambda$=%1.4g $\dot b$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[2][2:]), float(fm_spl[4][2:]))
        if float(fm_spl[3][3:]) > 0:
            pretty_fmts_oop[fmt] = r'%s+margin, $\dot\varepsilon_{xy}^{(\max)}$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[3][3:]))
        else:
            pretty_fmts_oop[fmt] = r'%s $\dot\varepsilon_{xy}^{(\max)}$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[3][3:]))

    linewidths = [2, 1, 0.5]
    linestyles = ['solid', 'dashed', 'dotted']

    if not debug:
        files = [glob('../stream/rstf/' + fmt + '_????.vtu') for fmt in fmts]
        inds = [np.argsort(np.array([float(fn[-8:-4]) for fn in filel])) for filel in files]
        fns = [[file_list[ind] for ind in inds_i] for inds_i, file_list in zip(inds, files)]

        times = np.hstack(([0.0, 0.1], np.arange(10.0, 100010.0, 10.0)))

        a_s = {}
        for i, fmt in enumerate(fmts):
            a_s[fmt] = get_vars_or_cache(x, y, fmt, fns[i], folder='../stream/rstf')
        timess = {fmt: times[:min(a_s[fmt]['fabric 1'].shape[1], len(times))] for fmt, fn in zip(fmts, fns)}
        taus = {name: val / (75000.0 / 50.0) for name, val in timess.items()}

    else:
        timess = {name: np.linspace(0, 4500, 100) for name in fmts}
        taus = {name: np.linspace(0, 3, 100) for name in fmts}
        vtu = {name: np.ones((9, 100)) for name in ['fabric 1', 'fabric 2', 'fabric 3', 'eigenv 1', 'eigenv 2', 'eigenv 3', 'eigenv 4']}
        a_s = {name: vtu for name in fmts}

    def make_fig(name, dist_num=0, taus=taus):
        fig = plt.figure(constrained_layout=True, figsize=(8, 3.5))
        if poster:
            leg = 0.3
        else:
            leg = 0.5
        gs = gridspec.GridSpec(nrows=4, ncols=2, hspace=0.14, wspace=0.12, left=0.075, bottom=0.0, top=0.95, right=0.98, height_ratios=(1, 1, 0.1, leg))
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
        ax3 = fig.add_subplot(gs[1, 0], sharex=ax1)
        ax4 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax3)
        ax_leg = fig.add_subplot(gs[3, :])
        ax_leg.axis('off')

        def do_plot(fmt, color, lw, ls, dist_num=0, taus=taus, ax1=ax1, ax2=ax2, ax3=ax3, ax4=ax4, fmts=pretty_fmts_oop):
            ax1.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][0 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][1 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax2.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][1 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][2 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

            ax3.plot(taus[fmt], a_s[fmt]['fabric 2'][1 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax4.plot(taus[fmt], a_s[fmt]['fabric 2'][2 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

        for ax in (ax1, ax2, ax3, ax4):
            ax.axhline(1. / 3., color='k', linestyle='dotted')

        # for ax in (ax7, ax8, ax9):
        #     ax.axhline(1., color='k', linestyle='dotted')

        for fmt, color in zip(ftw_names, colors):
            lw = linewidths[0]
            ls = linestyles[0]
            do_plot(fmt, color, lw, ls, taus=taus, dist_num=dist_num)

        if not poster:
            for fmt, color in zip(oop_names[1:], oop_colors[1:]):
                lw = linewidths[0]
                ls = linestyles[0]
                do_plot(fmt, color, lw, ls, taus=taus, dist_num=dist_num, fmts=pretty_fmts_oop)

        ax1.set_ylabel(r'$a^{(2)}_{yy}$', fontsize=fs)
        ax3.set_ylabel(r'$a^{(2)}_{zz}$', fontsize=fs)
        # ax7.set_ylabel(r'Woodcock $k$', fontsize=fs)
        ax1.set_title('Depth = {:d} m'.format(int(2000.0 - y[1])), fontsize=fs)
        ax2.set_title('Depth = {:d} m'.format(int(2000.0 - y[2])), fontsize=fs)

        ax1.set_ylim(0, 1.0)
        if np.max(taus[fmts[0]]) < 10:
            ax1.set_xlim(0, 3)
            ax3.set_xlabel(r'Time / $\tau_s$', fontsize=fs)
            ax4.set_xlabel(r'Time / $\tau_s$', fontsize=fs)
        else:
            ax1.set_xlim(0, 5000)
            ax3.set_xlabel('Time (yrs)', fontsize=fs)
            ax4.set_xlabel('Time (yrs)', fontsize=fs)
        ax3.set_ylim(0, 1)

        handles, labels = ax4.get_legend_handles_labels()
        ax_leg.legend(handles, labels, loc='upper left', frameon=False, ncol=3, fontsize=fs)

        for ax in (ax1, ax2):
            ax.xaxis.set_tick_params(labelbottom=False)

        for ax in (ax2, ax4):
            ax.yaxis.set_tick_params(labelleft=False)

        for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4)):
            ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
            ax.tick_params(axis='both', which='major', labelsize=fs)

        fig.tight_layout(pad=0.1)
        fig.savefig(name, dpi=300, transparent=poster)
    # make_fig('../plots/idealized_core_75_tau.pdf', dist_num=0, taus=taus)
    if poster:
        make_fig('../plots/poster_idealized_core_75_time.png', dist_num=0, taus=timess)
    else:
        make_fig('../plots/idealized_core_75_time.pdf', dist_num=0, taus=timess)


def backward(poster=False):
    x, y = np.array([75000., 75000., 75000., 50000., 50000., 50000., 25000., 25000., 25000.]), np.array([1500., 1000., 500., 1500., 1000., 500., 1500., 1000., 500.])

    ftws = ['ftw', 'uc', 'hiacc']
    oops = ['0.0', '1.0e-3', '1.0e-2', '1.0e-1']
    bms = ['0.0', '1.0e-2', '2.0e-2']
    rcs = ['2.0e-3', '2.0e-4', '2.0e-2']

    ftw_names = [b_name_fmt.format(ftw, rcs[0], oops[0], bms[0]) for ftw in ftws]
    oop_names = [b_name_fmt.format(ftws[0], rcs[0], oop, bms[0]) for oop in oops[1:]]

    if poster:
        fmts = ftw_names
    else:
        fmts = ftw_names + oop_names
    folders = ['rstf' for i in ftw_names] + ['stream_ftw' for i in oop_names]
    pretty_fmts = {}
    pretty_fmts_oop = {}
    for fmt in fmts:
        fm_spl = fmt.split('_')
        pretty_fmts[fmt] = r'%s $\lambda$=%1.4g $\dot b$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[2][2:]), float(fm_spl[4][2:]))
        if float(fm_spl[3][3:]) > 0:
            pretty_fmts_oop[fmt] = r'%s+margin, $\dot\varepsilon_{xy}^{(\max)}$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[3][3:]))
        else:
            pretty_fmts_oop[fmt] = r'%s $\dot\varepsilon_{xy}^{(\max)}$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[3][3:]))

    linewidths = [2, 1, 0.5]
    linestyles = ['solid', 'dashed', 'dotted']

    if not debug:
        files = [glob('../stream/{:s}/{:s}_????.vtu'.format(folder, fmt)) for fmt, folder in zip(fmts, folders)]
        inds = [np.argsort(np.array([float(fn[-8:-4]) for fn in filel])) for filel in files]
        fns = [[file_list[ind] for ind in inds_i] for inds_i, file_list in zip(inds, files)]

        times = np.hstack(([0.0, 0.1], np.arange(10.0, 100010.0, 10.0)))

        a_s = {}
        for i, (fmt, folder) in enumerate(zip(fmts, folders)):
            a_s[fmt] = get_vars_or_cache(x, y, fmt, fns[i], folder='../stream/' + folder)
        timess = {fmt: times[:min(a_s[fmt]['fabric 1'].shape[1], len(times))] for fmt, fn in zip(fmts, fns)}

    else:
        timess = {name: np.linspace(0, 4500, 100) for name in fmts}
        vtu = {name: np.ones((9, 100)) for name in ['fabric 1', 'fabric 2', 'fabric 3', 'eigenv 1', 'eigenv 2', 'eigenv 3', 'eigenv 4']}
        a_s = {name: vtu for name in fmts}

    def make_fig(name, dist_num=0, taus=timess):
        fig = plt.figure(constrained_layout=True, figsize=(8, 3.5))
        if poster:
            leg = 0.3
        else:
            leg = 0.5
        gs = gridspec.GridSpec(nrows=4, ncols=2, hspace=0.14, wspace=0.12, left=0.075, bottom=0.0, top=0.95, right=0.98, height_ratios=(1, 1, 0.1, leg))
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
        ax3 = fig.add_subplot(gs[1, 0], sharex=ax1)
        ax4 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax3)
        ax_leg = fig.add_subplot(gs[3, :])
        ax_leg.axis('off')

        def do_plot(fmt, color, lw, ls, dist_num=0, taus=taus, ax1=ax1, ax2=ax2, ax3=ax3, ax4=ax4, fmts=pretty_fmts_oop):
            ax1.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][0 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][1 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax2.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][1 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][2 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

            ax3.plot(taus[fmt], a_s[fmt]['fabric 2'][1 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax4.plot(taus[fmt], a_s[fmt]['fabric 2'][2 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

        for ax in (ax1, ax2, ax3, ax4):
            ax.axhline(1. / 3., color='k', linestyle='dotted')

        # for ax in (ax7, ax8, ax9):
        #     ax.axhline(1., color='k', linestyle='dotted')

        for fmt, color in zip(ftw_names, colors):
            lw = linewidths[0]
            ls = linestyles[0]
            do_plot(fmt, color, lw, ls, taus={k: v / 1000.0 for k, v in taus.items()}, dist_num=dist_num)

        if not poster:
            for fmt, color in zip(oop_names[1:], oop_colors[1:]):
                lw = linewidths[0]
                ls = linestyles[0]
                do_plot(fmt, color, lw, ls, taus={k: v / 1000.0 for k, v in taus.items()}, dist_num=dist_num, fmts=pretty_fmts_oop)

        ax1.set_ylabel(r'$a^{(2)}_{yy}$', fontsize=fs)
        ax3.set_ylabel(r'$a^{(2)}_{zz}$', fontsize=fs)
        # ax7.set_ylabel(r'Woodcock $k$', fontsize=fs)
        ax1.set_title('Depth = {:d} m'.format(int(2000.0 - y[1])), fontsize=fs)
        ax2.set_title('Depth = {:d} m'.format(int(2000.0 - y[2])), fontsize=fs)

        ax1.set_ylim(0, 1.0)
        ax1.set_xlim(0, 25)
        ax3.set_xlabel('Time (kyr)', fontsize=fs)
        ax4.set_xlabel('Time (kyr)', fontsize=fs)
        ax3.set_ylim(0, 1)

        handles, labels = ax4.get_legend_handles_labels()
        ax_leg.legend(handles, labels, loc='upper left', frameon=False, ncol=3, fontsize=fs)

        for ax in (ax1, ax2):
            ax.xaxis.set_tick_params(labelbottom=False)

        for ax in (ax2, ax4):
            ax.yaxis.set_tick_params(labelleft=False)

        for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4)):
            ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
            ax.tick_params(axis='both', which='major', labelsize=fs)

        fig.tight_layout(pad=0.1)
        fig.savefig(name, dpi=300, transparent=poster)
    if poster:
        make_fig('../plots/poster_idealized_core_75_time_backward.png', dist_num=0, taus=timess)
    else:
        make_fig('../plots/idealized_core_75_time_backward.pdf', dist_num=0, taus=timess)


def forward_supp():
    x, y = np.array([75000., 75000., 75000., 50000., 50000., 50000., 25000., 25000., 25000.]), np.array([1500., 1000., 500., 1500., 1000., 500., 1500., 1000., 500.])

    ftws = ['ftw', 'uc', 'hiacc']
    oops = ['0.0', '1.0e-3', '1.0e-2', '1.0e-1']
    bms = ['0.0', '1.0e-2', '2.0e-2']
    rcs = ['2.0e-3', '2.0e-4', '2.0e-2', '2.0e-1']

    ftw_names = [name_fmt.format(ftw, rcs[0], oops[0], bms[0]) for ftw in ftws]
    bm_names = [name_fmt.format(ftws[0], rcs[0], oops[0], bm) for bm in bms[1:]]
    rc_names = [name_fmt.format(ftws[0], rc, oops[0], bms[0]) for rc in rcs[1:]]
    oop_names = [name_fmt.format(ftws[0], rcs[0], oop, bms[0]) for oop in oops[1:]]

    fmts = ftw_names + bm_names + rc_names + oop_names
    pretty_fmts = {}
    pretty_fmts_oop = {}
    for fmt in fmts:
        fm_spl = fmt.split('_')
        pretty_fmts[fmt] = r'%s $\lambda$=%1.4g $\dot b$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[2][2:]), float(fm_spl[4][2:]))
        pretty_fmts_oop[fmt] = r'%s $\dot\varepsilon_{xy}^{(2)}$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[3][3:]))

    linewidths = [2, 1.5, 1, 0.5]
    linestyles = ['solid', 'dashed', 'dotted']

    if not debug:
        files = [glob('../stream/rstf/' + fmt + '_????.vtu') for fmt in fmts]
        inds = [np.argsort(np.array([float(fn[-8:-4]) for fn in filel])) for filel in files]
        fns = [[file_list[ind] for ind in inds_i] for inds_i, file_list in zip(inds, files)]

        times = np.hstack(([0.0, 0.1], np.arange(10.0, 100010.0, 10.0)))

        a_s = {}
        for i, fmt in enumerate(fmts):
            a_s[fmt] = get_vars_or_cache(x, y, fmt, fns[i], folder='../stream/rstf')
        timess = {fmt: times[:min(a_s[fmt]['fabric 1'].shape[1], len(times))] for fmt, fn in zip(fmts, fns)}
        taus = {name: val / (75000.0 / 50.0) for name, val in timess.items()}
        taus25 = {name: val / (25000. / 50.0) for name, val in timess.items()}
        taus50 = {name: val / (50000. / 50.0) for name, val in timess.items()}

    else:
        timess = {name: np.linspace(0, 4500, 100) for name in fmts}
        taus = {name: np.linspace(0, 3, 100) for name in fmts}
        vtu = {name: np.ones((9, 100)) for name in ['fabric 1', 'fabric 2', 'fabric 3', 'eigenv 1', 'eigenv 2', 'eigenv 3', 'eigenv 4']}
        a_s = {name: vtu for name in fmts}

    def make_fig(name, dist_num=0, taus=taus):
        fig = plt.figure(constrained_layout=True, figsize=(8, 5.5))
        gs = gridspec.GridSpec(nrows=5, ncols=3, hspace=0.14, wspace=0.12, left=0.075, bottom=0.0, top=0.965, right=0.98, height_ratios=(1, 1, 1, 0.06, 0.8))
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
        ax3 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
        ax4 = fig.add_subplot(gs[1, 0], sharex=ax1)
        ax5 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax4)
        ax6 = fig.add_subplot(gs[1, 2], sharex=ax1, sharey=ax4)
        ax7 = fig.add_subplot(gs[2, 0], sharex=ax1)
        ax8 = fig.add_subplot(gs[2, 1], sharex=ax1, sharey=ax7)
        ax9 = fig.add_subplot(gs[2, 2], sharex=ax1, sharey=ax7)
        ax_leg = fig.add_subplot(gs[4, :])
        ax_leg.axis('off')

        def do_plot(fmt, color, lw, ls, dist_num=0, taus=taus, ax1=ax1, ax2=ax2, ax3=ax3, ax4=ax4, ax5=ax5, ax6=ax6, ax7=ax7, ax8=ax8, ax9=ax9, fmts=pretty_fmts):
            ax1.plot(taus[fmt], a_s[fmt]['fabric 1'][0 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax2.plot(taus[fmt], a_s[fmt]['fabric 1'][1 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax3.plot(taus[fmt], a_s[fmt]['fabric 1'][2 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

            ax4.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][0 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][0 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax5.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][1 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][1 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax6.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][2 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][2 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

            ax7.plot(taus[fmt], a_s[fmt]['fabric 2'][0 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax8.plot(taus[fmt], a_s[fmt]['fabric 2'][1 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax9.plot(taus[fmt], a_s[fmt]['fabric 2'][2 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

            # ax7.plot(taus[fmt], np.abs(np.log10(np.abs(a_s[fmt]['eigenv 3'][0 + 3 * dist_num, :] / a_s[fmt]['eigenv 2'][0, :])) / np.log10(np.abs(a_s[fmt]['eigenv 2'][0, :] / a_s[fmt]['eigenv 1'][0, :]))), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            # ax8.plot(taus[fmt], np.abs(np.log10(np.abs(a_s[fmt]['eigenv 3'][1 + 3 * dist_num, :] / a_s[fmt]['eigenv 2'][1, :])) / np.log10(np.abs(a_s[fmt]['eigenv 2'][1, :] / a_s[fmt]['eigenv 1'][1, :]))), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            # ax9.plot(taus[fmt], np.abs(np.log10(np.abs(a_s[fmt]['eigenv 3'][2 + 3 * dist_num, :] / a_s[fmt]['eigenv 2'][2, :])) / np.log10(np.abs(a_s[fmt]['eigenv 2'][2, :] / a_s[fmt]['eigenv 1'][2, :]))), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

        for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9):
            ax.axhline(1. / 3., color='k', linestyle='dotted')

        # for ax in (ax7, ax8, ax9):
        #     ax.axhline(1., color='k', linestyle='dotted')

        for fmt, color in zip(oop_names[0:], oop_colors[0:]):
            lw = linewidths[0]
            ls = linestyles[0]
            do_plot(fmt, color, lw, ls, taus=taus, dist_num=dist_num, fmts=pretty_fmts_oop)

        for fmt, color in zip(ftw_names, colors):
            lw = linewidths[0]
            ls = linestyles[0]
            do_plot(fmt, color, lw, ls, taus=taus, dist_num=dist_num)

        for fmt, lw in zip(rc_names, linewidths[1:]):
            color = colors[0]
            ls = linestyles[0]
            do_plot(fmt, color, lw, ls, taus=taus, dist_num=dist_num)

        for fmt, ls in zip(bm_names, linestyles[1:]):
            color = colors[0]
            lw = linewidths[0]
            do_plot(fmt, color, lw, ls, taus=taus, dist_num=dist_num)


        ax1.set_ylabel(r'$a^{(2)}_{xx}$', fontsize=fs)
        ax4.set_ylabel(r'$a^{(2)}_{yy}$', fontsize=fs)
        ax7.set_ylabel(r'$a^{(2)}_{zz}$', fontsize=fs)
        # ax7.set_ylabel(r'Woodcock $k$', fontsize=fs)
        ax1.set_title('Depth = {:d} m'.format(int(2000.0 - y[0])), fontsize=fs)
        ax2.set_title('Depth = {:d} m'.format(int(2000.0 - y[1])), fontsize=fs)
        ax3.set_title('Depth = {:d} m'.format(int(2000.0 - y[2])), fontsize=fs)
        ax1.set_ylim(0, 0.5)
        if np.max(taus[fmts[0]]) < 10:
            ax1.set_xlim(0, 3)
            ax8.set_xlabel(r'Time / $\tau_s$', fontsize=fs)
        else:
            ax1.set_xlim(0, 10000)
            ax8.set_xlabel('Time (yrs)', fontsize=fs)
        ax4.set_ylim(0, 1)
        ax7.set_ylim(0, 1)
        # ax7.set_ylim(0, 2)

        handles, labels = ax7.get_legend_handles_labels()
        handles = handles[3:] + handles[:3]
        labels = labels[3:] + labels[:3]
        ax_leg.legend(handles, labels, loc='upper left', frameon=False, ncol=3, fontsize=fs)

        for ax in (ax1, ax2, ax3, ax4, ax5, ax6):
            ax.xaxis.set_tick_params(labelbottom=False)

        for ax in (ax2, ax3, ax5, ax6, ax8, ax9):
            ax.yaxis.set_tick_params(labelleft=False)

        for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)):
            ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
            ax.tick_params(axis='both', which='major', labelsize=fs)

        fig.tight_layout(pad=0.1)
        fig.savefig(name, dpi=300)
    # make_fig('../plots/idealized_core_75_tau.pdf', dist_num=0, taus=taus)
    make_fig('../plots/idealized_core_75_time_s.pdf', dist_num=0, taus=timess)
    # make_fig('../plots/idealized_core_25_tau.pdf', dist_num=2, taus=taus25)
    make_fig('../plots/idealized_core_25_time_s.pdf', dist_num=2, taus=timess)
    # make_fig('../plots/idealized_core_50_tau.pdf', dist_num=1, taus=taus50)
    make_fig('../plots/idealized_core_50_time_s.pdf', dist_num=1, taus=timess)


def backward_supp():
    x, y = np.array([75000., 75000., 75000., 50000., 50000., 50000., 25000., 25000., 25000.]), np.array([1500., 1000., 500., 1500., 1000., 500., 1500., 1000., 500.])

    ftws = ['ftw', 'uc', 'hiacc']
    oops = ['0.0', '1.0e-3', '1.0e-2', '1.0e-1']
    bms = ['0.0', '1.0e-2', '2.0e-2']
    rcs = ['2.0e-3', '2.0e-4', '2.0e-2', '2.0e-1']

    ftw_names = [b_name_fmt.format(ftw, rcs[0], oops[0], bms[0]) for ftw in ftws]
    bm_names = [b_name_fmt.format(ftws[0], rcs[0], oops[0], bm) for bm in bms[1:]]
    rc_names = [b_name_fmt.format(ftws[0], rc, oops[0], bms[0]) for rc in rcs[1:]]
    oop_names = [b_name_fmt.format(ftws[0], rcs[0], oop, bms[0]) for oop in oops[1:]]

    fmts = ftw_names + bm_names + rc_names
    pretty_fmts = {}
    pretty_fmts_oop = {}
    for fmt in fmts + oop_names:
        fm_spl = fmt.split('_')
        pretty_fmts[fmt] = r'%s $\lambda$=%1.4g $\dot b$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[2][2:]), float(fm_spl[4][2:]))
        pretty_fmts_oop[fmt] = r'%s $\dot\varepsilon_{xy}^{(2)}$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[3][3:]))

    linewidths = [2, 1.5, 1, 0.5]
    linestyles = ['solid', 'dashed', 'dotted']

    if not debug:
        files = [glob('../stream/rstf/' + fmt + '_????.vtu') for fmt in fmts]
        inds = [np.argsort(np.array([float(fn[-8:-4]) for fn in filel])) for filel in files]
        fns = [[file_list[ind] for ind in inds_i] for inds_i, file_list in zip(inds, files)]

        times = np.hstack(([0.0, 0.1], np.arange(10.0, 100010.0, 10.0)))

        a_s = {}
        for i, fmt in enumerate(fmts):
            a_s[fmt] = get_vars_or_cache(x, y, fmt, fns[i], folder='../stream/rstf')

        files_oop = [glob('../stream/stream_ftw/' + fmt + '_????.vtu') for fmt in oop_names]
        inds_oop = [np.argsort(np.array([float(fn[-8:-4]) for fn in filel])) for filel in files_oop]
        fns_oop = [[file_list[ind] for ind in inds_i] for inds_i, file_list in zip(inds_oop, files_oop)]
        for i, fmt in enumerate(oop_names):
            a_s[fmt] = get_vars_or_cache(x, y, fmt, fns_oop[i], folder='../stream/stream_ftw')

        timess = {fmt: times[:min(a_s[fmt]['fabric 1'].shape[1], len(times))] for fmt in fmts + oop_names}
        taus = {name: val / (75000.0 / 50.0) for name, val in timess.items()}
        taus25 = {name: val / (25000. / 50.0) for name, val in timess.items()}
        taus50 = {name: val / (50000. / 50.0) for name, val in timess.items()}

    else:
        timess = {name: np.linspace(0, 4500, 100) for name in fmts}
        taus = {name: np.linspace(0, 3, 100) for name in fmts}
        vtu = {name: np.ones((9, 100)) for name in ['fabric 1', 'fabric 2', 'fabric 3', 'eigenv 1', 'eigenv 2', 'eigenv 3', 'eigenv 4']}
        a_s = {name: vtu for name in fmts}

    def make_fig(name, dist_num=0, taus=taus):
        fig = plt.figure(constrained_layout=True, figsize=(8, 5.5))
        gs = gridspec.GridSpec(nrows=5, ncols=3, hspace=0.14, wspace=0.12, left=0.075, bottom=0.0, top=0.965, right=0.98, height_ratios=(1, 1, 1, 0.06, 0.8))
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
        ax3 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
        ax4 = fig.add_subplot(gs[1, 0], sharex=ax1)
        ax5 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax4)
        ax6 = fig.add_subplot(gs[1, 2], sharex=ax1, sharey=ax4)
        ax7 = fig.add_subplot(gs[2, 0], sharex=ax1)
        ax8 = fig.add_subplot(gs[2, 1], sharex=ax1, sharey=ax7)
        ax9 = fig.add_subplot(gs[2, 2], sharex=ax1, sharey=ax7)
        ax_leg = fig.add_subplot(gs[4, :])
        ax_leg.axis('off')

        def do_plot(fmt, color, lw, ls, dist_num=0, taus=taus, ax1=ax1, ax2=ax2, ax3=ax3, ax4=ax4, ax5=ax5, ax6=ax6, ax7=ax7, ax8=ax8, ax9=ax9, fmts=pretty_fmts):
            ax1.plot(taus[fmt], a_s[fmt]['fabric 1'][0 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax2.plot(taus[fmt], a_s[fmt]['fabric 1'][1 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax3.plot(taus[fmt], a_s[fmt]['fabric 1'][2 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

            ax4.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][0 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][0 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax5.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][1 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][1 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax6.plot(taus[fmt], np.abs(1.0 - a_s[fmt]['fabric 1'][2 + 3 * dist_num, :] - a_s[fmt]['fabric 2'][2 + 3 * dist_num, :]), color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])

            ax7.plot(taus[fmt], a_s[fmt]['fabric 2'][0 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax8.plot(taus[fmt], a_s[fmt]['fabric 2'][1 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])
            ax9.plot(taus[fmt], a_s[fmt]['fabric 2'][2 + 3 * dist_num, :], color=color, linewidth=lw, linestyle=ls, label=fmts[fmt])


        for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9):
            ax.axhline(1. / 3., color='k', linestyle='dotted')

        # for ax in (ax7, ax8, ax9):
        #     ax.axhline(1., color='k', linestyle='dotted')

        for fmt, color in zip(oop_names[0:], oop_colors[0:]):
            lw = linewidths[0]
            ls = linestyles[0]
            do_plot(fmt, color, lw, ls, taus=taus, dist_num=dist_num, fmts=pretty_fmts_oop)

        for fmt, color in zip(ftw_names, colors):
            lw = linewidths[0]
            ls = linestyles[0]
            do_plot(fmt, color, lw, ls, taus=taus, dist_num=dist_num)

        for fmt, lw in zip(rc_names, linewidths[1:]):
            color = colors[0]
            ls = linestyles[0]
            do_plot(fmt, color, lw, ls, taus=taus, dist_num=dist_num)

        for fmt, ls in zip(bm_names, linestyles[1:]):
            color = colors[0]
            lw = linewidths[0]
            do_plot(fmt, color, lw, ls, taus=taus, dist_num=dist_num)


        ax1.set_ylabel(r'$a^{(2)}_{xx}$', fontsize=fs)
        ax4.set_ylabel(r'$a^{(2)}_{yy}$', fontsize=fs)
        ax7.set_ylabel(r'$a^{(2)}_{zz}$', fontsize=fs)
        # ax7.set_ylabel(r'Woodcock $k$', fontsize=fs)
        ax1.set_title('Depth = {:d} m'.format(int(2000.0 - y[0])), fontsize=fs)
        ax2.set_title('Depth = {:d} m'.format(int(2000.0 - y[1])), fontsize=fs)
        ax3.set_title('Depth = {:d} m'.format(int(2000.0 - y[2])), fontsize=fs)
        ax1.set_ylim(0, 0.5)
        if False:  # np.max(taus[fmts[0]]) < 10:
            ax1.set_xlim(0, 3)
            ax8.set_xlabel(r'Time / $\tau_s$', fontsize=fs)
        else:
            ax1.set_xlim(0, 25)
            ax8.set_xlabel('Time (kyr)', fontsize=fs)
        ax4.set_ylim(0, 1)
        ax7.set_ylim(0, 1)
        # ax7.set_ylim(0, 2)

        handles, labels = ax7.get_legend_handles_labels()
        handles = handles[3:] + handles[:3]
        labels = labels[3:] + labels[:3]
        ax_leg.legend(handles, labels, loc='upper left', frameon=False, ncol=3, fontsize=fs)

        for ax in (ax1, ax2, ax3, ax4, ax5, ax6):
            ax.xaxis.set_tick_params(labelbottom=False)

        for ax in (ax2, ax3, ax5, ax6, ax8, ax9):
            ax.yaxis.set_tick_params(labelleft=False)

        for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9)):
            ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
            ax.tick_params(axis='both', which='major', labelsize=fs)

        fig.tight_layout(pad=0.1)
        fig.savefig(name, dpi=300)

    tms = {k: v / 1000.0 for k, v in timess.items()}
    # make_fig('../plots/idealized_core_75_tau_backward.pdf', dist_num=0, taus=taus)
    make_fig('../plots/idealized_core_75_time_backward_s.pdf', dist_num=0, taus=tms)
    # make_fig('../plots/idealized_core_25_tau_backward.pdf', dist_num=2, taus=taus25)
    make_fig('../plots/idealized_core_25_time_backward_s.pdf', dist_num=2, taus=tms)
    # make_fig('../plots/idealized_core_50_tau_backward.pdf', dist_num=1, taus=taus50)
    make_fig('../plots/idealized_core_50_time_backward_s.pdf', dist_num=1, taus=tms)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-slow', action='store_true')
    parser.add_argument('-forward', action='store_true')
    parser.add_argument('-backward', action='store_true')
    parser.add_argument('-oop', action='store_true')
    parser.add_argument('-supp', action='store_true')
    parser.add_argument('-poster', action='store_true')
    args = parser.parse_args()
    if args.slow:
        slow()
    if args.backward:
        backward(poster=args.poster)
    if args.forward:
        forward(poster=args.poster)
    if args.supp:
        forward_supp()
        backward_supp()
    if args.oop:
        oop(poster=args.poster)
