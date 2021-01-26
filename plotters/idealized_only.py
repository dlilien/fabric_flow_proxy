#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the MIT license.

"""

"""
from joblib import Parallel, delayed
import matplotlib.gridspec as gridspec
import numpy as np
from lib import anisolib, fabricplotlib
from modeltools.lib import fastvtulib
from matplotlib.tri import Triangulation
import matplotlib.pyplot as plt
from glob import glob
from matplotlib import rc as mplrc
mplrc('font', **{'family': 'sans-serif', 'sans-serif': ['Linux Biolinum O']})

# import os

# from mpl_toolkits.axes_grid1.inset_locator import inset_axes

fs = 9
bfs = 12

name_fmt = 's_{:s}_rc{:s}_oop{:s}_bm{:s}'
name_dict = {'ftw': 'Converge\nupstream', 'uc': 'Converge\neverywhere', 'hiacc': 'High\naccumulation'}


def woodcock(stuff_dict):
    norm = np.abs(stuff_dict['eigenv 3']) + np.abs(stuff_dict['eigenv 2']) + np.abs(stuff_dict['eigenv 1'])
    e1 = np.abs(stuff_dict['eigenv 1']) / norm
    e2 = np.abs(stuff_dict['eigenv 2']) / norm
    e3 = np.abs(stuff_dict['eigenv 3']) / norm
    out = np.log10(e3 / e2) / np.log10(e2 / e1)
    out[e1 == 0] = 1e8
    out[e2 == 0] = 1e8
    return out


def main():
    idealized_with_cartoon()
    comp_convergence()
    return
    idealized_init_zoom()
    idealized_init()
    comp_oop()
    comp_bm()
    comp_recryst()


def get_vtus(fmts, folder='rstf', sync=False):
    files = [glob('../stream/' + folder + '/' + fmt + '_????.vtu') for fmt in fmts]
    inds = [np.argsort(np.array([float(fn[-8:-4]) for fn in filel])) for filel in files]
    if not sync:
        print('Returning files:')
        print([filel[int(ind[-1])] for ind, filel in zip(inds, files)])
        vtus = [fastvtulib.get_structured_vtu(filel[int(ind[-1])]) for ind, filel in zip(inds, files)]
        maxtime = None
        return vtus, maxtime
    else:
        def get_times(fmt):
            times = []
            out_file = '../stream/' + folder + '/' + fmt + '.result'
            with open(out_file) as fout:
                for line in fout:
                    if line[:4] == 'Time':
                        times.append(float(line.rstrip('\n').split()[-1]))
            times = np.array(times)
            return np.unique(times - times[0])

        n_jobs = len(fmts)
        timess = Parallel(n_jobs=n_jobs)(delayed(get_times)(fmt) for fmt in fmts)
        maxtime = np.min([np.max(times) for times in timess])
        inds_times = [int(np.where(times == maxtime)[0][0]) for times in timess]
        vtus = [fastvtulib.get_structured_vtu(filel[int(ind[ind_t])]) for ind, ind_t, filel in zip(inds, inds_times, files)]

        for vtu in vtus:
            vtu.rawdata_dict['vel'] = np.sqrt(vtu.rawdata_dict['aiflow'][:, 0]**2. + vtu.rawdata_dict['aiflow'][:, 1]**2.)
        return vtus, maxtime


def comp_convergence(convergences=['uc', 'ftw', 'hiacc']):
    # We want to find the filename for the highest index where we have all 3
    fmts = [name_fmt.format(conv, '2.0e-3', '0.0', '0.0') for conv in convergences]
    vtus, time = get_vtus(fmts)

    tris = [Triangulation(np.array([rc[0] + 100000. for rc in vtu.raw_coords[:, 0]]) / 1000.,
                          np.array([rc[1] for rc in vtu.raw_coords[:, 0]]), vtu.simptt) for vtu in vtus]

    gs = gridspec.GridSpec(
        2,
        len(convergences) + 1,
        height_ratios=(
            1.,
            1.),
        width_ratios=(
            1,
            1,
            1,
            0.05),
        hspace=0.10,
        wspace=0.1,
        left=0.1,
        bottom=0.2,
        top=0.9,
        right=0.9)

    fig = plt.figure(figsize=(7.8, 4.0))
    a12_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    angles_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    # if len(axes) % 2 == 0:
    #     ax_c = fig.add_subplot(gs[1, len(axes) // 2 - 1:len(axes) // 2])
    # else:
    #     ax_c = fig.add_subplot(gs[1, (len(axes) - 1) // 2])
    ax_c_a12 = fig.add_subplot(gs[0, -1])
    ax_c_angles = fig.add_subplot(gs[1, -1])

    for axa, axang, tri, vtu, conv in zip(a12_axes, angles_axes, tris, vtus, convergences):

        cm1 = axa.tricontourf(tri, vtu.rawdata_dict['aiflow'][:, 0],
                              cmap='Reds', levels=np.linspace(0.0, 50, 101), extend='max')
        levels = [1000, 5000, 10000, 50000, 100000]
        CS = axa.tricontour(tri, vtu.rawdata_dict['age'], colors='k', levels=levels, linewidths=0.5)
        fmt = {level: '{:d} kyr'.format(level // 1000) for level in levels}
        axa.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=8)
        for c in cm1.collections:
            c.set_edgecolor("face")

        axa.set_title('{:s}'.format(name_dict[conv]), fontsize=fs)
        axa.set_xlim(0, 175)
        axa.set_xticks([0., 50., 100., 150.])
        axa.set_xticklabels(['' for tick in axa.get_xticklabels()])

        cm2 = axang.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 3'],
            cmap='summer',
            levels=np.linspace(
                0.333333,
                1,
                101),
            extend='neither')
        for c in cm2.collections:
            c.set_edgecolor("face")
        fabricplotlib.quiver(axang, vtu, scale=15, width=0.005)

    cbr = plt.colorbar(cm1, cax=ax_c_a12, orientation='vertical', ticks=(0, 25., 50.))
    cbr.set_label(label=r'$v_x$ (m a$^{-1})$', size=fs)
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)

    cbr = plt.colorbar(cm2, cax=ax_c_angles, orientation='vertical', ticks=(1. / 3., 2. / 3., 1.))
    cbr.set_label(label=r'$a^{(2)}_{1}$', size=fs)
    cbr.ax.set_yticklabels([r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)

    angles_axes[0].set_ylabel('Height (m)', fontsize=fs)
    for ax in a12_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in angles_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    angles_axes[1].set_xlabel('Distance (km)', fontsize=fs)

    angles_axes[0].scatter(-1000, -1000, marker=r'$\uparrow$', label='Single max. in x-z', color='k')
    angles_axes[0].legend(loc='upper left', bbox_to_anchor=(0.2, -0.3), ncol=2, fontsize=fs, framealpha=1.0)

    for letter, ax in zip('abcdefghijkl', a12_axes + angles_axes):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.set_ylim(0, 2250)
        ax.set_xlim(0, 175)
        ax.tick_params(axis='both', which='major', labelsize=fs)

    fig.savefig('../plots/idealized_convergence.png', dpi=300)
    fig.savefig('../plots/poster_idealized_convergence.png', dpi=300, transparent=True)

    fig = plt.figure(figsize=(7.5, 5.0))
    a12_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    angles_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    # if len(axes) % 2 == 0:
    #     ax_c = fig.add_subplot(gs[1, len(axes) // 2 - 1:len(axes) // 2])
    # else:
    #     ax_c = fig.add_subplot(gs[1, (len(axes) - 1) // 2])
    # ax_c_a12 = fig.add_subplot(gs[0, -1])
    ax_c_angles = fig.add_subplot(gs[1, -1])
    # ax_c_woodcock = fig.add_subplot(gs[2, -1])

    for axa, axang, tri, vtu, conv in zip(a12_axes, angles_axes, tris, vtus, convergences):
        cm1 = axa.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 1'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm1.collections:
            c.set_edgecolor("face")
        axa.set_title('{:s}'.format(name_dict[conv]), fontsize=fs)
        axa.set_xlim(0, 175)
        axa.set_xticks([0., 50., 100., 150.])
        axa.set_xticklabels(['' for tick in axa.get_xticklabels()])

        cm2 = axang.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 2'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm2.collections:
            c.set_edgecolor("face")
        axang.set_xlim(0, 175)
        axang.set_xticks([0., 50., 100., 150.])
        axang.set_xticklabels(['' for tick in axang.get_xticklabels()])
    cbr = plt.colorbar(
        cm2,
        cax=ax_c_angles,
        orientation='vertical',
        ticks=(
            0,
            1. / 3.,
            2. / 3.,
            1.0),
        label=r'$a^{(2)}_{zz}$')
    cbr.ax.set_yticklabels(['0', r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)

    a12_axes[0].set_ylabel('Height (m)')
    angles_axes[0].set_ylabel('Height (m)')
    for ax in a12_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in angles_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    angles_axes[1].set_xlabel('Distance (km)', fontsize=fs)


    for letter, ax in zip('abcdefghijkl', a12_axes + angles_axes):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.set_ylim(0, 2250)
        ax.tick_params(axis='both', which='major', labelsize=fs)
    fig.savefig('../plots/idealized_converge_eigs.png', dpi=300)
    fig.savefig('../plots/poster_idealized_converge_eigs.png', dpi=300, Transparent=True)


def comp_recryst(recrysts=['2.0e-2', '2.0e-3', '2.0e-4']):
    conv = 'ftw'
    bm = '0.0'
    oop = '0.0'
    # We want to find the filename for the highest index where we have all 3
    fmts = [name_fmt.format(conv, rc, oop, bm) for rc in recrysts]
    vtus, time = get_vtus(fmts)

    tris = [Triangulation(np.array([rc[0] + 100000. for rc in vtu.raw_coords[:, 0]]) / 1000.,
                          np.array([rc[1] for rc in vtu.raw_coords[:, 0]]), vtu.simptt) for vtu in vtus]

    gs = gridspec.GridSpec(3, 3 + 1, height_ratios=(1., 1., 1.), width_ratios=(1, 1, 1, 0.075),
                           hspace=0.10, wspace=0.1, left=0.1, bottom=0.1, top=0.9, right=0.9)

    fig = plt.figure(figsize=(7.5, 5.0))
    a12_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    angles_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    woodcock_axes = [fig.add_subplot(gs[2, i]) for i in range(3)]
    # if len(axes) % 2 == 0:
    #     ax_c = fig.add_subplot(gs[1, len(axes) // 2 - 1:len(axes) // 2])
    # else:
    #     ax_c = fig.add_subplot(gs[1, (len(axes) - 1) // 2])
    ax_c_a12 = fig.add_subplot(gs[0, -1])
    ax_c_angles = fig.add_subplot(gs[1, -1])
    ax_c_woodcock = fig.add_subplot(gs[2, -1])

    for axa, axang, axw, tri, vtu, rc in zip(a12_axes, angles_axes, woodcock_axes, tris, vtus, recrysts):

        cm1 = axa.tricontourf(tri, vtu.rawdata_dict['aiflow'][:, 0],
                              cmap='Reds', levels=np.linspace(0.0, 50, 101), extend='max')
        levels = [1000, 5000, 10000, 50000, 100000]
        CS = axa.tricontour(tri, vtu.rawdata_dict['age'], colors='k', levels=levels, linewidths=0.5)
        fmt = {level: '{:d} kyr'.format(level // 1000) for level in levels}
        axa.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=8)
        for c in cm1.collections:
            c.set_edgecolor("face")

        axa.set_title('{:s}'.format(rc))
        axa.set_xlim(0, 175)
        axa.set_xticks([0., 50., 100., 150.])
        axa.set_xticklabels(['' for tick in axa.get_xticklabels()])

        cm2 = axang.tricontourf(
            tri,
            vtu.rawdata_dict['fabric 2'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm2.collections:
            c.set_edgecolor("face")
        axang.set_xlim(0, 175)
        axang.set_xticks([0., 50., 100., 150.])
        axang.set_xticklabels(['' for tick in axang.get_xticklabels()])

        k = woodcock(vtu.rawdata_dict)
        k[np.logical_or(np.isinf(k), np.isnan(k))] = 1.0e8
        cm3 = axw.tricontourf(tri, k, cmap='PuOr', levels=np.linspace(0, 2, 201), extend='max')
        for c in cm3.collections:
            c.set_edgecolor("face")
        axw.set_xlim(0, 175)
        axw.set_xticks([0., 50., 100., 150.])

    cbr = plt.colorbar(cm1, cax=ax_c_a12, orientation='vertical', label=r'$v_x$ (m a$^{-1})$', ticks=(0, 25., 50.))
    cbr = plt.colorbar(
        cm2,
        cax=ax_c_angles,
        orientation='vertical',
        ticks=(
            0,
            1. / 3.,
            2. / 3.,
            1.),
        label=r'$a^{(2)}_{zz}$')
    cbr.ax.set_yticklabels(['0', r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])
    cbr = plt.colorbar(cm3, cax=ax_c_woodcock, orientation='vertical', ticks=(0, 1, 2), label='Woodcock k')

    angles_axes[0].set_ylabel('Height (m)')
    for ax in a12_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in angles_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in woodcock_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    woodcock_axes[1].set_xlabel('Distance (km)')

    for letter, ax in zip('abcdefghijkl', a12_axes + angles_axes + woodcock_axes):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.set_ylim(0, 2250)

    fig.savefig('../plots/idealized_recryst.png', dpi=300)

    fig = plt.figure(figsize=(7.5, 5.0))
    a12_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    angles_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    woodcock_axes = [fig.add_subplot(gs[2, i]) for i in range(3)]
    # if len(axes) % 2 == 0:
    #     ax_c = fig.add_subplot(gs[1, len(axes) // 2 - 1:len(axes) // 2])
    # else:
    #     ax_c = fig.add_subplot(gs[1, (len(axes) - 1) // 2])
    # ax_c_a12 = fig.add_subplot(gs[0, -1])
    ax_c_angles = fig.add_subplot(gs[1, -1])
    # ax_c_woodcock = fig.add_subplot(gs[2, -1])

    for axa, axang, axw, tri, vtu, rc in zip(a12_axes, angles_axes, woodcock_axes, tris, vtus, recrysts):
        cm1 = axa.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 1'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm1.collections:
            c.set_edgecolor("face")
        axa.set_title('{:s}'.format(rc))
        axa.set_xlim(0, 175)
        axa.set_xticks([0., 50., 100., 150.])
        axa.set_xticklabels(['' for tick in axa.get_xticklabels()])

        cm2 = axang.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 2'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm2.collections:
            c.set_edgecolor("face")
        axang.set_xlim(0, 175)
        axang.set_xticks([0., 50., 100., 150.])
        axang.set_xticklabels(['' for tick in axang.get_xticklabels()])

        cm3 = axw.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 3'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm3.collections:
            c.set_edgecolor("face")
        axw.set_xlim(0, 175)
        axw.set_xticks([0., 50., 100., 150.])

    cbr = plt.colorbar(
        cm2,
        cax=ax_c_angles,
        orientation='vertical',
        ticks=(
            0,
            1. / 3.,
            2. / 3.,
            1.0),
        label=r'$a^{(2)}_{zz}$')
    cbr.ax.set_yticklabels(['0', r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])

    angles_axes[0].set_ylabel('Height (m)')
    for ax in a12_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in angles_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in woodcock_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    woodcock_axes[1].set_xlabel('Distance (km)')

    for letter, ax in zip('abcdefghijkl', a12_axes + angles_axes + woodcock_axes):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.set_ylim(0, 2250)
    fig.savefig('../plots/idealized_recryst_eigs.png', dpi=300)


def comp_oop(oops=['1.0e-3', '1.0e-2', '1.0e-1']):
    conv = 'ftw'
    bm = '0.0'
    rc = '2.0e-3'
    # We want to find the filename for the highest index where we have all 3
    fmts = [name_fmt.format(conv, rc, oop, bm) for oop in oops]

    vtus, time = get_vtus(fmts, folder='stream_ftw')

    tris = [Triangulation(np.array([rc[0] + 100000. for rc in vtu.raw_coords[:, 0]]) / 1000.,
                          np.array([rc[1] for rc in vtu.raw_coords[:, 0]]), vtu.simptt) for vtu in vtus]

    gs = gridspec.GridSpec(3, 3 + 1, height_ratios=(1., 1., 1.), width_ratios=(1, 1, 1, 0.075),
                           hspace=0.10, wspace=0.1, left=0.1, bottom=0.1, top=0.9, right=0.9)

    fig = plt.figure(figsize=(7.5, 5.0))
    a12_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    angles_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    woodcock_axes = [fig.add_subplot(gs[2, i]) for i in range(3)]
    # if len(axes) % 2 == 0:
    #     ax_c = fig.add_subplot(gs[1, len(axes) // 2 - 1:len(axes) // 2])
    # else:
    #     ax_c = fig.add_subplot(gs[1, (len(axes) - 1) // 2])
    ax_c_a12 = fig.add_subplot(gs[0, -1])
    ax_c_angles = fig.add_subplot(gs[1, -1])
    ax_c_woodcock = fig.add_subplot(gs[2, -1])

    for axa, axang, axw, tri, vtu, oop in zip(a12_axes, angles_axes, woodcock_axes, tris, vtus, oops):
        ei6 = vtu.rawdata_dict['eigenv 6']
        ei6[ei6 < 0] = ei6[ei6 < 0] + 180.
        cm1 = axa.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 6'],
            cmap='twilight_shifted',
            levels=np.linspace(
                0,
                180,
                361),
            extend='neither')

        axa.set_title(r'$\dot\varepsilon_{xy}^{\max}$=%s a$^{-1}$' % oop)
        axa.set_xlim(0, 175)
        axa.set_xticks([0., 50., 100., 150.])
        axa.set_xticklabels(['' for tick in axa.get_xticklabels()])

        cm2 = axang.tricontourf(
            tri,
            vtu.rawdata_dict['fabric 2'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm2.collections:
            c.set_edgecolor("face")
        axang.set_xlim(0, 175)
        axang.set_xticks([0., 50., 100., 150.])
        axang.set_xticklabels(['' for tick in axang.get_xticklabels()])

        k = woodcock(vtu.rawdata_dict)
        cm3 = axw.tricontourf(tri, k, cmap='PuOr', levels=np.linspace(0, 2, 201), extend='max')
        for c in cm3.collections:
            c.set_edgecolor("face")
        axw.set_xlim(0, 175)
        axw.set_xticks([0., 50., 100., 150.])

    cbr = plt.colorbar(cm1, cax=ax_c_a12, orientation='vertical', label=r'$\phi (^o)$', ticks=(0., 180.))
    cbr = plt.colorbar(
        cm2,
        cax=ax_c_angles,
        orientation='vertical',
        ticks=(
            0,
            1. / 3.,
            2. / 3.,
            1.),
        label=r'$a^{(2)}_{zz}$')
    cbr.ax.set_yticklabels(['0', r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])
    cbr = plt.colorbar(cm3, cax=ax_c_woodcock, orientation='vertical', ticks=(0, 1, 2), label='Woodcock k')

    angles_axes[0].set_ylabel('Height (m)')
    for ax in a12_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in angles_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in woodcock_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    woodcock_axes[1].set_xlabel('Distance (km)')

    for letter, ax in zip('abcdefghijkl', a12_axes + angles_axes + woodcock_axes):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.set_ylim(0, 2250)

    fig.savefig('../plots/idealized_oop.png', dpi=300)

    fig = plt.figure(figsize=(7.5, 5.0))
    a12_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    angles_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    woodcock_axes = [fig.add_subplot(gs[2, i]) for i in range(3)]
    # if len(axes) % 2 == 0:
    #     ax_c = fig.add_subplot(gs[1, len(axes) // 2 - 1:len(axes) // 2])
    # else:
    #     ax_c = fig.add_subplot(gs[1, (len(axes) - 1) // 2])
    # ax_c_a12 = fig.add_subplot(gs[0, -1])
    ax_c_angles = fig.add_subplot(gs[1, -1])
    # ax_c_woodcock = fig.add_subplot(gs[2, -1])

    for axa, axang, axw, tri, vtu, oop in zip(a12_axes, angles_axes, woodcock_axes, tris, vtus, oops):
        cm1 = axa.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 1'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm1.collections:
            c.set_edgecolor("face")
        axa.set_title(r'$\dot\varepsilon_{xy}^{\max}$=%s' % oop)
        axa.set_xlim(0, 175)
        axa.set_xticks([0., 50., 100., 150.])
        axa.set_xticklabels(['' for tick in axa.get_xticklabels()])

        cm2 = axang.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 2'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm2.collections:
            c.set_edgecolor("face")
        axang.set_xlim(0, 175)
        axang.set_xticks([0., 50., 100., 150.])
        axang.set_xticklabels(['' for tick in axang.get_xticklabels()])

        cm3 = axw.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 3'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm3.collections:
            c.set_edgecolor("face")
        axw.set_xlim(0, 175)
        axw.set_xticks([0., 50., 100., 150.])

    cbr = plt.colorbar(
        cm2,
        cax=ax_c_angles,
        orientation='vertical',
        ticks=(
            0,
            1. / 3.,
            2. / 3.,
            1.0),
        label=r'$a^{(2)}_{zz}$')
    cbr.ax.set_yticklabels(['0', r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])

    angles_axes[0].set_ylabel('Height (m)')
    for ax in a12_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in angles_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in woodcock_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    woodcock_axes[1].set_xlabel('Distance (km)')

    for letter, ax in zip('abcdefghijkl', a12_axes + angles_axes + woodcock_axes):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.set_ylim(0, 2250)
    fig.savefig('../plots/idealized_oop_eigs.png', dpi=300)


def comp_bm(bms=['0.0', '1.0e-2', '2.0e-2']):
    conv = 'ftw'
    rc = '2.0e-3'
    oop = '0.0'
    # We want to find the filename for the highest index where we have all 3
    fmts = [name_fmt.format(conv, rc, oop, bm) for bm in bms]

    vtus, time = get_vtus(fmts)

    tris = [Triangulation(np.array([rc[0] + 100000. for rc in vtu.raw_coords[:, 0]]) / 1000.,
                          np.array([rc[1] for rc in vtu.raw_coords[:, 0]]), vtu.simptt) for vtu in vtus]

    gs = gridspec.GridSpec(3, 3 + 1, height_ratios=(1., 1., 1.), width_ratios=(1, 1, 1, 0.075),
                           hspace=0.10, wspace=0.1, left=0.1, bottom=0.1, top=0.9, right=0.9)

    fig = plt.figure(figsize=(7.5, 5.0))
    a12_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    angles_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    woodcock_axes = [fig.add_subplot(gs[2, i]) for i in range(3)]
    # if len(axes) % 2 == 0:
    #     ax_c = fig.add_subplot(gs[1, len(axes) // 2 - 1:len(axes) // 2])
    # else:
    #     ax_c = fig.add_subplot(gs[1, (len(axes) - 1) // 2])
    ax_c_a12 = fig.add_subplot(gs[0, -1])
    ax_c_angles = fig.add_subplot(gs[1, -1])
    ax_c_woodcock = fig.add_subplot(gs[2, -1])

    for axa, axang, axw, tri, vtu, bm in zip(a12_axes, angles_axes, woodcock_axes, tris, vtus, bms):

        cm1 = axa.tricontourf(tri, vtu.rawdata_dict['aiflow'][:, 0],
                              cmap='Reds', levels=np.linspace(0.0, 50, 101), extend='max')
        levels = [1000, 5000, 10000, 50000, 100000]
        CS = axa.tricontour(tri, vtu.rawdata_dict['age'], colors='k', levels=levels, linewidths=0.5)
        fmt = {level: '{:d} kyr'.format(level // 1000) for level in levels}
        axa.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=8)
        for c in cm1.collections:
            c.set_edgecolor("face")

        axa.set_title(r'$\dot b$=%s m a$^{-1}$' % bm)
        axa.set_xlim(0, 175)
        axa.set_xticks([0., 50., 100., 150.])
        axa.set_xticklabels(['' for tick in axa.get_xticklabels()])

        cm2 = axang.tricontourf(
            tri,
            vtu.rawdata_dict['fabric 2'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm2.collections:
            c.set_edgecolor("face")
        axang.set_xlim(0, 175)
        axang.set_xticks([0., 50., 100., 150.])
        axang.set_xticklabels(['' for tick in axang.get_xticklabels()])

        k = woodcock(vtu.rawdata_dict)
        cm3 = axw.tricontourf(tri, k, cmap='PuOr', levels=np.linspace(0, 2, 201), extend='max')
        for c in cm3.collections:
            c.set_edgecolor("face")
        axw.set_xlim(0, 175)
        axw.set_xticks([0., 50., 100., 150.])

    cbr = plt.colorbar(cm1, cax=ax_c_a12, orientation='vertical', label=r'$v_x$ (m a$^{-1})$', ticks=(0, 25., 50.))
    cbr = plt.colorbar(
        cm2,
        cax=ax_c_angles,
        orientation='vertical',
        ticks=(
            0,
            1. / 3.,
            2. / 3.,
            1.),
        label=r'$a^{(2)}_{zz}$')
    cbr.ax.set_yticklabels(['0', r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])
    cbr = plt.colorbar(cm3, cax=ax_c_woodcock, orientation='vertical', ticks=(0, 1, 2), label='Woodcock k')

    angles_axes[0].set_ylabel('Height (m)')
    for ax in a12_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in angles_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in woodcock_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    woodcock_axes[1].set_xlabel('Distance (km)')

    for letter, ax in zip('abcdefghijkl', a12_axes + angles_axes + woodcock_axes):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.set_ylim(0, 2250)

    fig.savefig('../plots/idealized_bm.png', dpi=300)

    fig = plt.figure(figsize=(7.5, 5.0))
    a12_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    angles_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]
    woodcock_axes = [fig.add_subplot(gs[2, i]) for i in range(3)]
    # if len(axes) % 2 == 0:
    #     ax_c = fig.add_subplot(gs[1, len(axes) // 2 - 1:len(axes) // 2])
    # else:
    #     ax_c = fig.add_subplot(gs[1, (len(axes) - 1) // 2])
    # ax_c_a12 = fig.add_subplot(gs[0, -1])
    ax_c_angles = fig.add_subplot(gs[1, -1])
    # ax_c_woodcock = fig.add_subplot(gs[2, -1])

    for axa, axang, axw, tri, vtu, bm in zip(a12_axes, angles_axes, woodcock_axes, tris, vtus, bms):
        cm1 = axa.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 1'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm1.collections:
            c.set_edgecolor("face")
        axa.set_title(r'$\dot b$=%s m a$^{-1}$' % bm)
        axa.set_xlim(0, 175)
        axa.set_xticks([0., 50., 100., 150.])
        axa.set_xticklabels(['' for tick in axa.get_xticklabels()])

        cm2 = axang.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 2'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm2.collections:
            c.set_edgecolor("face")
        axang.set_xlim(0, 175)
        axang.set_xticks([0., 50., 100., 150.])
        axang.set_xticklabels(['' for tick in axang.get_xticklabels()])

        cm3 = axw.tricontourf(
            tri,
            vtu.rawdata_dict['eigenv 3'],
            cmap='BrBG',
            norm=anisolib.fab_norm,
            levels=np.linspace(
                0,
                1,
                101),
            extend='neither')
        for c in cm3.collections:
            c.set_edgecolor("face")
        axw.set_xlim(0, 175)
        axw.set_xticks([0., 50., 100., 150.])

    cbr = plt.colorbar(
        cm2,
        cax=ax_c_angles,
        orientation='vertical',
        ticks=(
            0,
            1. / 3.,
            2. / 3.,
            1.0),
        label=r'$a^{(2)}_{zz}$')
    cbr.ax.set_yticklabels(['0', r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])

    angles_axes[0].set_ylabel('Height (m)')
    for ax in a12_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in angles_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    for ax in woodcock_axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    woodcock_axes[1].set_xlabel('Distance (km)')

    for letter, ax in zip('abcdefghijkl', a12_axes + angles_axes + woodcock_axes):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.set_ylim(0, 2250)
    fig.savefig('../plots/idealized_bm_eigs.png', dpi=300)


def idealized_init_zoom():
    # We want to find the filename for the highest index where we have all 3
    files = glob('../stream/rstf/deform2e-5_????.vtu')
    inds = np.argsort(np.array([float(fn[-8:-4]) for fn in files]))
    files = [files[i] for i in inds]

    vtu = fastvtulib.get_structured_vtu(files[-1])
    vtu.rawdata_dict['vel'] = np.sqrt(vtu.rawdata_dict['aiflow'][:, 0]**2. + vtu.rawdata_dict['aiflow'][:, 1]**2.)

    tri = Triangulation(np.array([rc[0] + 100000. for rc in vtu.raw_coords[:, 0]]) / 1000.,
                        np.array([rc[1] for rc in vtu.raw_coords[:, 0]]), vtu.simptt)
    gs = gridspec.GridSpec(
        2,
        3,
        width_ratios=(
            1.,
            1.,
            1.),
        height_ratios=(
            1,
            0.1),
        hspace=0.55,
        wspace=0.1,
        left=0.1,
        bottom=0.2,
        top=0.9,
        right=0.98)

    fig = plt.figure(figsize=(7.5, 3.0))
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    cbr_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]

    cm1 = axes[0].tricontourf(tri, vtu.rawdata_dict['strainrate 4'], cmap='bwr',
                              levels=np.linspace(-2.5e-4, 2.5e-4, 101), extend='both')
    levels = [1000, 5000, 10000, 50000, 100000]
    axes[0].set_xlim(0, 20)
    axes[0].set_ylim(0, 2250)
    CS = axes[0].tricontour(tri, vtu.rawdata_dict['age'], colors='k', levels=levels, linewidths=0.5)
    fmt = {level: '{:d} kyr'.format(level // 1000) for level in levels}
    axes[0].set_xlim(0, 20)
    axes[0].set_ylim(0, 2250)
    axes[0].clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=8)

    for c in cm1.collections:
        c.set_edgecolor("face")
    axes[0].set_title(r'$\dot\varepsilon_{xz}$')
    axes[1].set_title('Vertical fabric')
    axes[2].set_title('Angle')
    axes[0].set_xlim(0, 20)
    axes[0].set_xticks([0., 10., 20.])

    cm2 = axes[1].tricontourf(
        tri,
        vtu.rawdata_dict['fabric 2'],
        cmap='BrBG',
        norm=anisolib.fab_norm,
        levels=np.linspace(
            0,
            1,
            101),
        extend='neither')
    for c in cm2.collections:
        c.set_edgecolor("face")
    axes[1].set_xlim(0, 20)
    axes[1].set_xticks([0., 10., 20.])

    cm3 = axes[2].tricontourf(tri, vtu.rawdata_dict['eigenv 4'], cmap='PuOr',
                              levels=np.linspace(-5, 5, 201), extend='both')
    for c in cm3.collections:
        c.set_edgecolor("face")
    axes[2].set_xlim(0, 20)
    axes[2].set_xticks([0., 10., 20.])

    cbr = plt.colorbar(cm1, cax=cbr_axes[0], orientation='horizontal',
                       label=r'$\dot\varepsilon_{xz}$', ticks=(-2.5e-4, 0, 2.5e-4))
    cbr.ax.set_xticklabels(['-2.5e-4', '0', '2.5e-4'])
    cbr = plt.colorbar(
        cm2,
        cax=cbr_axes[1],
        orientation='horizontal',
        ticks=(
            0,
            1. / 3.,
            2. / 3.,
            1.0),
        label=r'$a^{(2)}_{zz}$')
    cbr.ax.set_xticklabels(['0', r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])
    cbr = plt.colorbar(cm3, cax=cbr_axes[2], orientation='horizontal', ticks=(-5, 0, 5), label=r'Angle ($^o$)')

    axes[0].set_ylabel('Height (m)')
    for ax in axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    axes[1].set_xlabel('Distance (km)')

    for letter, ax in zip('abc', axes):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)

        ax.set_ylim(0, 2250)
    fig.savefig('../plots/idealized_init_zoom.png', dpi=300)


def idealized_init():
    # We want to find the filename for the highest index where we have all 3
    files = glob('../stream/rstf/deform2e-3_????.vtu')
    inds = np.argsort(np.array([float(fn[-8:-4]) for fn in files]))
    files = [files[i] for i in inds]

    out_file = '../stream/relaxed_stream_tall_fabric/deform_diffuse.result'
    vtu = fastvtulib.get_structured_vtu(files[-1])

    vtu.rawdata_dict['vel'] = np.sqrt(vtu.rawdata_dict['aiflow'][:, 0]**2. + vtu.rawdata_dict['aiflow'][:, 1]**2.)

    tri = Triangulation(np.array([rc[0] + 100000. for rc in vtu.raw_coords[:, 0]]) /
                        1000., np.array([rc[1] for rc in vtu.raw_coords[:, 0]]), vtu.simptt)
    gs = gridspec.GridSpec(
        2,
        3,
        width_ratios=(
            1.,
            1.,
            1.),
        height_ratios=(
            1,
            0.1),
        hspace=0.55,
        wspace=0.1,
        left=0.1,
        bottom=0.2,
        top=0.9,
        right=0.98)

    fig = plt.figure(figsize=(7.5, 3.0))
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    cbr_axes = [fig.add_subplot(gs[1, i]) for i in range(3)]

    cm1 = axes[0].tricontourf(tri, vtu.rawdata_dict['aiflow'][:, 0], cmap='Reds',
                              levels=np.linspace(0.0, 10., 101), extend='max')
    levels = [1000, 5000, 10000, 50000, 100000]
    CS = axes[0].tricontour(tri, vtu.rawdata_dict['age'], colors='k', levels=levels, linewidths=0.5)
    fmt = {level: '{:d} kyr'.format(level // 1000) for level in levels}
    axes[0].clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=8)

    for c in cm1.collections:
        c.set_edgecolor("face")
    # axes[0].set_title('Speed', fontsize=bfs)
    # axes[1].set_title('Vertical fabric', fontsize=bfs)
    # axes[2].set_title('Woodcock Parameter', fontsize=bfs)
    axes[0].set_xlim(0, 175)
    axes[0].set_xticks([0., 50., 100., 150.])

    cm2 = axes[1].tricontourf(
        tri,
        vtu.rawdata_dict['fabric 2'],
        cmap='BrBG',
        norm=anisolib.fab_norm,
        levels=np.linspace(
            0,
            1,
            101),
        extend='neither')
    for c in cm2.collections:
        c.set_edgecolor("face")
    axes[1].set_xlim(0, 175)
    axes[1].set_xticks([0., 50., 100., 150.])

    k = woodcock(vtu.rawdata_dict)
    cm3 = axes[2].tricontourf(tri, k, cmap='PuOr', levels=np.linspace(0, 2, 201), extend='max')
    for c in cm3.collections:
        c.set_edgecolor("face")
    axes[2].set_xlim(0, 175)
    axes[2].set_xticks([0., 50., 100., 150.])

    cbr = plt.colorbar(cm1, cax=cbr_axes[0], orientation='horizontal', ticks=(0, 5., 10.))
    cbr.set_label(label=r'$v_x$ (m a$^{-1})$', size=fs)
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)
    cbr = plt.colorbar(cm2, cax=cbr_axes[1], orientation='horizontal', ticks=(0, 1. / 3., 2. / 3., 1.0))
    cbr.set_label(label=r'$a^{(2)}_{zz}$', size=fs)
    cbr.ax.set_xticklabels(['0', r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)
    cbr = plt.colorbar(cm3, cax=cbr_axes[2], orientation='horizontal', ticks=(0, 1, 2))
    cbr.set_label(label=r'Woodcock $k$', size=fs)
    cbr.ax.set_xticklabels(['0', '1', '2'])
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)

    axes[0].set_ylabel('Height (m)', fontsize=fs)
    for ax in axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    axes[1].set_xlabel('Distance (km)', fontsize=fs)

    for letter, ax in zip('abc', axes):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis='both', which='major', labelsize=fs)

        ax.set_ylim(0, 2250)
    fig.savefig('../plots/idealized_init.png', dpi=300)
    fig.savefig('../plots/poster_idealized_init.png', dpi=300, transparent=True)


def idealized_with_cartoon():
    # We want to find the filename for the highest index where we have all 3
    files = glob('../stream/rstf/deform2e-3_????.vtu')
    inds = np.argsort(np.array([float(fn[-8:-4]) for fn in files]))
    files = [files[i] for i in inds]

    vtu = fastvtulib.get_structured_vtu(files[-1])

    vtu.rawdata_dict['vel'] = np.sqrt(vtu.rawdata_dict['aiflow'][:, 0]**2. + vtu.rawdata_dict['aiflow'][:, 1]**2.)

    tri = Triangulation(np.array([rc[0] + 100000. for rc in vtu.raw_coords[:, 0]]) /
                        1000., np.array([rc[1] for rc in vtu.raw_coords[:, 0]]), vtu.simptt)
    gs = gridspec.GridSpec(
        7,
        4,
        width_ratios=(
            1.,
            1.,
            0.5,
            0.035),
        height_ratios=(
            0.4,
            0.05,
            0.4,
            0.05,
            0.2,
            0.05,
            0.15),
        hspace=0.00,
        wspace=0.1,
        left=0.07,
        bottom=0.07,
        top=0.93,
        right=0.92)

    fig = plt.figure(figsize=(7.9, 3.1))
    axes = [fig.add_subplot(gs[0:4, i]) for i in range(2)]
    cbr_axes = [fig.add_subplot(gs[5, i]) for i in range(2)]
    
    ball_cbr_axis = fig.add_subplot(gs[0:5, 3])

    cartoon_axes = [fig.add_subplot(gs[0, 2], projection=fabricplotlib.PRJ),
                    fig.add_subplot(gs[2, 2], projection=fabricplotlib.PRJ),
                    fig.add_subplot(gs[4:, 2], projection=fabricplotlib.PRJ)]

    cm1 = axes[0].tricontourf(tri, vtu.rawdata_dict['aiflow'][:, 0], cmap='Reds',
                              levels=np.linspace(0.0, 10., 101), extend='max')
    levels = [1000, 5000, 10000, 50000, 100000]
    CS = axes[0].tricontour(tri, vtu.rawdata_dict['age'], colors='k', levels=levels, linewidths=0.5)
    fmt = {level: '{:d} kyr'.format(level // 1000) for level in levels}
    axes[0].clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=8)

    for c in cm1.collections:
        c.set_edgecolor("face")
    # axes[0].set_title('Speed', fontsize=bfs)
    # axes[1].set_title('Vertical fabric', fontsize=bfs)
    # axes[2].set_title('Woodcock Parameter', fontsize=bfs)
    axes[0].set_xlim(0, 175)
    axes[0].set_xticks([0., 50., 100., 150.])

    # k = woodcock(vtu.rawdata_dict)
    vtu.rawdata_dict['eigenv 3'][np.isnan(vtu.rawdata_dict['eigenv 3'])] = 1.0
    vtu.rawdata_dict['eigenv 3'][np.isinf(vtu.rawdata_dict['eigenv 3'])] = 1.0
    cm3 = axes[1].tricontourf(tri, vtu.rawdata_dict['eigenv 3'], cmap='summer', extend='neither', levels=np.linspace(0.33333, 1, 101))
    for c in cm3.collections:
        c.set_edgecolor("face")
    axes[1].set_xlim(0, 175)
    axes[1].set_xticks([0., 50., 100., 150.])
    fabricplotlib.quiver(axes[1], vtu, scale=15, width=0.005)
    axes[1].scatter(-1000, -1000, marker=r'$\uparrow$', label='Single max. in x-z', color='k')
    axes[1].legend(loc='lower left', bbox_to_anchor=(0.15, 0.85), fontsize=fs, framealpha=1.0)

    cbr = plt.colorbar(cm1, cax=cbr_axes[0], orientation='horizontal', ticks=(0, 5., 10.))
    cbr.set_label(label=r'$v_x$ (m a$^{-1})$', size=fs)
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)

    cbr = plt.colorbar(cm3, cax=cbr_axes[1], orientation='horizontal', ticks=(1. / 3., 2. / 3., 1))
    cbr.set_label(label=r'$a_1^{(2)}$', size=fs)
    cbr.ax.set_xticklabels([r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)

    axes[0].set_ylabel('Height (m)', fontsize=fs)
    for ax in axes[1:]:
        ax.set_yticklabels(['' for tick in ax.get_yticklabels()])
    axes[0].set_xlabel('Distance (km)', fontsize=fs)
    axes[1].set_xlabel('Distance (km)', fontsize=fs)

    for letter, ax in zip('abc', axes):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis='both', which='major', labelsize=fs)

        ax.set_ylim(0, 2250)

    for letter, ax in zip('cde', cartoon_axes):
        ax.text(0.00, 0.9, letter, transform=ax.transAxes, fontsize=bfs)


    x, y = np.array([70000, 70000, 70000]), np.array([1900, 1200, 200])
    axes[1].text(x[0] / 1000. + 100, y[0], 'c', fontsize=fs, ha='center', va='center', bbox=dict(boxstyle='square,pad=0.1', facecolor='white', alpha=0.75))
    axes[1].text(x[1] / 1000. + 100, y[1], 'd', fontsize=fs, ha='center', va='center', bbox=dict(boxstyle='square,pad=0.1', facecolor='white', alpha=0.75))
    axes[1].text(x[2] / 1000. + 100, y[2], 'e', fontsize=fs, ha='center', va='center', bbox=dict(boxstyle='square,pad=0.1', facecolor='white', alpha=0.75))
    fab_at_pts = vtu.get_pts_2d(anisolib.fabs, x, y)
    a2 = anisolib.fabric_dict_to_a2(fab_at_pts)

    for i in range(3):
        cm = fabricplotlib.a2plot(a2[:, :, i], ax=cartoon_axes[i], cbr=False, show=False, levels=13)
    cbr = plt.colorbar(cm, cax=ball_cbr_axis, orientation='vertical')
    cbr.set_label(label=r'ODF($\theta,\phi$)', size=fs)
    # cbr.ax.set_xticklabels(['0', '1', '2'])
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)

    fig.savefig('../plots/idealized_init_cartoon.png', dpi=300)
    fig.savefig('../plots/poster_idealized_init_cartoon.png', dpi=300, transparent=True)


if __name__ == '__main__':
    main()
