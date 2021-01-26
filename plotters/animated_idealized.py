#! /usr/bin/env python
"""
===========
MovieWriter
===========

This example uses a MovieWriter directly to grab individual frames and write
them to a file. This avoids any event loop integration, but has the advantage
of working with even the Agg backend. This is not recommended for use in an
interactive setting.

"""
# -*- noplot -*-
import os.path
home = os.path.expanduser('~')

from glob import glob

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
import matplotlib
mpl_ver = matplotlib.__version__.split('.')

from modeltools.lib import fastvtulib
from matplotlib.tri import Triangulation
import matplotlib.gridspec as gridspec

from lib import anisolib, fabricplotlib


fabs = ['fabric 1', 'fabric 2', 'fabric 3', 'fabric 4', 'fabric 5']
fs = 8
bfs = 12


def single_plot(fl_n=1, fps=2, fn='slide', panel2='fabric 4', countby=1, folder='rstf'):
    files = glob('../stream/' + folder + '/' + fn + '????.vtu')
    if fn[0] == 'd':
        maxn = 2502
    else:
        maxn = 1002
    if len(files) > maxn + 1:
        files = files[1:maxn]
    else:
        files = files[1:]
    inds = np.argsort(np.array([float(fn[-8:-4]) for fn in files]))
    fns = [files[ind] for ind in inds]
    FFMpegWriter = manimation.writers['ffmpeg']

    title = '{:s}'.format(fn.title())
    metadata = dict(title=title, artist='Matplotlib', comment='Elmer/Ice Model Results')
    writer = FFMpegWriter(fps=fps, metadata=metadata, extra_args=['-vcodec', 'rawvideo'])

    times = np.arange(0.0, 10.0 * len(files), 10.0)
    lens = [len(fn) for fn in fns if len(fn) > 0]

    if len(lens) == 0:
        return
    steps = min(lens)

    vtus = [fastvtulib.get_structured_vtu(tf) if tf is not None else None for tf in fns[1::countby]]

    tris = [Triangulation(np.array([rc[0] for rc in vtu.raw_coords[:, 0]]) / 1000. + 100, np.array([rc[1] for rc in vtu.raw_coords[:, 0]]), vtu.simptt) if vtu is not None else None for vtu in vtus]

    fig, ((ax1), (ax2)) = plt.subplots(nrows=2, ncols=1, figsize=(8, 6), sharex=True, sharey=True)
    ax1.set_title(r't={:5.1f}'.format(0.0))
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
        top=0.90,
        right=0.92)

    fig = plt.figure(figsize=(7.9, 3.1))
    axes = [None, fig.add_subplot(gs[0:4, 0:2])]
    cbr_axes = [None, fig.add_subplot(gs[5, 0])]
    
    ball_cbr_axis = fig.add_subplot(gs[0:5, 3])

    cartoon_axes = [fig.add_subplot(gs[0, 2], projection=fabricplotlib.PRJ),
                    fig.add_subplot(gs[2, 2], projection=fabricplotlib.PRJ),
                    fig.add_subplot(gs[4:, 2], projection=fabricplotlib.PRJ)]

    # k = woodcock(vtu.rawdata_dict)
    vtus[0].rawdata_dict['eigenv 3'][np.isnan(vtus[0].rawdata_dict['eigenv 3'])] = 1.0
    vtus[0].rawdata_dict['eigenv 3'][np.isinf(vtus[0].rawdata_dict['eigenv 3'])] = 1.0
    cm3 = axes[1].tricontourf(tris[0], vtus[0].rawdata_dict['eigenv 3'], cmap='summer', extend='neither', levels=np.linspace(0.33333, 1, 101))
    for c in cm3.collections:
        c.set_edgecolor("face")
    axes[1].set_xlim(0, 175)
    axes[1].set_xticks([0., 50., 100., 150.])
    quiv = fabricplotlib.quiver(axes[1], vtus[0], scale=25, width=0.002)
    axes[1].scatter(-1000, -1000, marker=r'$\uparrow$', label='Single max. in x-z', color='k')
    axes[1].legend(loc='lower left', bbox_to_anchor=(0.15, 0.85), fontsize=fs, framealpha=1.0)

    cbr = plt.colorbar(cm3, cax=cbr_axes[1], orientation='horizontal', ticks=(1. / 3., 2. / 3., 1))
    cbr.set_label(label=r'$a_1^{(2)}$', size=fs)
    cbr.ax.set_xticklabels([r'$\frac{1}{3}$', r'$\frac{2}{3}$', '1'])
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)

    axes[1].set_xlabel('Distance (km)', fontsize=fs)
    axes[1].set_ylabel('Height (m)', fontsize=fs)

    axes[1].text(0.05, 0.85, 'a', transform=axes[1].transAxes, fontsize=bfs)
    axes[1].tick_params(axis='both', which='major', labelsize=fs)
    axes[1].set_ylim(0, 2250)

    t = axes[1].text(-5, 2400, 'Time={:d}'.format(0))

    for letter, ax in zip('bcde', cartoon_axes):
        ax.text(0.00, 0.9, letter, transform=ax.transAxes, fontsize=bfs)

    x, y = np.array([70000, 70000, 70000]), np.array([1900, 1200, 200])
    axes[1].text(x[0] / 1000. + 100, y[0], 'b', fontsize=fs, ha='center', va='center', bbox=dict(boxstyle='square,pad=0.1', facecolor='white', alpha=0.75))
    axes[1].text(x[1] / 1000. + 100, y[1], 'c', fontsize=fs, ha='center', va='center', bbox=dict(boxstyle='square,pad=0.1', facecolor='white', alpha=0.75))
    axes[1].text(x[2] / 1000. + 100, y[2], 'd', fontsize=fs, ha='center', va='center', bbox=dict(boxstyle='square,pad=0.1', facecolor='white', alpha=0.75))
    fab_at_pts = vtus[0].get_pts_2d(anisolib.fabs, x, y)
    a2 = anisolib.fabric_dict_to_a2(fab_at_pts)

    cms = [None for i in range(3)]
    for i in range(3):
        cms[i] = fabricplotlib.a2plot(a2[:, :, i], ax=cartoon_axes[i], cbr=False, show=False, levels=13)
    cbr = plt.colorbar(cms[0], cax=ball_cbr_axis, orientation='vertical')
    cbr.set_label(label=r'ODF($\theta,\phi$)', size=fs)
    # cbr.ax.set_xticklabels(['0', '1', '2'])
    cbr.ax.tick_params(axis='both', which='major', labelsize=fs)

    if folder == 'rstf':
        pref = ''
    else:
        pref = 'oop'
    with writer.saving(fig, '{:s}_fabric_evolution.avi'.format(pref + fn), 150):
        for i, (svtu, tri, time) in enumerate(zip(vtus, tris, times[0::countby])):
            ax1.set_title(r't={:5.1f}'.format(time))
            if svtu is not None:
                # Note that this is not fully accurate since fabric 2 only approximates the maximum eigenvalue
                for j in range(3):
                    for c in cms[j].collections:
                        c.remove()
                for c in cm3.collections:
                    c.remove()
                for p in quiv:
                    try:
                        p.remove()
                    except:
                        try:
                            for c in p.collections:
                                c.remove()
                        except AttributeError:
                            # I think this is the dots for OOP
                            for c in p:
                                c.remove()

                t.set_text('Time={:d}'.format(int(time)))

                cm3 = axes[1].tricontourf(tri, svtu.rawdata_dict['eigenv 3'], cmap='summer', extend='neither', levels=np.linspace(0.33333, 1, 101))
                quiv = fabricplotlib.quiver(axes[1], svtu, scale=25, width=0.002)

                fab_at_pts = svtu.get_pts_2d(anisolib.fabs, x, y)
                a2 = anisolib.fabric_dict_to_a2(fab_at_pts)
                for j in range(3):
                    cms[j] = fabricplotlib.a2plot(a2[:, :, j], ax=cartoon_axes[j], cbr=False, show=False, levels=13)
            writer.grab_frame()

fwd_countby = 1
fwd_fps = 15
back_countby = 1
back_fps = 30

# Already done forward
# single_plot(panel2='fabric 1', fn='s_midacc_rc2.0e-3_oop0.0_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_uc_rc2.0e-3_oop0.0_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_hiacc_rc2.0e-3_oop0.0_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_ftw_rc2.0e-3_oop0.0_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_ftwlc_rc2.0e-3_oop0.0_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_nuclc_rc2.0e-3_oop0.0_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_ftw_rc2.0e-3_oop1.0e-3_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_ftw_rc2.0e-1_oop0.0_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_ftw_rc2.0e-2_oop0.0_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_ftw_rc2.0e-4_oop0.0_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_ftw_rc2.0e-3_oop1.0e-2_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_nuc_rc2.0e-3_oop0.0_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_ftw_rc2.0_oop0.0_bm0.0_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_ftw_rc2.0e-3_oop0.0_bm1.0e-2_', countby=fwd_countby, fps=fwd_fps)
# single_plot(panel2='fabric 1', fn='s_ftw_rc2.0e-3_oop0.0_bm1.0e-3_', countby=fwd_countby, fps=fwd_fps)

# Already done backward
# single_plot(panel2='fabric 1', fn='d_ftw_rc2.0e-3_oop0.0_bm0.0_', countby=back_countby, fps=back_fps)
# single_plot(panel2='fabric 1', fn='d_uc_rc2.0e-3_oop0.0_bm0.0_', countby=back_countby, fps=back_fps)
# single_plot(panel2='fabric 1', fn='d_hiacc_rc2.0e-3_oop0.0_bm0.0_', countby=back_countby, fps=back_fps)
# single_plot(panel2='fabric 1', fn='d_ftw_rc2.0e-4_oop0.0_bm0.0_', countby=back_countby, fps=back_fps)
# single_plot(panel2='fabric 1', fn='d_ftw_rc2.0e-1_oop0.0_bm0.0_', countby=back_countby, fps=back_fps)
# single_plot(panel2='fabric 1', fn='d_ftw_rc2.0e-3_oop1.0e-2_bm0.0_', countby=back_countby, fps=back_fps, folder='stream_ftw')
# single_plot(panel2='fabric 1', fn='d_ftw_rc2.0e-2_oop0.0_bm0.0_', countby=back_countby, fps=back_fps)
single_plot(panel2='fabric 1', fn='d_ftw_rc2.0e-3_oop0.0_bm2.0e-2_', countby=back_countby, fps=back_fps)

raise
# Still running forward
single_plot(panel2='fabric 1', fn='s_ftw_rc2.0e-3_oop0.0_bm2.0e-2_', countby=fwd_countby, fps=fwd_fps)

# Still running backward
single_plot(panel2='fabric 1', fn='d_ftw_rc2.0e-3_oop1.0e-3_bm0.0_', countby=back_countby, fps=back_fps, folder='stream_ftw')
single_plot(panel2='fabric 1', fn='d_ftw_rc2.0e-3_oop1.0e-1_bm0.0_', countby=back_countby, fps=back_fps, folder='stream_ftw')

single_plot(panel2='fabric 1', fn='d_ftw_rc2.0_oop0.0_bm0.0_', countby=back_countby, fps=back_fps)

single_plot(panel2='fabric 1', fn='d_ftw_rc2.0e-3_oop0.0_bm1.0e-2_', countby=back_countby, fps=back_fps)

