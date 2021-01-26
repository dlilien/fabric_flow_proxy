#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 dlilien <dlilien@hozideh>
#
# Distributed under terms of the MIT license.

"""
Drill a hypothetical core at the 175 km mark
"""

from matplotlib import rc
rc('font',**{'family':'sans-serif', 'sans-serif':['Linux Biolinum O']})
from glob import glob

import matplotlib.pyplot as plt
from modeltools.lib import fastvtulib
import numpy as np
# import matplotlib.gridspec as gridspec
from joblib import Parallel, delayed

fs = 9
bfs = 12

name_fmt = 's_{:s}_rc{:s}_oop{:s}_bm{:s}'
name_dict = {'ftw': 'Converge\nupstream', 'uc': 'Converge\neverywhere', 'hiacc': 'High\naccumulation'}
short_name_dict = {'ftw': 'Conv. upstream', 'uc': 'Conv. everywhere', 'hiacc': 'High acc.'}

ftws = ['ftw', 'uc', 'hiacc']
oops = ['0.0', '1.0e-3', '1.0e-2', '1.0e-1']
bms = ['0.0', '1.0e-2', '2.0e-2']
rcs = ['2.0e-3', '2.0e-4', '2.0e-2']

ftw_names = [name_fmt.format(ftw, rcs[0], oops[0], bms[0]) for ftw in ftws]
bm_names = [name_fmt.format(ftws[0], rcs[0], oops[0], bm) for bm in bms[1:]]
rc_names = [name_fmt.format(ftws[0], rc, oops[0], bms[0]) for rc in rcs[1:]]
oop_names = [name_fmt.format(ftws[0], rcs[0], oop, bms[0]) for oop in oops[1:]]

fmts = ftw_names + bm_names + rc_names + oop_names
folders = ['rstf' for i in ftw_names + bm_names + rc_names] + ['stream_ftw' for i in oop_names]
pretty_fmts = {}
pretty_fmts_oop = {}
for fmt in fmts:
    fm_spl = fmt.split('_')
    pretty_fmts[fmt] = r'%s $\lambda$=%1.4g $\dot b$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[2][2:]), float(fm_spl[4][2:]))
    pretty_fmts_oop[fmt] = r'%s $\dot\varepsilon_{xy}^{(2)}$=%1.4g' % (short_name_dict[fm_spl[1]], float(fm_spl[3][3:]))

files = [glob('../stream/' + folder + '/' + fmt + '_????.vtu') for fmt, folder in zip(fmts, folders)]
inds = [np.argsort(np.array([float(fn[-8:-4]) for fn in filel])) for filel in files]


def get_times(fmt, folder):
    times = []
    out_file = '../stream/' + folder + '/' + fmt + '.result'
    with open(out_file) as fout:
        for line in fout:
            if line[:4] == 'Time':
                times.append(float(line.rstrip('\n').split()[-1]))
    times = np.array(times)
    return np.unique(times - times[0])


n_jobs = len(fmts)
timess = Parallel(n_jobs=n_jobs)(delayed(get_times)(fmt, folder) for fmt, folder in zip(fmts, folders))

for f, t in zip(fmts, timess):
    print(f, t)

maxtime = np.min([np.max(times) for times in timess[:-3]])
maxtime_oop = np.min([np.max(times) for times in timess])

try:
    inds_times = [int(np.where(times == maxtime)[0][0]) for times in timess[:-3]]
    inds_times_oop = [int(np.where(timess[-1] == maxtime_oop)[0][0]) for times in timess]
except:
    inds_times = [301 for i in timess[:-3]]
    inds_times_oop = [301 for i in timess]

vtus = {fmt: fastvtulib.get_structured_vtu(filel[int(ind[ind_t])]) for fmt, ind, ind_t, filel in zip(fmts[:-3], inds, inds_times, files)}
vtus_oop = {fmt: fastvtulib.get_structured_vtu(filel[int(ind[ind_t])]) for fmt, ind, ind_t, filel in zip(fmts, inds, inds_times_oop, files)}

print('With OOP:', maxtime_oop)
print('Without OOP:', maxtime)

deform_files = glob('../stream/rstf/deform_diffuse????.vtu')
deform_inds = np.argsort(np.array([float(fn[-8:-4]) for fn in deform_files]))
deform_vtu = fastvtulib.get_structured_vtu(deform_files[deform_inds[-1]])

colors = ['C0', 'C1', 'C2']
linewidths = [2, 1, 0.5]
linestyles = ['solid', 'dashed', 'dotted']

topy = 2020.
dist = 75.0
x, y = np.ones(200) * 1000.0 * dist, np.linspace(0, topy, 200)


##########
# No OOP
###########
eigs = {fmt: vtus[fmt].get_pts_2d(['eigenv 1', 'eigenv 2', 'eigenv 3'], x, y) for fmt in fmts[:-3]}
a_s = {fmt: vtus[fmt].get_pts_2d(['fabric 1', 'fabric 2', 'fabric 3'], x, y) for fmt in fmts[:-3]}
eigs['deform'] = deform_vtu.get_pts_2d(['eigenv 1', 'eigenv 2', 'eigenv 3'], x, y)
a_s['deform'] = deform_vtu.get_pts_2d(['fabric 1', 'fabric 2', 'fabric 3'], x, y)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 8))
ax1.set_title(r'$S_3$')
ax2.set_title(r'$S_2$')
ax3.set_title(r'$S_1$')
for ax in (ax1, ax2, ax3):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 2000)
ax1.legend(loc='lower left', fontsize=10, frameon=False)
ax2.set_xlabel('Eigenvalue')
ax1.set_ylabel('Elevation (m)')
fig.tight_layout(pad=0.5)

ax1.plot(eigs['deform']['eigenv 3'], y, color='k', linewidth=3, label='Internal deformation')
ax2.plot(eigs['deform']['eigenv 2'], y, color='k', linewidth=3, label='Internal deformation')
ax3.plot(eigs['deform']['eigenv 1'], y, color='k', linewidth=3, label='Internal deformation')
for fmt, color in zip(ftw_names, colors):
    lw = linewidths[0]
    ls = linestyles[0]
    ax1.plot(eigs[fmt]['eigenv 3'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
    ax2.plot(eigs[fmt]['eigenv 2'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
    ax3.plot(eigs[fmt]['eigenv 1'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
ax1.legend(loc='lower left', fontsize=10, frameon=False)
fig.savefig('../plots/idealized_core_eigs_ftw_only.pdf', dpi=300)

for fmt, lw in zip(rc_names, linewidths[1:]):
    color = colors[0]
    ls = linestyles[0]
    ax1.plot(eigs[fmt]['eigenv 3'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
    ax2.plot(eigs[fmt]['eigenv 2'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
    ax3.plot(eigs[fmt]['eigenv 1'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
for fmt, ls in zip(bm_names, linestyles[1:]):
    color = colors[0]
    lw = linewidths[0]
    ax1.plot(eigs[fmt]['eigenv 3'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
    ax2.plot(eigs[fmt]['eigenv 2'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
    ax3.plot(eigs[fmt]['eigenv 1'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
ax1.legend()

fig.savefig('../plots/idealized_core_eigs.pdf', dpi=300)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 8))
ax1.plot(a_s['deform']['fabric 1'], y, color='k', linewidth=3, label='Internal deformation')
ax2.plot(a_s['deform']['fabric 2'], y, color='k', linewidth=3, label='Internal deformation')
ax3.plot(a_s['deform']['fabric 3'], y, color='k', linewidth=3, label='Internal deformation')
for fmt, color in zip(ftw_names, colors):
    lw = linewidths[0]
    ls = linestyles[0]
    ax1.plot(a_s[fmt]['fabric 1'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
    ax2.plot(a_s[fmt]['fabric 2'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
    ax3.plot(a_s[fmt]['fabric 3'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
for fmt, lw in zip(rc_names, linewidths[1:]):
    color = colors[0]
    ls = linestyles[0]
    ax1.plot(a_s[fmt]['fabric 1'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
    ax2.plot(a_s[fmt]['fabric 2'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
    ax3.plot(a_s[fmt]['fabric 3'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
for fmt, ls in zip(bm_names, linestyles[1:]):
    color = colors[0]
    lw = linewidths[0]
    ax1.plot(a_s[fmt]['fabric 1'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
    ax2.plot(a_s[fmt]['fabric 2'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])
    ax3.plot(a_s[fmt]['fabric 3'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts[fmt])

ax1.set_title(r'$a^{(2)}_{xx}$')
ax2.set_title(r'$a^{(2)}_{zz}$')
ax3.set_title(r'$a^{(2)}_{xz}$')
for ax in (ax1, ax2):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 2000)
ax3.set_xlim(-0.15, 0.15)
ax3.set_xticks([-0.1, 0.0, 0.1])
ax3.set_ylim(0, 2000)
ax1.legend(loc='lower right', fontsize=10, frameon=False)
ax2.set_xlabel('Fabric component')
ax1.set_ylabel('Elevation (m)')
fig.tight_layout(pad=0.5)
fig.savefig('../plots/idealized_core_as.pdf', dpi=300)


##########
# OOP
###########
eigs = {fmt: vtus_oop[fmt].get_pts_2d(['eigenv 1', 'eigenv 2', 'eigenv 3'], x, y) for fmt in fmts}
a_s = {fmt: vtus_oop[fmt].get_pts_2d(['fabric 1', 'fabric 2', 'fabric 3'], x, y) for fmt in fmts}
eigs['deform'] = deform_vtu.get_pts_2d(['eigenv 1', 'eigenv 2', 'eigenv 3'], x, y)
a_s['deform'] = deform_vtu.get_pts_2d(['fabric 1', 'fabric 2', 'fabric 3'], x, y)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 8))
ax1.plot(eigs['deform']['eigenv 3'], y, color='k', linewidth=3, label='Internal deformation')
ax2.plot(eigs['deform']['eigenv 2'], y, color='k', linewidth=3, label='Internal deformation')
ax3.plot(eigs['deform']['eigenv 1'], y, color='k', linewidth=3, label='Internal deformation')
for fmt, lw in zip([ftw_names[0]] + oop_names, linewidths):
    color = colors[0]
    ls = linestyles[0]
    ax1.plot(eigs[fmt]['eigenv 3'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts_oop[fmt])
    ax2.plot(eigs[fmt]['eigenv 2'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts_oop[fmt])
    ax3.plot(eigs[fmt]['eigenv 1'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts_oop[fmt])
ax1.set_title(r'$S_3$')
ax2.set_title(r'$S_2$')
ax3.set_title(r'$S_1$')
for ax in (ax1, ax2, ax3):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 2000)
ax1.legend(loc='lower left', fontsize=10, frameon=False)
ax2.set_xlabel('Eigenvalue')
ax1.set_ylabel('Elevation (m)')
fig.tight_layout(pad=0.5)
fig.savefig('../plots/idealized_core_oop_eigs.pdf', dpi=300)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 8))
ax1.plot(a_s['deform']['fabric 1'], y, color='k', linewidth=3, label='Internal deformation')
ax2.plot(a_s['deform']['fabric 2'], y, color='k', linewidth=3, label='Internal deformation')
ax3.plot(a_s['deform']['fabric 3'], y, color='k', linewidth=3, label='Internal deformation')
for fmt, lw in zip([ftw_names[0]] + oop_names, linewidths):
    color = colors[0]
    ls = linestyles[0]
    ax1.plot(a_s[fmt]['fabric 1'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts_oop[fmt])
    ax2.plot(a_s[fmt]['fabric 2'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts_oop[fmt])
    ax3.plot(a_s[fmt]['fabric 3'], y, color=color, linewidth=lw, linestyle=ls, label=pretty_fmts_oop[fmt])
ax1.set_title(r'$a^{(2)}_{xx}$')
ax2.set_title(r'$a^{(2)}_{zz}$')
ax3.set_title(r'$a^{(2)}_{xz}$')
for ax in (ax1, ax2):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 2000)
ax3.set_xlim(-0.15, 0.15)
ax3.set_xticks([-0.1, 0.0, 0.1])
ax3.set_ylim(0, 2000)
ax1.legend(loc='lower right', fontsize=10, frameon=False)
ax2.set_xlabel('Fabric component')
ax1.set_ylabel('Elevation (m)')
fig.tight_layout(pad=0.5)
fig.savefig('../plots/idealized_core_oop_as.pdf', dpi=300)
