#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run the Fujita model on output from Elmer/Ice.

Created on Mon Jul 13 15:48:56 2020

@author: David Lilien (dlilien90@gmail.com)
"""
from matplotlib import rc
rc('font',**{'family':'sans-serif', 'sans-serif':['Linux Biolinum O']})

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import fujita_model
from modeltools.lib import fastvtulib


π = np.pi
fs = 9
bfs = 12


def main():
    """Run the model once or twice."""
    td_plots()
    plot_different_streaming()
    typical_plots()


def td_plots():
    fns = ['../stream/stream_ftw/s_ftw_rc2.0e-3_oop1.0e-3_bm0.0_0001.vtu',
           '../stream/stream_ftw/s_ftw_rc2.0e-3_oop1.0e-3_bm0.0_0102.vtu',
           '../stream/stream_ftw/s_ftw_rc2.0e-3_oop1.0e-3_bm0.0_0202.vtu',
           '../stream/stream_ftw/s_ftw_rc2.0e-3_oop1.0e-3_bm0.0_0302.vtu']

    obs_angles = np.linspace(0, π, 181)
    gs = GridSpec(2, 5, left=0.1, bottom=0.1, right=0.9, top=0.95,
                  wspace=0.14, hspace=0.08, width_ratios=(1, 1, 1, 1, 0.1))
    fig = plt.figure(figsize=(8, 5))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(gs[0, 3], sharex=ax1, sharey=ax1)
    axc = fig.add_subplot(gs[:, 4])
    ax5 = fig.add_subplot(gs[1, 0], sharex=ax1, sharey=ax1)
    ax6 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax1)
    ax7 = fig.add_subplot(gs[1, 2], sharex=ax1, sharey=ax1)
    ax8 = fig.add_subplot(gs[1, 3], sharex=ax1, sharey=ax1)

    y, power, cm = plot_vtu(fns[0],
        ax1, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[1],
        ax2, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[2],
        ax3, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[3],
        ax4, 75000.0, 1000, obs_angles)

    y, power, cm = plot_vtu(fns[0],
        ax5, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[1],
        ax6, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[2],
        ax7, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[3],
        ax8, 75000.0, 1000, obs_angles, scatter_aniso=1.0)

    ax1.set_title('Streaming', fontsize=bfs)
    ax2.set_title('1000 yrs', fontsize=bfs)
    ax3.set_title('2000 yrs', fontsize=bfs)
    ax4.set_title('3000 yrs', fontsize=bfs)
    ax1.set_ylabel('Birefr. and scat.', fontsize=bfs)
    ax5.set_ylabel('Birefr. only', fontsize=bfs)
    plt.colorbar(cm, cax=axc, extend='both')
    axc.set_ylabel('Relative power (dB)', fontsize=fs)
    axc.tick_params(axis='both', which='major', labelsize=fs)
    fig.text(0.02, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize=fs)
    fig.text(0.45, 0.02, r'Observation Angle ($^o$)', va='center', fontsize=fs)

    for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8):
        ax.set_xticks(np.linspace(0, 1, 5) * π)
        ax.set_xticklabels(['0', '45', '90', '135', '180'])
    for ax in (ax1, ax2, ax3, ax4):
        ax.xaxis.set_tick_params(labelbottom=False)

    for ax in (ax2, ax3, ax4, ax6, ax7, ax8):
        ax.yaxis.set_tick_params(labelleft=False)

    for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8)):
        ax.set_ylim(2000, 0)
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis='both', which='major', labelsize=fs)
    plt.savefig('../plots/pres_oop_td_fwd.pdf', dpi=300)


    fns = ['../stream/rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_0001.vtu',
           '../stream/rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_0102.vtu',
           '../stream/rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_0202.vtu',
           '../stream/rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_0302.vtu']

    obs_angles = np.linspace(0, π, 181)
    gs = GridSpec(2, 5, left=0.1, bottom=0.1, right=0.9, top=0.95,
                  wspace=0.14, hspace=0.08, width_ratios=(1, 1, 1, 1, 0.1))
    fig = plt.figure(figsize=(8, 5))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(gs[0, 3], sharex=ax1, sharey=ax1)
    axc = fig.add_subplot(gs[:, 4])
    ax5 = fig.add_subplot(gs[1, 0], sharex=ax1, sharey=ax1)
    ax6 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax1)
    ax7 = fig.add_subplot(gs[1, 2], sharex=ax1, sharey=ax1)
    ax8 = fig.add_subplot(gs[1, 3], sharex=ax1, sharey=ax1)


    y, power, cm = plot_vtu(fns[0],
        ax1, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[1],
        ax2, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[2],
        ax3, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[3],
        ax4, 75000.0, 1000, obs_angles)

    y, power, cm = plot_vtu(fns[0],
        ax5, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[1],
        ax6, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[2],
        ax7, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[3],
        ax8, 75000.0, 1000, obs_angles, scatter_aniso=1.0)

    ax1.set_title('Deforming', fontsize=bfs)
    ax2.set_title('1000 yrs', fontsize=bfs)
    ax3.set_title('2000 yrs', fontsize=bfs)
    ax4.set_title('3000 yrs', fontsize=bfs)
    ax1.set_ylabel('Birefr. and scat.', fontsize=bfs)
    ax5.set_ylabel('Birefr. only', fontsize=bfs)
    plt.colorbar(cm, cax=axc, extend='both')
    axc.set_ylabel('Relative power (dB)', fontsize=fs)
    axc.tick_params(axis='both', which='major', labelsize=fs)
    fig.text(0.02, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize=fs)
    fig.text(0.45, 0.02, r'Observation Angle ($^o$)', va='center', fontsize=fs)

    for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8):
        ax.set_xticks(np.linspace(0, 1, 5) * π)
        ax.set_xticklabels(['0', '45', '90', '135', '180'])
    for ax in (ax1, ax2, ax3, ax4):
        ax.xaxis.set_tick_params(labelbottom=False)

    for ax in (ax2, ax3, ax4, ax6, ax7, ax8):
        ax.yaxis.set_tick_params(labelleft=False)

    for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8)):
        ax.set_ylim(2000, 0)
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis='both', which='major', labelsize=fs)
    plt.savefig('../plots/pres_ftw_td_fwd.pdf', dpi=300)

    fns = ['../stream/rstf/d_ftw_rc2.0e-3_oop0.0_bm0.0_0001.vtu',
           '../stream/rstf/d_ftw_rc2.0e-3_oop0.0_bm0.0_0502.vtu',
           '../stream/rstf/d_ftw_rc2.0e-3_oop0.0_bm0.0_1002.vtu',
           '../stream/rstf/d_ftw_rc2.0e-3_oop0.0_bm0.0_1502.vtu']

    obs_angles = np.linspace(0, π, 181)
    gs = GridSpec(2, 5, left=0.1, bottom=0.1, right=0.9, top=0.95,
                  wspace=0.14, hspace=0.08, width_ratios=(1, 1, 1, 1, 0.1))
    fig = plt.figure(figsize=(8, 5))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(gs[0, 3], sharex=ax1, sharey=ax1)
    axc = fig.add_subplot(gs[:, 4])
    ax5 = fig.add_subplot(gs[1, 0], sharex=ax1, sharey=ax1)
    ax6 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax1)
    ax7 = fig.add_subplot(gs[1, 2], sharex=ax1, sharey=ax1)
    ax8 = fig.add_subplot(gs[1, 3], sharex=ax1, sharey=ax1)


    y, power, cm = plot_vtu(fns[0],
        ax1, 50000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[1],
        ax2, 50000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[2],
        ax3, 50000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[3],
        ax4, 50000.0, 1000, obs_angles)

    y, power, cm = plot_vtu(fns[0],
        ax5, 50000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[1],
        ax6, 50000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[2],
        ax7, 50000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[3],
        ax8, 50000.0, 1000, obs_angles, scatter_aniso=1.0)

    ax1.set_title('Streaming', fontsize=bfs)
    ax2.set_title('5000 yrs', fontsize=bfs)
    ax3.set_title('10000 yrs', fontsize=bfs)
    ax4.set_title('15000 yrs', fontsize=bfs)
    ax1.set_ylabel('Birefr. and scat.', fontsize=bfs)
    ax5.set_ylabel('Birefr. only', fontsize=bfs)
    plt.colorbar(cm, cax=axc, extend='both')
    axc.set_ylabel('Relative power (dB)', fontsize=fs)
    axc.tick_params(axis='both', which='major', labelsize=fs)
    fig.text(0.02, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize=fs)
    fig.text(0.45, 0.02, r'Observation Angle ($^o$)', va='center', fontsize=fs)

    for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8):
        ax.set_xticks(np.linspace(0, 1, 5) * π)
        ax.set_xticklabels(['0', '45', '90', '135', '180'])
    for ax in (ax1, ax2, ax3, ax4):
        ax.xaxis.set_tick_params(labelbottom=False)

    for ax in (ax2, ax3, ax4, ax6, ax7, ax8):
        ax.yaxis.set_tick_params(labelleft=False)

    for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8)):
        ax.set_ylim(2000, 0)
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis='both', which='major', labelsize=fs)
    plt.savefig('../plots/pres_ftw_td_bck_50.pdf', dpi=300)


    obs_angles = np.linspace(0, π, 181)
    gs = GridSpec(2, 5, left=0.1, bottom=0.1, right=0.9, top=0.95,
                  wspace=0.14, hspace=0.08, width_ratios=(1, 1, 1, 1, 0.1))
    fig = plt.figure(figsize=(8, 5))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(gs[0, 3], sharex=ax1, sharey=ax1)
    axc = fig.add_subplot(gs[:, 4])
    ax5 = fig.add_subplot(gs[1, 0], sharex=ax1, sharey=ax1)
    ax6 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax1)
    ax7 = fig.add_subplot(gs[1, 2], sharex=ax1, sharey=ax1)
    ax8 = fig.add_subplot(gs[1, 3], sharex=ax1, sharey=ax1)


    y, power, cm = plot_vtu(fns[0],
        ax1, 25000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[1],
        ax2, 25000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[2],
        ax3, 25000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[3],
        ax4, 25000.0, 1000, obs_angles)

    y, power, cm = plot_vtu(fns[0],
        ax5, 25000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[1],
        ax6, 25000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[2],
        ax7, 25000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[3],
        ax8, 25000.0, 1000, obs_angles, scatter_aniso=1.0)

    ax1.set_title('Streaming', fontsize=bfs)
    ax2.set_title('5000 yrs', fontsize=bfs)
    ax3.set_title('10000 yrs', fontsize=bfs)
    ax4.set_title('15000 yrs', fontsize=bfs)
    ax1.set_ylabel('Birefr. and scat.', fontsize=bfs)
    ax5.set_ylabel('Birefr. only', fontsize=bfs)
    plt.colorbar(cm, cax=axc, extend='both')
    axc.set_ylabel('Relative power (dB)', fontsize=fs)
    axc.tick_params(axis='both', which='major', labelsize=fs)
    fig.text(0.02, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize=fs)
    fig.text(0.45, 0.02, r'Observation Angle ($^o$)', va='center', fontsize=fs)

    for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8):
        ax.set_xticks(np.linspace(0, 1, 5) * π)
        ax.set_xticklabels(['0', '45', '90', '135', '180'])
    for ax in (ax1, ax2, ax3, ax4):
        ax.xaxis.set_tick_params(labelbottom=False)

    for ax in (ax2, ax3, ax4, ax6, ax7, ax8):
        ax.yaxis.set_tick_params(labelleft=False)

    for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8)):
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis='both', which='major', labelsize=fs)
        ax.set_ylim(2000, 0)
    plt.savefig('../plots/pres_ftw_td_bck_25.pdf', dpi=300)

    obs_angles = np.linspace(0, π, 181)
    gs = GridSpec(2, 5, left=0.1, bottom=0.1, right=0.9, top=0.95,
                  wspace=0.14, hspace=0.08, width_ratios=(1, 1, 1, 1, 0.1))
    fig = plt.figure(figsize=(8, 5))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(gs[0, 3], sharex=ax1, sharey=ax1)
    axc = fig.add_subplot(gs[:, 4])
    ax5 = fig.add_subplot(gs[1, 0], sharex=ax1, sharey=ax1)
    ax6 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax1)
    ax7 = fig.add_subplot(gs[1, 2], sharex=ax1, sharey=ax1)
    ax8 = fig.add_subplot(gs[1, 3], sharex=ax1, sharey=ax1)


    y, power, cm = plot_vtu(fns[0],
        ax1, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[1],
        ax2, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[2],
        ax3, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[3],
        ax4, 75000.0, 1000, obs_angles)

    y, power, cm = plot_vtu(fns[0],
        ax5, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[1],
        ax6, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[2],
        ax7, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[3],
        ax8, 75000.0, 1000, obs_angles, scatter_aniso=1.0)

    ax1.set_title('Streaming', fontsize=bfs)
    ax2.set_title('5000 yrs', fontsize=bfs)
    ax3.set_title('10000 yrs', fontsize=bfs)
    ax4.set_title('15000 yrs', fontsize=bfs)
    ax1.set_ylabel('Birefr. and scat.', fontsize=bfs)
    ax5.set_ylabel('Birefr. only', fontsize=bfs)
    plt.colorbar(cm, cax=axc, extend='both')
    axc.set_ylabel('Relative power (dB)', fontsize=fs)
    axc.tick_params(axis='both', which='major', labelsize=fs)
    fig.text(0.02, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize=fs)
    fig.text(0.45, 0.02, r'Observation Angle ($^o$)', va='center', fontsize=fs)

    for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8):
        ax.set_xticks(np.linspace(0, 1, 5) * π)
        ax.set_xticklabels(['0', '45', '90', '135', '180'])
    for ax in (ax1, ax2, ax3, ax4):
        ax.xaxis.set_tick_params(labelbottom=False)

    for ax in (ax2, ax3, ax4, ax6, ax7, ax8):
        ax.yaxis.set_tick_params(labelleft=False)

    for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8)):
        ax.set_ylim(2000, 0)
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis='both', which='major', labelsize=fs)
    plt.savefig('../plots/pres_ftw_td_bck.pdf', dpi=300)

    fns = ['../stream/rstf/d_ftw_rc2.0e-3_oop0.0_bm0.0_0001.vtu',
           '../stream/rstf/d_ftw_rc2.0e-3_oop0.0_bm0.0_0502.vtu',
           '../stream/rstf/d_ftw_rc2.0e-3_oop0.0_bm0.0_1002.vtu',
           '../stream/rstf/d_ftw_rc2.0e-3_oop0.0_bm0.0_1502.vtu',
           '../stream/rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_0001.vtu']
    gs = GridSpec(1, 6, left=0.1, bottom=0.15, right=0.9, top=0.93,
                  wspace=0.14, hspace=0.08, width_ratios=(1, 1, 1, 1, 1, 0.1))
    fig = plt.figure(figsize=(8, 3))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(gs[0, 3], sharex=ax1, sharey=ax1)
    ax5 = fig.add_subplot(gs[0, 4], sharex=ax1, sharey=ax1)
    axc = fig.add_subplot(gs[:, 5])


    y, power, cm = plot_vtu(fns[0],
        ax1, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[1],
        ax2, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[2],
        ax3, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[3],
        ax4, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[4],
        ax5, 75000.0, 1000, obs_angles)

    ax1.set_title('Streaming', fontsize=bfs)
    ax2.set_title('5000 yrs', fontsize=bfs)
    ax3.set_title('10000 yrs', fontsize=bfs)
    ax4.set_title('15000 yrs', fontsize=bfs)
    ax5.set_title('Pure deform.', fontsize=bfs)
    ax1.set_ylabel('Depth (m)', fontsize=bfs)
    plt.colorbar(cm, cax=axc, extend='both')
    axc.set_ylabel('Relative power (dB)', fontsize=fs)
    axc.tick_params(axis='both', which='major', labelsize=fs)
    fig.text(0.45, 0.02, r'Observation Angle ($^o$)', va='center', fontsize=fs)

    for ax in (ax1, ax2, ax3, ax4, ax5):
        ax.set_xticks(np.linspace(0, 1, 5) * π)
        ax.set_xticklabels(['0', '45', '90', '135', '180'])

    for ax in (ax2, ax3, ax4, ax5):
        ax.yaxis.set_tick_params(labelleft=False)

    for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4, ax5)):
        ax.set_ylim(2000, 0)
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis='both', which='major', labelsize=fs)
    plt.savefig('../plots/poster_pres_ftw_td_bck.png', transparent=True, dpi=300)


def typical_plots():

    fns = ['../stream/rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_0001.vtu',
           '../stream/rstf/s_uc_rc2.0e-3_oop0.0_bm0.0_0502.vtu',
           '../stream/rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_0502.vtu',
           '../stream/stream_ftw/s_ftw_rc2.0e-3_oop1.0e-2_bm0.0_0302.vtu']

    obs_angles = np.linspace(0, π, 181)
    gs = GridSpec(2, 5, left=0.1, bottom=0.1, right=0.9, top=0.95,
                  wspace=0.14, hspace=0.08, width_ratios=(1, 1, 1, 1, 0.1))
    fig = plt.figure(figsize=(8, 5))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(gs[0, 3], sharex=ax1, sharey=ax1)
    axc = fig.add_subplot(gs[:, 4])
    ax5 = fig.add_subplot(gs[1, 0], sharex=ax1, sharey=ax1)
    ax6 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax1)
    ax7 = fig.add_subplot(gs[1, 2], sharex=ax1, sharey=ax1)
    ax8 = fig.add_subplot(gs[1, 3], sharex=ax1, sharey=ax1)


    y, power, cm = plot_vtu(fns[0],
        ax1, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[1],
        ax2, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[2],
        ax3, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[3],
        ax4, 75000.0, 1000, obs_angles)

    y, power, cm = plot_vtu(fns[0],
        ax5, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[1],
        ax6, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[2],
        ax7, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[3],
        ax8, 75000.0, 1000, obs_angles, scatter_aniso=1.0)

    ax1.set_title('Deforming', fontsize=bfs)
    ax2.set_title('Conv. upstream', fontsize=bfs)
    ax3.set_title('Conv. everywhere', fontsize=bfs)
    ax4.set_title(r'$\dot\varepsilon_{xy}^{max}=0.01$ a$^{-1}$', fontsize=bfs)
    ax1.set_ylabel('Birefr. and scat.', fontsize=bfs)
    ax5.set_ylabel('Birefr. only', fontsize=bfs)
    plt.colorbar(cm, cax=axc, extend='both')
    axc.set_ylabel('Relative power (dB)', fontsize=fs)
    axc.tick_params(axis='both', which='major', labelsize=fs)
    fig.text(0.02, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize=fs)
    fig.text(0.45, 0.02, r'Observation Angle ($^o$)', va='center', fontsize=fs)

    for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8):
        ax.set_xticks(np.linspace(0, 1, 5) * π)
        ax.set_xticklabels(['0', '45', '90', '135', '180'])
    for ax in (ax1, ax2, ax3, ax4):
        ax.xaxis.set_tick_params(labelbottom=False)

    for ax in (ax2, ax3, ax4, ax6, ax7, ax8):
        ax.yaxis.set_tick_params(labelleft=False)

    for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8)):
        ax.set_ylim(2000, 0)
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis='both', which='major', labelsize=fs)
    plt.savefig('../plots/pres_modeling.pdf', dpi=300)

    fns = ['../stream/rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_0502.vtu',
           '../stream/stream_ftw/s_ftw_rc2.0e-3_oop1.0e-3_bm0.0_0302.vtu',
           '../stream/stream_ftw/s_ftw_rc2.0e-3_oop1.0e-2_bm0.0_0302.vtu',
           '../stream/stream_ftw/s_ftw_rc2.0e-3_oop1.0e-1_bm0.0_0302.vtu']

    obs_angles = np.linspace(0, π, 181)
    gs = GridSpec(2, 5, left=0.1, bottom=0.1, right=0.9, top=0.95,
                  wspace=0.14, hspace=0.08, width_ratios=(1, 1, 1, 1, 0.1))
    fig = plt.figure(figsize=(8, 5))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(gs[0, 3], sharex=ax1, sharey=ax1)
    axc = fig.add_subplot(gs[:, 4])
    ax5 = fig.add_subplot(gs[1, 0], sharex=ax1, sharey=ax1)
    ax6 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax1)
    ax7 = fig.add_subplot(gs[1, 2], sharex=ax1, sharey=ax1)
    ax8 = fig.add_subplot(gs[1, 3], sharex=ax1, sharey=ax1)


    y, power, cm = plot_vtu(fns[0],
        ax1, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[1],
        ax2, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[2],
        ax3, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[3],
        ax4, 75000.0, 1000, obs_angles)

    y, power, cm = plot_vtu(fns[0],
        ax5, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[1],
        ax6, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[2],
        ax7, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[3],
        ax8, 75000.0, 1000, obs_angles, scatter_aniso=1.0)

    ax1.set_title('Deforming', fontsize=bfs)
    ax2.set_title(r'$\dot\varepsilon_{xy}^{max}=0.001$ a$^{-1}$', fontsize=bfs)
    ax3.set_title(r'$\dot\varepsilon_{xy}^{max}=0.01$ a$^{-1}$', fontsize=bfs)
    ax4.set_title(r'$\dot\varepsilon_{xy}^{max}=0.1$ a$^{-1}$', fontsize=bfs)
    ax1.set_ylabel('Birefr. and scat.', fontsize=bfs)
    ax5.set_ylabel('Birefr. only', fontsize=bfs)
    plt.colorbar(cm, cax=axc, extend='both')
    axc.set_ylabel('Relative power (dB)', fontsize=fs)
    axc.tick_params(axis='both', which='major', labelsize=fs)
    fig.text(0.02, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize=fs)
    fig.text(0.45, 0.02, 'Observation Angle (Degrees)', va='center', fontsize=fs)

    for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8):
        ax.set_xticks(np.linspace(0, 1, 5) * π)
        ax.set_xticklabels(['0', '45', '90', '135', '180'])
    for ax in (ax1, ax2, ax3, ax4):
        ax.xaxis.set_tick_params(labelbottom=False)

    for ax in (ax2, ax3, ax4, ax6, ax7, ax8):
        ax.yaxis.set_tick_params(labelleft=False)

    for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8)):
        ax.set_ylim(2000, 0)
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis='both', which='major', labelsize=fs)
    plt.savefig('../plots/pres_oop_steady.pdf', dpi=300)



def plot_different_streaming():
    fns = ['../stream/rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_0001.vtu',
           '../stream/rstf/s_uc_rc2.0e-3_oop0.0_bm0.0_0502.vtu',
           '../stream/rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_0502.vtu',
           '../stream/rstf/s_hiacc_rc2.0e-3_oop0.0_bm0.0_0502.vtu']

    obs_angles = np.linspace(0, π, 181)
    gs = GridSpec(2, 5, left=0.1, bottom=0.1, right=0.9, top=0.95,
                  wspace=0.14, hspace=0.08, width_ratios=(1, 1, 1, 1, 0.1))
    fig = plt.figure(figsize=(8, 5))
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(gs[0, 2], sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(gs[0, 3], sharex=ax1, sharey=ax1)
    axc = fig.add_subplot(gs[:, 4])
    ax5 = fig.add_subplot(gs[1, 0], sharex=ax1, sharey=ax1)
    ax6 = fig.add_subplot(gs[1, 1], sharex=ax1, sharey=ax1)
    ax7 = fig.add_subplot(gs[1, 2], sharex=ax1, sharey=ax1)
    ax8 = fig.add_subplot(gs[1, 3], sharex=ax1, sharey=ax1)


    y, power, cm = plot_vtu(fns[0],
        ax1, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[1],
        ax2, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[2],
        ax3, 75000.0, 1000, obs_angles)
    y, power, cm = plot_vtu(fns[3],
        ax4, 75000.0, 1000, obs_angles)

    y, power, cm = plot_vtu(fns[0],
        ax5, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[1],
        ax6, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[2],
        ax7, 75000.0, 1000, obs_angles, scatter_aniso=1.0)
    y, power, cm = plot_vtu(fns[3],
        ax8, 75000.0, 1000, obs_angles, scatter_aniso=1.0)

    ax1.set_title('Deforming', fontsize=bfs)
    ax3.set_title('Conv. Upstream', fontsize=bfs)
    ax2.set_title('Conv. Everywhere', fontsize=bfs)
    ax4.set_title('High acc.', fontsize=bfs)
    ax1.set_ylabel('Birefr. and scat.', fontsize=bfs)
    ax5.set_ylabel('Birefr. only', fontsize=bfs)
    plt.colorbar(cm, cax=axc, extend='both')
    axc.set_ylabel('Relative power (dB)', fontsize=fs)
    axc.tick_params(axis='both', which='major', labelsize=fs)
    fig.text(0.02, 0.5, 'Depth (m)', va='center', rotation='vertical', fontsize=fs)
    fig.text(0.45, 0.02, r'Observation Angle ($^o$)', va='center', fontsize=fs)

    for ax in (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8):
        ax.set_xticks(np.linspace(0, 1, 5) * π)
        ax.set_xticklabels(['0', '45', '90', '135', '180'])
    for ax in (ax1, ax2, ax3, ax4):
        ax.xaxis.set_tick_params(labelbottom=False)

    for ax in (ax2, ax3, ax4, ax6, ax7, ax8):
        ax.yaxis.set_tick_params(labelleft=False)

    for letter, ax in zip('abcdefghijklmnop', (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8)):
        ax.set_ylim(2000, 0)
        ax.text(0.05, 0.85, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis='both', which='major', labelsize=fs)
    plt.savefig('../plots/pres_modeling_streams.pdf', dpi=300)


def plot_vtu(fn, ax, x, y_pts, obs_angles,
             E_transmit=np.array([1.0e3, 0]), scatter_aniso=1.0e-1,
             freq=300.0e6, vmin=-30, vmax=10):
    """Load elmer results and plot them on the given axis."""
    y, power = vtu_to_φ_and_ε(fn, x, y_pts, obs_angles,
                              E_transmit=E_transmit,
                              scatter_aniso=scatter_aniso, freq=freq)
    cm = ax.imshow(fujita_model.power_anomaly(power[:, :, 0]),
                   extent=[np.min(obs_angles), np.max(obs_angles),
                           y[0] - y[-1], 0],
                   vmin=vmin, vmax=vmax)
    ax.set_aspect(aspect="auto")
    return y, power, cm


def vtu_to_φ_and_ε(fn, x, y_pts, obs_angles,
                   E_transmit=np.array([1.0e3, 0]), scatter_aniso=1.0e-1,
                   freq=300.0e6):
    vtu = fastvtulib.get_structured_vtu(fn)
    y = np.linspace(np.min(vtu.coords[:, 1]), np.max(vtu.coords[:, 1]), y_pts)
    x = np.ones((y_pts,)) * x
    test_p = vtu.get_pts_2d(['fabric 1'], x, y)
    y_valid = y[~np.isnan(test_p['fabric 1'])]
    y = np.linspace(np.min(y_valid), np.max(y_valid), y_pts)[::-1]
    depths = y[0] - y
    pts = vtu.get_pts_2d(['fabric 1', 'fabric 2', 'fabric 3',
                          'fabric 4', 'fabric 5', 'temp'], x, y)
    A_xy = np.zeros((len(y), 2, 2), dtype=np.complex128)
    A_xy[:, 0, 0] = pts['fabric 1']
    A_xy[:, 1, 1] = 1.0 - pts['fabric 1'] - pts['fabric 2']
    A_xy[:, 1, 0] = pts['fabric 5']
    A_xy[:, 0, 1] = pts['fabric 5']
    φ, ε = fujita_model.φ_and_ε(A_xy, T=pts['temp'])

    if fn == '../stream/rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_0001.vtu':
        with open('deforming_fabric.csv', 'w') as fout:
            fout.write('depth[m] epsilon_x epsilon_y\n')
            for i in range(len(depths)):
                fout.write('{:f} {:f} {:f}\n'.format(depths[i], np.real(ε[i, 0]), np.real(ε[i, 1])))


    # Anisotropic scattering
    S_aniso = np.zeros((y.shape[0], 2, 2))
    S_aniso[:, 0, 0] = 1.0e0
    S_aniso[:, 1, 1] = scatter_aniso
    S_aniso = S_aniso * 1.0e-4
    return y, fujita_model.fujita_model_εS(E_transmit, depths, obs_angles,
                                           ε, S_aniso, φ, freq)


if __name__ == '__main__':
    main()
