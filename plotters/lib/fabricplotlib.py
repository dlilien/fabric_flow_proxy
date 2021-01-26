#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © Nicholas Rathmann, 2020

import numpy as np
import scipy.special as sp
import copy

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
# from .anisolib import cartoon_norm

inclination = 60
rot = -90 * 1
PRJ = ccrs.Orthographic(rot, 90 - inclination)

HOOP = [(0.25 * np.cos(x), np.sin(x)) for x in np.linspace(0, 2 * np.pi, 50)]


def dumb_woodcock(stuff_dict):
    norm = np.abs(stuff_dict['eigenv 3']) + np.abs(stuff_dict['eigenv 2']) + np.abs(stuff_dict['eigenv 1'])
    e1 = np.abs(stuff_dict['eigenv 1']) / norm
    e2 = np.abs(stuff_dict['eigenv 2']) / norm
    e3 = np.abs(stuff_dict['eigenv 3']) / norm
    out = e3 / e2
    out[e1 == 0] = 1e8
    out[e2 == 0] = 1e8
    return out > 5.0


def a2plot(a2, ax=None, cbr=True, show=True,
           latres=20, lonres=40, levels=8):
    phi = np.linspace(0, 2 * np.pi, lonres)  # LON
    theta = np.linspace(0, np.pi, latres)  # CO-LAT
    phi, theta = np.meshgrid(phi, theta)
    lon, colat = phi, theta
    lat = np.pi / 2 - colat

    x = 0
    y = 1
    z = 2
    k = np.sqrt(15.0 / (2.0 * np.pi))

    # ODF expansion coefs from a^2
    c00  =  (a2[x,x]+a2[y,y]+a2[z,z])/np.sqrt(4*np.pi)
    c20  = -1/4*np.sqrt(5/np.pi)*(a2[x,x]+a2[y,y]-2*a2[z,z])
    c22m =  k/4*(a2[x,x]+2j*a2[x,y]-a2[y,y]) # c_2^{-2}
    c22p =  k/4*(a2[x,x]-2j*a2[x,y]-a2[y,y]) # c_2^{+2}
    c21m = +k/2*(a2[x,z]+1j*a2[y,z])
    c21p = -k/2*(a2[x,z]-1j*a2[y,z])

    # Discretize ODF
    lm  = [(0, 0), (2, -2), (2, -1), (2, 0), (2, 1), (2, 2)]
    clm = [c00, c22m, c21m, c20, c21p, c22p]
    ODF = np.sum([clm[ii] * sp.sph_harm(m, l, phi,theta) for ii,(l,m) in enumerate(lm) ], axis=0)
    ODF = np.real(ODF) # Any imag values are rounding errors.

    if ax is None:
        plt.figure()
        gs = gridspec.GridSpec(1, 1)
        ax = plt.subplot(gs[0, 0], projection=PRJ)
        ax.set_global()
    cmap = copy.copy(mpl.cm.get_cmap("Greys"))
    # cmap.set_under('#fee0d2')
    # ODF[ODF < 0.0] = 0.0001
    # h1 = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), ODF, transform=ccrs.PlateCarree(), norm=mpl.colors.LogNorm(0.01, 1), levels=np.logspace(-2, 0, 11), extend='min', cmap="Greys")
    ODF[ODF < 0.0] = 0.0
    h1 = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), ODF, transform=ccrs.PlateCarree(), levels=np.linspace(0, 0.3, levels), extend='max', cmap="Greys")
    ax.gridlines(crs=ccrs.PlateCarree())
    if cbr:
        cb1 = plt.colorbar(h1, ax=ax, fraction=0.035, aspect=19,
                           orientation='horizontal', pad=0.1)
        cb1.set_label(
            r'ODF$(\theta,\phi)$ - mean(ODF) (i.e. isotropic contr. removed)')
    if show:
        print(np.mean(ODF))
        plt.show()
    return h1


def fabric_to_hor_rot(f1, f2, f5):
    """Skip the overly complicated calculations, just go horizontal.

    Essentially assumes that one of the E_i's is vertical.
    Result in degrees.
    """
    A_xy = np.zeros((len(f1), 2, 2))
    A_xy[:, 0, 0] = f1
    A_xy[:, 0, 1] = f5
    A_xy[:, 1, 0] = f5
    A_xy[:, 1, 1] = 1.0 - f1 - f2
    A_xy[np.isnan(A_xy)] = 1.0
    A_xy[np.isinf(A_xy)] = 1.0
    E, rot_mat = np.linalg.eig(A_xy)

    norm = (A_xy[:, 0, 0] + A_xy[:, 1, 1]) / (E[:, 0] + E[:, 1])
    E[:, 0] = E[:, 0] * norm
    E[:, 1] = E[:, 1] * norm
    φ = np.arccos(-rot_mat[:, 0, 0]) * 180.0 / np.pi
    φ[rot_mat[:, 0, 1] < 0.0] = -φ[rot_mat[:, 0, 1] < 0.0]
    φ[φ < -90.0] = φ[φ < -90.0] + 180.0
    # φ[φ < -45.0] = φ[φ < -45.0] + 90.0
    φ[φ > 90.0] = φ[φ > 90.0] - 180.0
    # φ[φ > 45.0] = φ[φ > 45.0] - 90.0

    return φ


def fabric_to_ver_rot(f1, f2, f3):
    """Skip the overly complicated calculations, just go horizontal.

    Essentially assumes that one of the E_i's is vertical.
    Result in degrees.
    """
    A_xy = np.zeros((len(f1), 2, 2))
    A_xy[:, 0, 0] = f2
    A_xy[:, 0, 1] = f3
    A_xy[:, 1, 0] = f3
    A_xy[:, 1, 1] = f1
    A_xy[np.isnan(A_xy)] = 1.0
    A_xy[np.isinf(A_xy)] = 1.0
    E, rot_mat = np.linalg.eig(A_xy)

    φ = np.arccos(-rot_mat[:, 0, 0]) * 180.0 / np.pi
    # φ[rot_mat[:, 0, 1] < 0.0] = -φ[rot_mat[:, 0, 1] < 0.0]
    φ[φ < -90.0] = φ[φ < -90.0] + 180.0
    # φ[φ < -45.0] = φ[φ < -45.0] + 90.0
    φ[φ > 90.0] = φ[φ > 90.0] - 180.0
    # φ[φ > 45.0] = φ[φ > 45.0] - 90.0

    return φ


def quiver(ax, vtu, s=75, x=np.arange(-83333, 100000, 33333.), y=np.arange(100, 2000, 400.), xshift=lambda x: x / 1000. + 100., scale=1.0, width=None):
    X, Y = np.meshgrid(x, y)
    dat = vtu.get_pts_2d(['eigenv 1', 'eigenv 2', 'eigenv 3', 'eigenv 4', 'fabric 1', 'fabric 2', 'fabric 3', 'fabric 5'], X.flatten(), Y.flatten())
    return quiver_dict(ax, dat, s=s, X=xshift(X), Y=Y, scale=scale, width=width)


def quiver_dict(ax, dat, s=75, X=None, Y=None, inx=True, scale=1.0, width=None):
    hang = fabric_to_ver_rot(dat['fabric 1'], dat['fabric 2'], dat['fabric 5'])
    if inx:
        ang = fabric_to_ver_rot(dat['fabric 1'], dat['fabric 2'], dat['fabric 3'])
        ang[ang < -45] = ang[ang < -45] + 90.
        ang[ang > 45] = ang[ang > 45] - 90.
        u = np.sin(ang / 180.0 * np.pi) * dat['eigenv 3']
        v = np.cos(ang / 180.0 * np.pi) * dat['eigenv 3']
        units = 'width'
    else:
        u = np.zeros_like(dat['eigenv 3'])
        v = np.ones_like(dat['eigenv 3'])
        units = 'height'
        ang = np.zeros_like(dat['fabric 1'])

    singlemax = dumb_woodcock(dat)
    vertical = dat['fabric 2'] > (1.0 - dat['fabric 2'] - dat['fabric 1'])
    vertical_sm = np.logical_and(vertical, singlemax)
    hor_sm = np.logical_and(~vertical, singlemax)
    quiv = ax.quiver(X.flatten()[vertical_sm], Y.flatten()[vertical_sm], u[vertical_sm], v[vertical_sm], units=units, scale=scale, width=width)
    if np.any(hor_sm):
        hq = [ax.plot(X.flatten()[hor_sm],
                      Y.flatten()[hor_sm],
                      marker='.',
                      markersize=2,
                      linestyle='none',
                      color='0.4',
                      label='Single max. partly into page')]
    else:
        hq = []

    planlabel = False
    ooplabel = False
    other_pts = []
    for i in range(np.sum(~singlemax)):
        t = MarkerStyle(marker=HOOP)
        t._transform = t.get_transform().rotate_deg(-ang[~singlemax][i])
        if np.isnan(dat['fabric 1'].flatten()[~singlemax][i]):
            continue
        if np.abs(hang.flatten()[~singlemax][i]) > 1:
            color = '0.4'
            if not ooplabel:
                label = 'Vert. girdle, normal out of x-z'
                ooplabel = True
            else:
                label = None
        else:
            color = 'k'
            if not planlabel:
                label = 'Vert. girdle, normal in x-z'
                planlabel = True
            else:
                label = None

        other_pts.append(ax.scatter(X.flatten()[~singlemax][i],
                         Y.flatten()[~singlemax][i],
                         marker=t,
                         s=s,
                         linewidth=0.5,
                         c='none',
                         edgecolors=color,
                         label=label))
    return [quiv] + hq + other_pts


if __name__ == '__main__':
    # Demos
    # a2 = [[0, 0, 0], [0, 0, 0], [0, 0, 1]]  # max in z
    # a2 = [[1, 0, 0], [0, 0, 0], [0, 0, 0]]  # max in x
    # a2 = [[0, 0, 0], [0, 1, 0], [0, 0, 0]]  # max in y
    a2 = [[1 / np.sqrt(2), 0, 0], [0, 1 / np.sqrt(2), 0], [0, 0, 0]]  # x-y girdle
    # a2 = [[1 / np.sqrt(2), 0, 0], [0, 0, 0], [0, 0, 1 / np.sqrt(2)]]  # x-z girdle
    # a2 = [[1. / 3., 0, 0], [0, 1. / 3., 0], [0, 0, 1. / 3.]]  # iso

    a2plot(np.matrix(a2))
