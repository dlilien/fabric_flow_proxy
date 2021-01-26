#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 dlilien <dlilien@hozideh>
#
# Distributed under terms of the MIT license.

"""

"""

import vtk
from vtk.util.numpy_support import vtk_to_numpy
from matplotlib.tri import Triangulation, LinearTriInterpolator
import numpy as np


class VTU(object):
    def __init__(self, filename=None):
        """Creates a vtu object by reading the specified file."""
        if filename is None:
            self.ugrid = vtk.vtkUnstructuredGrid()
        else:
            self.gridreader = None
            if filename[-4:] == ".vtu":
                self.gridreader = vtk.vtkXMLUnstructuredGridReader()
            elif filename[-5:] == ".pvtu":
                self.gridreader = vtk.vtkXMLPUnstructuredGridReader()
            else:
                raise Exception("ERROR: don't recognise file extension" + filename)
            self.gridreader.SetFileName(filename)
            self.gridreader.Update()
            self.ugrid = self.gridreader.GetOutput()
            self.points = self.ugrid.GetPoints()
            self.pointdata = self.ugrid.GetPointData()
            self.cells = self.ugrid.GetCells()
            self.celldata = self.ugrid.GetCellData()
            if self.ugrid.GetNumberOfPoints() + self.ugrid.GetNumberOfCells() == 0:
                raise Exception("ERROR: No points or cells found after loading vtu " + filename)
        self._coords = None
        self.filename = filename
        self._data_dict = {}
        self._rawdata_dict = {}
        self._rawcelldata_dict = {}
        self._raw_coords = None

    @property
    def coords(self):
        if self._coords is None:
            self._coords, self._un_ind = self._get_coords()
        return self._coords

    @property
    def raw_coords(self):
        if self._raw_coords is None:
            data = vtk_to_numpy(self.ugrid.GetPoints().GetData()).round(decimals=8)
            ncols = data.shape[1]
            dtype = data.dtype.descr * ncols
            self._raw_coords = data.view(dtype)
        return self._raw_coords

    @coords.setter
    def coords(self, value):
        self._coords = value

    @property
    def un_ind(self):
        if self._un_ind is None:
            self._coords, self._un_ind = self._get_coords()
        return self._un_ind

    @un_ind.setter
    def un_ind(self, value):
        self._un_ind = value

    def _get_coords(self):
        data = vtk_to_numpy(self.ugrid.GetPoints().GetData()).round(decimals=1)
        ncols = data.shape[1]
        dtype = data.dtype.descr * ncols
        structr = data.view(dtype)
        c, ui = np.unique(structr, return_index=True)
        return c.view(data.dtype).reshape(-1, ncols), ui

    @property
    def data_dict(self):
        if self._data_dict == {}:
            self._data_dict = self._get_data_dict()
        return self._data_dict

    def _get_data_dict(self):
        return {self.pointdata.GetArrayName(i): vtk_to_numpy(self.pointdata.GetArray(i))[self.un_ind] for i in range(self.pointdata.GetNumberOfArrays())}

    @property
    def rawdata_dict(self):
        if self._rawdata_dict == {}:
            self._rawdata_dict = self._get_rawdata_dict()
        return self._rawdata_dict

    def _get_rawdata_dict(self):
        dd = {self.pointdata.GetArrayName(i): vtk_to_numpy(self.pointdata.GetArray(i)) for i in range(self.pointdata.GetNumberOfArrays())}
        for i, c in enumerate(['coordsX', 'coordsY', 'coordsZ']):
            dd[c] = self.raw_coords['f{:d}'.format(i)].flatten()
        return dd

    def _get_rawcelldata_dict(self):
        dd = {self.celldata.GetArrayName(i): vtk_to_numpy(self.celldata.GetArray(i)) for i in range(self.celldata.GetNumberOfArrays())}
        return dd

    @property
    def rawcelldata_dict(self):
        if self._rawcelldata_dict == {}:
            self._rawcelldata_dict = self._get_rawcelldata_dict()
        return self._rawcelldata_dict

    def get_geometry_pt_ids(self, geo_id, order=None):
        rel_cell_ids = np.argwhere(self.rawcelldata_dict['GeometryIds'] == geo_id).flatten()
        if len(rel_cell_ids) == 0:
            print('No Id {:d}'.format(geo_id))
            return None
        rel_pts = []
        for cell_id in rel_cell_ids:
            cell = self.ugrid.GetCell(cell_id)
            pids = cell.GetPointIds()
            for i in range(pids.GetNumberOfIds()):
                pid = pids.GetId(i)
                if pid not in rel_pts:
                    rel_pts.append(pid)
        pts = np.array(rel_pts, dtype=int)
        if order is not None:
            if order.upper() not in ['X', 'Y', 'Z']:
                raise ValueError('Order must be x, y or z')
            else:
                crds = self.rawdata_dict['coords' + order.upper()][pts]
                return pts[np.argsort(crds)]
        else:
            return pts

    def get_geometry_data(self, geo_id, order=None):
        pts = self.get_geometry_pt_ids(geo_id, order=order)
        return {key: val[pts] for key, val in self.rawdata_dict.items()}
