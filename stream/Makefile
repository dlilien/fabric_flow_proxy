#! /bin/bash
#
# comp_bc.sh
# Copyright (C) 2019 dal22 <dal22@loki.ess.washington.edu>
#
# Distributed under terms of the MIT license.
#


all: read_routines.o SIA_BC.so input_data.so temp_read.so input_data_oop.so

clean:
	rm -f BCs.so input_data.so input_data_oop.so

%.so: %.f90
		elmerf90 -I. read_routines.o -o $@ $<

%.o: %.f90
		elmerf90 -o $@ $<

mesh: stream_ftw/mesh.nodes rstf/mesh.nodes

stream_ftw/mesh.nodes: stream_ftw.msh
	ElmerGrid 14 2 stream_ftw.msh -autoclean

rstf/mesh.nodes: rstf.msh
	ElmerGrid 14 2 rstf.msh -autoclean
	python shrink_tall.py 5.0 rstf

%.msh: %.geo
	gmsh $< -o $@ -1 -2 8
