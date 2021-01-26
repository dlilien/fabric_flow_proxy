#! /bin/bash
#
# comp_bc.sh
# Copyright (C) 2019 dal22 <dal22@loki.ess.washington.edu>
#
# Distributed under terms of the MIT license.
#


all: SIA_BC.so input_data.so temp_read.so input_data_oop.so

clean:
	rm -f BCs.so input_data.so input_data_oop.so

%.so: %.f90
		elmerf90 -I/home/dal22/sw/ElmerReplacements /home/dal22/sw/ElmerReplacements/read_routines.o -o $@ $<

mesh: stream/mesh.nodes relaxed_stream_tall/mesh.nodes relaxed_stream_tall_fabric/mesh.nodes  rstf/mesh.nodes

stream_ftw: stream_ftw/mesh.nodes
	echo stream_ftw

stream_ftw/mesh.nodes: stream_ftw.msh
	ElmerGrid 14 2 stream_ftw.msh -autoclean
	python shrink_tall.py 5.0 stream_ftw

stream/mesh.nodes: stream.msh
	ElmerGrid 14 2 stream.msh -autoclean
	ElmerGrid 2 2 stream -metis 5 -halo
	ElmerGrid 2 2 stream -metis 8 -halo

relaxed_stream_tall/mesh.nodes: relaxed_stream_tall.msh
	ElmerGrid 14 2 relaxed_stream_tall.msh -autoclean
	python shrink_tall.py 10.0 relaxed_stream_tall
	ElmerGrid 2 2 relaxed_stream_tall -metis 5 -halo
	ElmerGrid 2 2 relaxed_stream_tall -metis 8 -halo

relaxed_stream_tall_fabric/mesh.nodes: relaxed_stream_tall_fabric.msh
	ElmerGrid 14 2 relaxed_stream_tall_fabric.msh -autoclean
	python shrink_tall.py 5.0 relaxed_stream_tall_fabric

rstf/mesh.nodes: rstf.msh
	ElmerGrid 14 2 rstf.msh -autoclean
	python shrink_tall.py 5.0 rstf

%.msh: %.geo
	gmsh $< -o $@ -1 -2 8

relaxed_stream_tall.geo: stream/deform_merged0076.vtu
	python recover_relaxed_surf.py 10.0 stream

relaxed_stream_tall_fabric.geo: relaxed_stream_tall/deform0151.vtu
	python recover_relaxed_surf.py 5.0 relaxed_stream_tall

rstf.geo: rstf/deform2e-3_0102.vtu
	python recover_relaxed_surf.py 5.0 rstf2

stream_ftw.geo: rstf/s_ftw_rc2.0e-3_oop0.0_bm0.0_1002.vtu
	python recover_relaxed_surf.py 5.0 stream_ftw
