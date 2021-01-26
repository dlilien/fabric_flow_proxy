#! /bin/bash
#
# pull_vars.sh
# Copyright (C) 2020 dal22 <dal22@loki.ess.washington.edu>
#
# Distributed under terms of the MIT license.
#

nodewisedata --vars 'age,temp,temp homologous,fabric 1,fabric 2,fabric 3' --reg b,b,b,b,b,b -full stream_ftw/s_ftw_rc2.0e-3_oop0.0_bm0.0_????.vtu


nodewisedata --vars 'age,temp,temp homologous,fabric 1,fabric 2,fabric 3' --reg b,b,b,b,b,b -full rstf/deform2e-3_????.vtu

exit
# In round 1, we recover the results from the first relaxation
nodewisedata --vars 'age,temp,temp homologous' --reg b,b,b -full relaxed_stream_tall/deform????.vtu
# I think we can just do these 3 since 4 and 5 should be zero
nodewisedata --vars 'fabric 1' -full -xavg relaxed_stream_tall/deform????.vtu
nodewisedata --vars 'fabric 2' -full -xavg relaxed_stream_tall/deform????.vtu
nodewisedata --vars 'fabric 3' -full -xavg relaxed_stream_tall/deform????.vtu

nodewisedata --vars 'age,temp,temp homologous,fabric 1,fabric 2,fabric 3' --reg b,b,b,b,b,b -full relaxed_stream_tall_fabric/deform_diffuse_rr????.vtu

