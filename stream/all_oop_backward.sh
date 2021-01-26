#! /bin/bash
#
# oop_forward.sh
# Copyright (C) 2020 dal22 <dal22@fenrir.ess.washington.edu>
#
# Distributed under terms of the MIT license.
#

for oop in 1.0e-3 1.0e-2 1.0e-1; do
    ./backward_oop.sh ftw 2.0e-3 $oop 0.0
done

