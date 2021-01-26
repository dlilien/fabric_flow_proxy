#! /bin/bash
#
# back_to_deform.sh
# Copyright (C) 2020 dal22 <dal22@fenrir.ess.washington.edu>
#
# Distributed under terms of the MIT license.
#

conv=$1
rc=$2
oop=$3
bm=$4

python prepare_forward.py $conv $rc $oop $bm
echo s_"$conv"_rc"$rc"_oop"$oop"_bm"$bm".sif
./run_and_log.sh s_"$conv"_rc"$rc"_oop"$oop"_bm"$bm".sif
