#! /bin/bash
#
# backward_oop.sh
# Copyright (C) 2020 dal22 <dal22@fenrir.ess.washington.edu>
#
# Distributed under terms of the MIT license.
#

conv=$1
rc=$2
oop=$3
bm=$4

python prepare_backward_oop.py $conv $rc $oop $bm

# tmux split-window -b -l 2 "echo a_"$conv"_rc"$rc"_oop"$oop"_bm"$bm".sif; ./run_and_log.sh a_"$conv"_rc"$rc"_oop"$oop"_bm"$bm".sif"
# tmux select-pane -D
tmux split-window -b -l 2 "echo d_"$conv"_rc"$rc"_oop"$oop"_bm"$bm".sif; ./run_and_log.sh d_"$conv"_rc"$rc"_oop"$oop"_bm"$bm".sif"
tmux select-pane -D
