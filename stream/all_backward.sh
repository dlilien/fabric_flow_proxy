#! /bin/bash
#
# all_backward.sh
# Copyright (C) 2020 dal22 <dal22@fenrir.ess.washington.edu>
#
# Distributed under terms of the MIT license.
#

conv=ftw
rc=2.0e-3
oop=0.0
for bm in 1.0e-3; do
    fn=s_"$conv"_rc"$rc"_oop"$oop"_bm"$bm".sif
    tmux split-window -b -l 2 "./back_to_deform.sh $fn"
    tmux select-pane -D
done
exit


conv=ftw
oop=0.0
bm=0.0
for rc in 2.0 2.0e-1; do
    fn=s_"$conv"_rc"$rc"_oop"$oop"_bm"$bm".sif
    tmux split-window -b -l 2 "./back_to_deform.sh $fn"
    tmux select-pane -D
done
exit
rc=2.0e-3
oop=0.0
bm=0.0
for conv in uc hiacc; do
    fn=s_"$conv"_rc"$rc"_oop"$oop"_bm"$bm".sif
    tmux split-window -b -l 2 "./back_to_deform.sh $fn"
    tmux select-pane -D
done
exit

conv=ftw
rc=2.0e-3
oop=0.0
for bm in 0.0 1.0e-2 2.0e-2; do
    fn=s_"$conv"_rc"$rc"_oop"$oop"_bm"$bm".sif
    tmux split-window -b -l 2 "./back_to_deform.sh $fn"
    tmux select-pane -D
done


conv=ftw
oop=0.0
bm=0.0
for rc in 2.0e-4 2.0e-2; do
    fn=s_"$conv"_rc"$rc"_oop"$oop"_bm"$bm".sif
    tmux split-window -b -l 2 "./back_to_deform.sh $fn"
    tmux select-pane -D
done

