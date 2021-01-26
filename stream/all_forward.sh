#! /bin/bash
#
# all_forward.sh
# Copyright (C) 2020 dal22 <dal22@fenrir.ess.washington.edu>
#
# Distributed under terms of the MIT license.
#


for oop in 1.0e-3 1.0e-2 1.0e-1; do
    tmux split-window -b -l 2 "./forward.sh ftw 2.0e-3 $oop 0.0"
    tmux select-pane -D
done

exit
for conv in ftw uc hiacc; do
    tmux split-window -b -l 2 "./forward.sh $conv 2.0e-3 0.0 0.0"
    tmux select-pane -D
done

for rc in 2.0e-4 2.0e-2; do
    tmux split-window -b -l 2 "./forward.sh ftw $rc 0.0 0.0"
    tmux select-pane -D
done

for bm in 1.0e-2 2.0e-2; do
    tmux split-window -b -l 2 "./forward.sh ftw 2.0e-3 0.0 $bm"
    tmux select-pane -D
done


for bm in 1.0e-3; do
    tmux split-window -b -l 2 "./forward.sh ftw 2.0e-3 0.0 $bm"
    tmux select-pane -D
done

for rc in 2.0e-1 2.0; do
    tmux split-window -b -l 2 "./forward.sh ftw $rc 0.0 0.0"
    tmux select-pane -D
done

for conv in ftwlc nuc nuclc midacc; do
    tmux split-window -b -l 2 "./forward.sh $conv 2.0e-3 0.0 0.0"
    tmux select-pane -D
done
