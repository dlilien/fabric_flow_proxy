#! /bin/bash
#
# run_and_log.sh
# Copyright (C) 2020 dal22 <dal22@fenrir.ess.washington.edu>
#
# Distributed under terms of the MIT license.
#


sif=$1
pref=${sif%.*}

out=logs/"$pref".log
if [[ `hostname` =~ sleipnir* ]]; then
    # Avoid writing to a RAID on a different machine
    out=/home/dal22/logs/"$pref".log
fi

ElmerSolver $sif > $out
