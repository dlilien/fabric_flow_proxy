#! /bin/bash
#
# back_to_deform.sh
# Copyright (C) 2020 dal22 <dal22@fenrir.ess.washington.edu>
#
# Distributed under terms of the MIT license.
#

input_sif=$1
name=${input_sif%.*}
name=${name#*_}

python prepare_reverse.py $name
echo d_$name
./run_and_log.sh d_"$name".sif
