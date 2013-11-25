#! /usr/bin/env fish

set OUTPUT "output/"
set EXPERIMENT  $argv[1]
set STDOUT $OUTPUT$EXPERIMENT

if not mkdir $STDOUT
    exit 1
end


for max_nodes in (seq $argv[2] $argv[3])
  qsub utils/submit_mpi.sh -l nodes=$max_nodes:ppn=8 -e $STDOUT -o $STDOUT
end

