#!/bin/bash

./cluster_stability.r

wait

echo "Stability done!"

./cluster_cooccurence.r

wait

echo "Cooccurence done!"

# Kill multiprocess generate by nohup:
# https://unix.stackexchange.com/questions/279642/why-i-couldnt-kill-the-nohup-process
# ps aux | grep r | grep ./cluster_stability.r # This will show you the process you're looking for
# kill $(ps aux | grep r | grep ./cluster_stability.r | awk '{print $2}') # This will kill it.
