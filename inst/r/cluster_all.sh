#!/bin/bash

nohup ./cluster_stability.r &
echo "Stability done!"

nohup ./cluster_cooccurence.r &
echo "Cooccurence done!"

