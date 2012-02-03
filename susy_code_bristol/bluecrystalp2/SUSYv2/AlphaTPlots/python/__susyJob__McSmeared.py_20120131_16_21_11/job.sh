#!/bin/sh
source /gpfs/cluster/phys/cl7359/susy_work/SUSYv2/setup.sh
cd /gpfs/cluster/phys/cl7359/susy_work/SUSYv2/AlphaTPlots/python
./McSmeared.py -j __susyJob__McSmeared.py_20120131_16_21_11/job.json -J ${SGE_TASK_ID} -S __susyJob__McSmeared.py_20120131_16_21_11/status/${SGE_TASK_ID}