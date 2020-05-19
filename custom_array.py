import sys
import os
import time
import platform
import random
import subprocess

# How many nodes are we going to distribute the
# work to?
nodes = 1

out  = subprocess.run(['split', '-l', '2', 'd', 'test_20.fasta', 'peptide_'], \
                      capture_output=True, text=True)

exit()


# Create a custom sbatch file for this run
sbatch="""#!/bin/bash
#SBATCH --output=%%A-%%j-%%a.out
#SBATCH --error=%%A-%%j-%%a.err
#SBATCH --job-name=custom_array_job
#SBATCH --ntasks=1
#SBATCH --partition=broadwl

#SBATCH --array=0-%s

date
hostname
pwd
echo "$SLURM_ARRAY_JOB_ID $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"
echo "Run you script with my_custom_script.py temp_%s"
# Here we would run our program with the database temp_??
sleep 5

Scwrl4 -i temp_%%s -s peptide_%%s -o temp_%%s.pdb

date
""" % (chunks,i)

# Write the sbatch string to a file
sbatch_file = "custom_array.sbatch"
with open(sbatch_file, "w") as outfile:
    outfile.write(sbatch)
exit()

# Submit the file using sbatch and collect the 
# job id from the output
out  = subprocess.run(['sbatch', '-A', 'mpcs56420', sbatch_file], capture_output=True, text=True)
job_id = out.stdout.strip("Submitted batch job ")
job_id = job_id.strip()
print("Main Job: " + job_id)

# Create a folder to hold all our files prefixed  with the job id
workspace_folder = job_id + "-workspace"

# Create a sbatch file for a "cleanup" job that will move all the output files
sbatch="""#!/bin/bash
#SBATCH --output=%%j-cleanup.out
#SBATCH --error=%%j-cleanup.err
#SBATCH --job-name=custom_array_job_cleanup
#SBATCH --ntasks=1
#SBATCH --partition=broadwl
#SBATCH --dependency=afterok:%s

date
hostname
pwd

echo "Cleanup"
sleep 5
mkdir %s
mv %s-* %s/
mv custom_array.sbatch temp_* %s/


date
""" % (job_id, workspace_folder, job_id, workspace_folder, workspace_folder)

# Write it to a file
sbatch_file = "cleanup.%s.sbatch" % job_id
with open(sbatch_file, "w") as outfile:
    outfile.write(sbatch)

# Submit it
out  = subprocess.run(['sbatch', '-A', 'mpcs56420', sbatch_file], capture_output=True, text=True)
#print(out)
job_id = out.stdout.strip("Submitted batch job ")
job_id = job_id.strip()
print("Cleanup: %s" % job_id)

