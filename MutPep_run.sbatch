#!/bin/bash
#SBATCH --job-name=serial_job_test    # Job name
#SBATCH --mail-type=ALL               # Mail events - NONE, BEGIN, END, FAIL, ALL
#SBATCH --mail-user=caseys@uchicago.edu    # Where to send mail
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=1             # Number of threads per task
#SBATCH --mem=1gb                     # Job memory request
#SBATCH --time=00:15:00               # Time limit hrs:min:sec
#SBATCH --output=%j.serial_job.out
#SBATCH --error=%j.serial_job.err

echo "> Running plot script on a single CPU core"
pwd; hostname; date

echo "> Staring the run"
python3 MutatePeptides.py small_test.fasta 1,2,3 /scratch/midway2/caseys/SmallTestOutput
echo "> Done Running"

