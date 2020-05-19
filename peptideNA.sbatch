#!/bin/bash
#SBATCH --output=%A-%j-%a.out
#SBATCH --error=%A-%j-%a.err
#SBATCH --job-name=peptide_generation
#SBATCH --ntasks=1
#SBATCH --partition=broadwl
#SBATCH --array=1-500


JOB_LABEL=$1
echo "Submittin job with lable: "$JOB_LABEL
JOB_DIR=PEPTIDE-$JOB_LABEL

# Database name
DB="test_1000.fasta"

date
hostname
pwd

PADDED_NUMBER=$(printf %04d $SLURM_ARRAY_TASK_ID)

# Make a temp directory
mkdir $JOB_DIR
split -l 2 -a 4 --additional-suffix=".fasta" -d $DB $JOB_DIR/peptide_
Scwrl4 -i template.pdb -s $JOB_DIR/peptide_$PADDED_NUMBER.fasta -o $JOB_DIR/peptide_$PADDED_NUMBER.pdb
python filter.py $JOB_DIR/peptide_$PADDED_NUMBER.pdb $PADDED_NUMBER | grep ">" > $JOB_DIR/peptide_$PADDED_NUMBER.metrics.txt

# Cleanup
mv $SLURM_JOB_ID* $JOB_DIR
cp $DB $JOB_DIR





#python
echo "$SLURM_ARRAY_JOB_ID $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID"
echo "Run you script with my_custom_script.py temp_9"
sleep 5

date
