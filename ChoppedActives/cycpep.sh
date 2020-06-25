#!/bin/bash

NATIVE=$1
SEQUENCE=$2
DESIGN=$3
ID=$4


echo /project2/mpcs56420/rosetta_bin_linux_2020.08.61146_bundle/main/source/bin/simple_cycpep_predict.static.linuxgccrelease \
-in:file:native $NATIVE \
-cyclic_peptide:sequence_file $SEQUENCE \
-cyclic_peptide:design_peptide false \
-cyclic_peptide:allowed_residues_by_position $DESIGN \
-cyclic_peptide:genkic_closure_attempts 1000 \
-cyclic_peptide:min_genkic_hbonds 2 -mute all \
-unmute protocols.cyclic_peptide_predict.SimpleCycpepPredictApplication \
-out:nstruct 250 \
-out:file:o design-$ID.pdb
