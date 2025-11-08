#!/bin/bash

num_threads=10
data_path="/home/pw0032/Data/FUS-RCT"

sub_name=sub-018-RCT
session=ses-2-00

T1_VOL=$(ls $data_path/$sub_name/$session/anat/${sub_name}_${session}_T1*.nii.gz | head -n 1)

recon-all -all -parallel -openmp $num_threads -i $T1_VOL -s ${sub_name}_$session
