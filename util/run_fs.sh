#!/bin/bash

num_threads=16
data_path="/mnt/msi/Data/naviFUS"

sub_name=sub-009-NAVI
session=ses-00

T1_VOL=$(ls $data_path/$sub_name/$session/anat/${sub_name}_${session}_T1*.nii.gz | head -n 1)

recon-all -all -parallel -openmp $num_threads -i $T1_VOL -s $sub_name_$session
