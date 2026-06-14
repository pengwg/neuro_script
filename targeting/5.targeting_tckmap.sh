#!/bin/bash

threads=12

# Absolute or relative path of the data folder to where the script located
data_path=/home/pw0032/Data/FUS-RCT
subject=sub-018-RCT
session=ses-2-00

num_tracks=300k

mrtrix_path=$data_path/$subject/$session/dwi/mrtrix3

cd $(dirname %0)
basedir=$(pwd)

cd $mrtrix_path/tckmap

T1_AC=$(ls $data_path/$subject/$session/anat/${subject}_${session}_AC.nii.gz | head -n 1)
antsRegistrationSyNQuick.sh -d 3 -t r -f $T1_AC -m $mrtrix_path/T1_FS_coreg.nii.gz -o dwi2AC_ -n $threads

# Map from global tracks
# tckmap $in_path/tracks_10M.tck $out_path/TDI_weighted.mif -vox 1.0 -tck_weights_in $in_path/sift_10M.txt -nthreads $threads
# tckmap $in_path/tracks_10M.tck $out_path/endpoint_density_weighted.mif -vox 1.0 -ends_only -tck_weights_in $in_path/sift_10M.txt -nthreads $threads

# Map from tracks generated around NAc and selective regions
tckmap sphere_10_select_${num_tracks}_NAc.tck TDI_NAc.mif -vox 1.0 -force -nthreads $threads
mrconvert TDI_NAc.mif TDI_NAc.nii.gz -force
matlab -batch "addpath('$basedir'); apply_rigid_transform('TDI_NAc.nii.gz', 'TDI_NAc_AC.nii.gz', 'dwi2AC_0GenericAffine.mat')"

tckmap sphere_10_select_${num_tracks}_NAc.tck endpoint_density_NAc.mif -vox 1.0 -ends_only -force -nthreads $threads
mrconvert endpoint_density_NAc.mif endpoint_density_NAc.nii.gz -force
matlab -batch "addpath('$basedir'); apply_rigid_transform('endpoint_density_NAc.nii.gz', 'endpoint_density_NAc_AC.nii.gz', 'dwi2AC_0GenericAffine.mat')"


# Using sift2 weights
tcksift2 -act $mrtrix_path/5tt_coreg_hsvs.mif -out_mu sift2_mu.txt -out_coeffs sift2_coeffs.txt \
    sphere_10_select_${num_tracks}_NAc.tck $mrtrix_path/wmfod_norm.mif sift2_${num_tracks}.txt -force -nthreads $threads
         
tckmap sphere_10_select_${num_tracks}_NAc.tck TDI_NAc_weighted.mif -vox 1.0 -tck_weights_in sift2_${num_tracks}.txt -force -nthreads $threads
mrconvert TDI_NAc_weighted.mif TDI_NAc_weighted.nii.gz -force
matlab -batch "addpath('$basedir'); apply_rigid_transform('TDI_NAc_weighted.nii.gz', 'TDI_NAc_weighted_AC.nii.gz', 'dwi2AC_0GenericAffine.mat')"

tckmap sphere_10_select_${num_tracks}_NAc.tck endpoint_density_NAc_weighted.mif -vox 1.0 -ends_only -tck_weights_in sift2_${num_tracks}.txt -force -nthreads $threads
mrconvert endpoint_density_NAc_weighted.mif endpoint_density_NAc_weighted.nii.gz -force
matlab -batch "addpath('$basedir'); apply_rigid_transform('endpoint_density_NAc_weighted.nii.gz', 'endpoint_density_NAc_weighted_AC.nii.gz', 'dwi2AC_0GenericAffine.mat')"

        


