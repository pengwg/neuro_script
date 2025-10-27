#!/bin/bash

threads=12

# Absolute or relative path of the data folder to where the script located
data_path=/mnt/msi/Data/naviFUS
subject=sub-009-NAVI
session=ses-00

num_tracks=300k

in_path=$data_path/$subject/$session/dwi/mrtrix3
out_path=$data_path/$subject/$session/anat

# Map from global tracks
# tckmap $in_path/tracks_10M.tck $out_path/TDI_weighted.mif -vox 1.0 -tck_weights_in $in_path/sift_10M.txt -nthreads $threads
# tckmap $in_path/tracks_10M.tck $out_path/endpoint_density_weighted.mif -vox 1.0 -ends_only -tck_weights_in $in_path/sift_10M.txt -nthreads $threads

# Map from tracks generated around NAc and selective regions
tckmap $in_path/targeting_tracks_AC/sphere_10_select_${num_tracks}_NAc.tck $out_path/TDI_NAc.mif -vox 1.0 -nthreads $threads
tckmap $in_path/targeting_tracks_AC/sphere_10_select_${num_tracks}_NAc.tck $out_path/endpoint_density_NAc.mif -vox 1.0 -ends_only -nthreads $threads

   # Using sift2 weights
tcksift2 -act $in_path/5tt_coreg_hsvs.mif -out_mu $out_path/sift_mu.txt -out_coeffs $out_path/sift_coeffs.txt \
         $in_path/targeting_tracks_AC/sphere_10_select_${num_tracks}_NAc.tck $in_path/wmfod_norm.mif $out_path/sift2_${num_tracks}.txt -nthreads $threads
         
tckmap $in_path/targeting_tracks_AC/sphere_10_select_${num_tracks}_NAc.tck $out_path/TDI_NAc_weighted.mif -vox 1.0 -tck_weights_in $out_path/sift2_${num_tracks}.txt -nthreads $threads
tckmap $in_path/targeting_tracks_AC/sphere_10_select_${num_tracks}_NAc.tck $out_path/endpoint_density_NAc_weighted.mif -vox 1.0 -ends_only -tck_weights_in $out_path/sift2_${num_tracks}.txt -nthreads $threads
