#!/bin/bash

data_path=/mnt/evo/FUS-RCT

r=3
postfix=reduced

threads=12
rerun=0

subject=sub-025-RCT

dwi_path=$(find $data_path/$subject -type d -name dwi | head -n 1)   

out_path=${dwi_path}/../targeting

mkdir -p $out_path

if [[ ! -f "$dwi_path/mrtrix3/sift2_10M_aligned.txt" ]]; then
    tcksift2 -act $dwi_path/mrtrix3/5tt_hsvs_aligned.mif $dwi_path/mrtrix3/tracks_10M_aligned.tck $dwi_path/mrtrix3/wmfod_norm_aligned.mif $dwi_path/mrtrix3/sift2_10M_aligned.txt -nthreads $threads
fi

# Generate tracks and tckmap from both left and right seeds    
tracks_file="$out_path/target2moFC_r${r}_aligned_$postfix.tck"
weight_file="$out_path/tckedit_r${r}_aligned.txt"
include_spec="$out_path/include_spec_aligned_reduced.nii.gz"

matlab -batch "gen_targets_mask_aligned('$out_path', ${r})"

if [[ ! -f "$tracks_file" ]] || (( rerun )); then
    tckedit_cmd="tckedit $dwi_path/mrtrix3/tracks_10M_aligned.tck $tracks_file \
        -include $out_path/targets_mask_aligned.nii.gz \
        -include $include_spec \
        -tck_weights_in $dwi_path/mrtrix3/sift2_10M_aligned.txt \
        -tck_weights_out $weight_file \
        -nthreads $threads -force"
        
    echo $tckedit_cmd
    $tckedit_cmd      
    # mrconvert $out_path/TDI_target2moFC_seedr${r}_sift2_bilateral.mif $out_path/TDI_target2moFC_seedr${r}_sift2_bilateral.nii.gz -datatype float32 -force
fi

tckmap $tracks_file $out_path/TDI_target2moFC_r${r}_sift2_bilateral_aligned.nii.gz -template $out_path/include_spec_aligned.nii.gz -tck_weights_in $weight_file -force -nthreads $threads

