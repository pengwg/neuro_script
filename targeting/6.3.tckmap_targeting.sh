#!/bin/bash

data_path=/mnt/evo/FUS-NV

r=3.5
num_tracks=6k
num_tracks_per_side=3k
postfix=reduced

threads=12
rerun=0

subject=sub-021-NAV

dwi_path=$(find $data_path/$subject -type d -name dwi | head -n 1)

T1_aligned=$(find ${dwi_path}/../anat -name "*_aligned.nii.gz")
T1_dwi=${dwi_path}/mrtrix3/mean_b0_preprocessed.nii.gz

if ! [[ -f "$T1_dwi" ]]; then
    echo "$subject run BATMAN first!"
    continue
fi   

out_path=${dwi_path}/../targeting

mkdir -p $out_path

# Register target volume to DWI       
if ! [[ -f "$out_path/T1w_aligned_warp2dwi.lta" ]]; then
    echo "mri_synthmorph register -m rigid -t $out_path/T1w_aligned_warp2dwi.lta -T $out_path/T1w_aligned_warp2dwi.inv.lta \
        -o $out_path/T1w_aligned_morphed2dwi.nii.gz $T1_aligned  $T1_dwi -g -j $threads"
    mri_synthmorph register -m rigid -t $out_path/T1w_aligned_warp2dwi.lta -T $out_path/T1w_aligned_warp2dwi.inv.lta \
        -o $out_path/T1w_aligned_morphed2dwi.nii.gz $T1_aligned  $T1_dwi -g -j $threads
    lta_convert --inlta $out_path/T1w_aligned_warp2dwi.inv.lta --outitk $out_path/T1w_aligned_warp2dwi_itk.inv.txt
fi

# Generate tracks and tckmap from both left and right seeds    
tracks_file="$out_path/target2moFC_r${r}_tckedit_$postfix.tck"
weight_file="$out_path/tckedit_r${r}.txt"

if [[ ! -f "$tracks_file" ]] || (( rerun )); then
    # Transform target coordinates to DWI space
    antsApplyTransformsToPoints -d 3 -i ${dwi_path}/../anat/targets_LPS.csv -o $out_path/targets_dwi_LPS.csv -t $out_path/T1w_aligned_warp2dwi_itk.inv.txt
    matlab -batch "gen_targets_mask('$out_path', ${r})"
    
    tckedit_cmd="tckedit $dwi_path/mrtrix3/tracks_10M.tck $tracks_file \
        -include $out_path/targets_mask.nii.gz \
        -include ${dwi_path}/../anat/include_mofc_reduced.nii.gz \
        -tck_weights_in $dwi_path/mrtrix3/sift2_10M.txt \
        -tck_weights_out $weight_file \
        -nthreads $threads -force"
        
    echo $tckedit_cmd
    $tckedit_cmd
                
    tckmap $tracks_file $out_path/TDI_target2moFC_r${r}_sift2_bilateral_$postfix.nii.gz -template $T1_dwi -tck_weights_in $weight_file -force -nthreads $threads
    # mrconvert $out_path/TDI_target2moFC_seedr${r}_sift2_bilateral.mif $out_path/TDI_target2moFC_seedr${r}_sift2_bilateral.nii.gz -datatype float32 -force
fi
     
# Transform tckmaps into aligned T1 space

if [[ ! -f "$out_path/TDI_target2moFC_aligned_r${r}_sift2_bilateral_$postfix.nii.gz" ]] || (( rerun )); then
    mri_vol2vol --mov $out_path/TDI_target2moFC_r${r}_sift2_bilateral_$postfix.nii.gz \
        --o $out_path/TDI_target2moFC_aligned_r${r}_sift2_bilateral_$postfix.nii.gz --lta $out_path/T1w_aligned_warp2dwi.inv.lta --keep-precision
fi

