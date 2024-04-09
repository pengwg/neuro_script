#!/bin/bash

check_registraion_only=1

session_path=~/Work/fusOUD/FUS/sub-220-FUS/ses-30

if ! [ -d $session_path ]; then
    echo "Path not exists: $session_path"
    exit 1
fi

IFS='/' read -ra parts <<< $session_path
N=${#parts[@]}
sub_name="${parts[N-2]}_${parts[N-1]}"
    
cd $session_path/dwi

if ! [ -d mrtrix3 ]; then
    echo "$session_path/dwi/mrtrix3 folder not exists, run BATMAN script first!"
    exit 1
fi

cd mrtrix3

# Check freesurfer registration
mrview mean_b0_preprocessed.mif -overlay.load T1_FS_coreg.nii.gz -overlay.load aparc+aseg_coreg.nii.gz &

if [ $check_registraion_only ]
    exit 0
fi

# Check 3 PAs alignment
mapfile -t PA_files < <(find . -type f -name \*2mm_PA\*)
mrview ${PA_files[0]} -overlay.load ${PA_files[1]} -overlay.load ${PA_files[2]} &

# Check bias correction
mrview ${sub_name}_den_unr_preproc_unbiased.mif -overlay.load bias.mif &

# Check response function(s) for spherical deconvolution
shview wm.txt &

# Check volume fraction
mrview ${sub_name}_den_unr_preproc_unbiased.mif -overlay.load vf.mif &
mrview vf.mif -odf.load_sh wmfod.mif &

# Check 5tt volume
mrview T1_FS_coreg.nii.gz -overlay.load 5tt_coreg.mif &



