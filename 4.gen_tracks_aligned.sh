#!/bin/bash

data_path=/mnt/evo/FUS-NV
SUBJECTS_DIR=/mnt/evo/FS

num_tracks=10M

threads=12
rerun=0

sub_nums=(004 006 008 010 011 016 017 018 019 020 021)

for num in "${sub_nums[@]}"; do
    sub_name=sub-OUDNV0$num

    dwi_path=$(find $data_path/${sub_name} -type d -name dwi | head -n 1)
    T1_aligned=$(find ${dwi_path}/../anat -name "*-acpc.nii*")
    
    if [[ ! -f "${dwi_path}/mrtrix3/mean_b0_preprocessed.nii.gz" ]]; then
        echo "${sub_name} run BATMAN first!"
        exit
    fi   

    cd ${dwi_path}/mrtrix3

    if [[ ! -f "T1_aligned_brain.nii.gz" ]]; then
        mri_synthstrip -i $T1_aligned -o T1_aligned_brain.nii.gz -t $threads
    fi

    # Register DWI volume to AC-PC aligned volume       
    if [[ ! -f "wmfod_norm_aligned.mif" ]]; then
        mri_synthmorph register -m rigid -t dwi2aligned.lta -T dwi2aligned.inv.lta mean_b0_preprocessed.nii.gz T1_aligned_brain.nii.gz -g -j $threads
        lta_convert --inlta dwi2aligned.lta --outitk dwi2aligned.itk.txt
        transformconvert dwi2aligned.itk.txt itk_import dwi2aligned.matrix.txt -force
        mrtransform wmfod_norm.mif -linear dwi2aligned.matrix.txt wmfod_norm_aligned.mif -reorient_fod yes -force
        mrtransform mean_b0_preprocessed.mif -linear dwi2aligned.matrix.txt mean_b0_preprocessed_aligned.mif -force
    fi

    # Register freesurfer volume to AC-PC aligned
    if [[ ! -f "5tt_hsvs_aligned.mif" ]]; then       
        mri_synthmorph register -m rigid -t FS2aligned.lta -T FS2aligned.inv.lta T1_FS.nii.gz T1_aligned_brain.nii.gz -g -j $threads
        lta_convert --inlta FS2aligned.lta --outitk FS2aligned.itk.txt
        transformconvert FS2aligned.itk.txt itk_import FS2aligned.matrix.txt
        mrtransform 5tt_hsvs.mif -linear FS2aligned.matrix.txt 5tt_hsvs_aligned.mif -force
    fi

    # Tracks generation
    if [[ ! -f "tracks_${num_tracks}_aligned.tck" ]]; then
        5tt2gmwmi 5tt_hsvs_aligned.mif gmwmSeed_aligned.mif -force
        tckgen -act 5tt_hsvs_aligned.mif -backtrack -seed_gmwmi gmwmSeed_aligned.mif -select ${num_tracks} \
            wmfod_norm_aligned.mif tracks_${num_tracks}_aligned.tck -cutoff 0.08 -nthreads $threads
            #-crop_at_gmwmi
        
        # Clean up and exit if user presses Ctrl->C
        if ! [[ $? -eq 0 ]]  ; then
            rm tracks_${num_tracks}_aligned.tck
            exit 1
        fi
    fi

    if [[ ! -f "sift2_${num_tracks}_aligned.txt" ]]; then
        tcksift2 -act 5tt_hsvs_aligned.mif tracks_${num_tracks}_aligned.tck wmfod_norm_aligned.mif sift2_${num_tracks}_aligned.txt -nthreads $threads
    fi
done
