#!/bin/bash

cores=10

# Absolute or relative path of the data folder to where the script located
data_path=FUS/
subject=sub-215-FUS
session=ses-00

num_tracks=100k

YELLOW='\033[0;33m'
GREEN='\033[0;32m'
NC='\033[0m'

cd $(dirname %0)
if [ -d "$data_path/$subject/$session/dwi/mrtrix" ]; then
    cd $data_path/$subject/$session/dwi/mrtrix
else
    echo -e "${YELLOW}$subject/$session run BATMAN script first.$NC"
    exit 1
fi

printf "\n${GREEN}Entering $subject/$session/dwi/mrtrix...$NC\n"

# Alway use treatment volume from ses-00 folder 
if [ -d "../../../ses-00/anat" ]; then
    REF_nii=$(find ../../../ses-00/anat \( -name "${subject}_ses-00_treatment.nii" -o -name "${subject}_ses-00_treatment.nii.gz" \) | head -n 1)
fi

if [ -z "$REF_nii" ]; then
    echo -e "${YELLOW}${subject}_ses-00_treatment.nii(.gz) not found.$NC"
    exit 1
fi

# Alwasy use mask volume from ses-00 folder
if ! [ -f "../../../ses-00/anat/${subject}_ses-00_mask.nii.gz" ]; then
    echo -e "${YELLOW}${subject}_ses-00_mask.nii.gz not found.$NC"
    exit 1
fi

# Coregistered anatomic volume from dwi process
if ! [ -f "T1_coreg.nii.gz" ]; then
    echo -e "${YELLOW}T1_coreg.nii.gz not found.$NC"
    exit 1
fi

if ! [ -d "tracks_from_mask" ]; then
    mkdir tracks_from_mask
fi
cd tracks_from_mask

# Generate transform matrix from treatment to dwi registration
if ! [ -f "treatment2dwi_0GenericAffine.mat" ]; then
    antsRegistrationSyNQuick.sh -d 3 -t r -f ../T1_coreg.nii.gz -m "../$REF_nii" -o treatment2dwi_ -n $cores
    
    # Check registration
    mrview treatment2dwi_Warped.nii.gz -overlay.load ../fs_parcels_coreg.mif -overlay.opacity 0.5 -mode 2 &
fi

# Apply the registration transform to the mask volume
antsApplyTransforms -d 3 -i "../../../../ses-00/anat/${subject}_ses-00_mask.nii.gz" -o mask_to_dwi.nii.gz -r ../T1_coreg.nii.gz -t treatment2dwi_0GenericAffine.mat
mrconvert mask_to_dwi.nii.gz mask_to_dwi.mif -force

# Generating tracks from the registered mask volume as seed
if ! [ -f "tracks_${num_tracks}_from_mask.tck" ]; then
    tckgen -act ../5tt_coreg.mif -backtrack \
           -seed_image mask_to_dwi.mif -select ${num_tracks} ../wmfod_norm.mif \
           tracks_${num_tracks}_from_mask.tck -nthreads $cores -force
fi

tcksift2 -act ../5tt_coreg.mif -out_mu sift_mu_${num_tracks}.txt -out_coeffs sift_coeffs_${num_tracks}.txt \
         tracks_${num_tracks}_from_mask.tck ../wmfod_norm.mif \
         sift_${num_tracks}_from_mask.txt -nthreads $cores -force

# Computing the histogram of tracks lenghth
tckstats tracks_${num_tracks}_from_mask.tck -histogram tracks_length_${num_tracks}_hist.csv -dump tracks_length_${num_tracks}.csv -force

# Computing fractional anisotropy
dwi2tensor ../${subject}_${session}_den_unr_preproc_unbiased.mif tensor.mif -force -nthreads $cores
tensor2metric tensor.mif -fa FA.mif -force -nthreads $cores

# Computing the mean FA of tracks
tcksample -stat_tck mean tracks_${num_tracks}_from_mask.tck FA.mif tracks_meanFA_${num_tracks}.csv -force -nthreads $cores
 
# Computing mean FA connectome 
tck2connectome -symmetric -zero_diagonal \
               -tck_weights_in sift_${num_tracks}_from_mask.txt tracks_${num_tracks}_from_mask.tck \
               ../fs_parcels_coreg.mif mask_meanFA_${num_tracks}_connectome.csv \
               -scale_file tracks_meanFA_${num_tracks}.csv -stat_edge mean -force 

tck2connectome -tck_weights_in sift_${num_tracks}_from_mask.txt tracks_${num_tracks}_from_mask.tck \
               ../fs_parcels_coreg.mif mask_meanFA_${num_tracks}_vector_connectome.csv -vector \
               -scale_file tracks_meanFA_${num_tracks}.csv -stat_edge mean -force
