#!/bin/bash

cores=10

# Absolute or relative path of the data folder to where the script located
data_path=FUS/
subject=sub-216-FUS
session=ses-00

YELLOW='\033[0;33m'
GREEN='\033[0;32m'
NC='\033[0m'

cd $(dirname %0)
if [ -d "$data_path/$subject/$session/dwi/mrtrix" ]; then
    cd $data_path/$subject/$session/dwi/mrtrix
else
    echo -e "${YELLOW}$subject/$session run mrtrix first.$NC"
    exit 1
fi

printf "\n${GREEN}Entering $subject/$session/dwi/mrtrix...$NC\n"

chmod a+x *

if [ -d "../../../ses-00/anat" ]; then
    REF_nii=$(find ../../../ses-00/anat \( -name "${subject}_ses-00_treatment.nii" -o -name "${subject}_ses-00_treatment.nii.gz" \) | head -n 1)
fi
    
if [ -z "$REF_nii" ]; then
    echo -e "${YELLOW}${subject}_ses-00_treatment.nii(.gz) not found.$NC"
    exit 1
fi

if ! [ -f "T1_coreg.nii.gz" ]; then
    echo -e "${YELLOW}T1_coreg.nii.gz not found.$NC"
    exit 1
fi

if ! [ -f "../../../ses-00/anat/${subject}_ses-00_mask.nii.gz" ]; then
    echo -e "${YELLOW}${subject}_ses-00_mask.nii.gz not found.$NC"
    exit 1
fi

if ! [ -d "tracks_from_mask" ]; then
    mkdir tracks_from_mask
fi
cd tracks_from_mask

if ! [ -f "treatment2dwi_0GenericAffine.mat" ]; then
    antsRegistrationSyNQuick.sh -d 3 -t r -f ../T1_coreg.nii.gz -m "../$REF_nii" -o treatment2dwi_ -n $cores
fi

antsApplyTransforms -d 3 -i "../../../../ses-00/anat/${subject}_ses-00_mask.nii.gz" -o mask_to_dwi.nii.gz -r ../T1_coreg.nii.gz -t treatment2dwi_0GenericAffine.mat
mrconvert mask_to_dwi.nii.gz mask_to_dwi.mif -force

tckgen -act ../5tt_coreg.mif -backtrack -seed_image mask_to_dwi.mif -select 100k ../wmfod_norm.mif tracks_100k_from_mask.tck -nthreads $cores -force

chmod a+x *
