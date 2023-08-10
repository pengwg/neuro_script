#!/bin/bash

cores=10

# Absolute or relative path of the data folder to where the script located
data_path=FUS/
subject=sub-216-FUS
session=ses-00

# Choose one of the following reference type which corresponds to different sub-name_seeds_{ref_type}.csv files and reference nifti volumes.
# treatment: Subspot LPS coordinate in treatment volume
# AC or PC: LPS coordinates in the AC-PC-Midline coordinate system where AC or PC is (0,0,0)

ref_type='treatment'
# ref_type='AC'
# ref_type='PC'

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

printf "\n${GREEN}Entering $subject/$session/dwi/mrtrix/...$NC\n"

# Always use the reference volume and seed files from ses-00
if [ -d "../../../ses-00/anat" ]; then
    REF_nii=$(find ../../../ses-00/anat \( -name "${subject}_ses-00_$ref_type.nii" -o -name "${subject}_ses-00_$ref_type.nii.gz" \) | head -n 1)
fi

if [ -z "$REF_nii" ]; then
    echo -e "${YELLOW}${subject}_${session}_$ref_type.nii(.gz) not found.$NC"
    exit 1
fi

if ! [ -f "../../../ses-00/anat/${subject}_seeds_$ref_type.csv" ]; then
    echo -e "${YELLOW}${subject}_seeds_$ref_type.csv not found.$NC"
    exit 1
fi

# T1 volume in the dwi space
if ! [ -f "T1_coreg.nii.gz" ]; then
    echo -e "${YELLOW}T1_coreg.nii.gz not found.$NC"
    exit 1
fi

if ! [ -d "tracks_from_seeds" ]; then
    mkdir tracks_from_seeds
fi

cd tracks_from_seeds

# Register the ref volume to T1_coreg.nii.gz
# ${ref_type}2dwi_0GenericAffine.mat and ${ref_type}2dwi_Warped.nii.gz will be produced
if ! [ -f "${ref_type}2dwi_0GenericAffine.mat" ]; then
    antsRegistrationSyNQuick.sh -d 3 -t r -f ../T1_coreg.nii.gz -m "../$REF_nii" -o ${ref_type}2dwi_ -n $cores
fi 

# Apply the registration transformation to the seed points in the csv file.
# Note that ANTs software assume the coordinates in LPS system.
antsApplyTransformsToPoints -d 3 -i "../../../../ses-00/anat/${subject}_seeds_$ref_type.csv" -o "${subject}_seeds_to_dwi.csv" -t [${ref_type}2dwi_0GenericAffine.mat, 1]

exec 3< <(tail -n +2 "../../../../ses-00/anat/${subject}_seeds_$ref_type.csv")
exec 4< <(tail -n +2 "${subject}_seeds_to_dwi.csv")

# x0, y0, z0 are LPS coordinates in the reference volume
# x, y, z are transformed LPS coordinates in the dwi space
# r0 and r are radius of the seed sphere
while IFS=',' read -r x0 y0 z0 r0 <&3 && IFS=',' read -r x y z r <&4; do
    # Skip the nan lines (corresponding to blank lines in _seeds_$ref_type.csv) in the csv files
    if [ "$x" = "nan" ]; then
        continue
    fi
    
    # tckgen work with RAS coordinates. The following is conversion from LPS to RAS
    x=$(echo "$x * -1" | bc)
    y=$(echo "$y * -1" | bc)
    x0=$(echo "$x0 * -1" | bc)
    y0=$(echo "$y0 * -1" | bc)
    tckgen -act ../5tt_coreg.mif -backtrack -seed_sphere $x,$y,$z,$r -select 1k ../wmfod_norm.mif "tracks_1k_${ref_type}_${x0}_${y0}_${z0}_RAS_${x}_${y}_${z}.tck" -nthreads $cores -force
done

# Close the file descriptors
exec 3<&-
exec 4<&-

