#!/bin/bash

cores=10

# Absolute or relative path of the data folder to where the script located
data_path=FUS/
subject=sub-219-FUS
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


if [ -d "../../anat" ]; then
    T1_ACPC_nii=$(find ../../anat \( -name "*T1w_ACPC.nii" -o -name "*T1w_ACPC.nii.gz" \) | head -n 1)
fi
    
if [ -z "$T1_ACPC_nii" ]; then
    echo -e "${YELLOW}AC-PC aligned T1 volume not found.$NC"
    exit 1
fi

if ! [ -f "mean_b0_preprocessed.nii.gz" ]; then
    echo -e "${YELLOW}mean_b0_preprocessed.nii.gz not found.$NC"
    exit 1
fi

if ! [ -f "ACPCtodwi_0GenericAffine.mat" ]; then
    antsRegistrationSyNQuick.sh -d 3 -t r -f mean_b0_preprocessed.nii.gz -m "$T1_ACPC_nii" -o ACPCtodwi_
    antsApplyTransforms -d 3 -i "$T1_ACPC_nii"  -o T1_ACPC_coreg.nii.gz -r "$T1_ACPC_nii" -t ACPCtodwi_0GenericAffine.mat
fi

if ! [ -f "../../anat/${subject}_seeds.csv" ]; then
    echo -e "${YELLOW}${subject}_seeds.csv not found.$NC"
    exit 1
fi

antsApplyTransformsToPoints -d 3 -i "../../anat/${subject}_seeds.csv" -o "${subject}_seeds_to_dwi.csv" -t [ACPCtodwi_0GenericAffine.mat, 1]

exec 3< <(tail -n +2 "../../anat/${subject}_seeds.csv")
exec 4< <(tail -n +2 "${subject}_seeds_to_dwi.csv")

while IFS=',' read -r x0 y0 z0 r0 label comment <&3 && IFS=',' read -r x y z r label2 comment2 <&4; do
    x=$(echo "$x * -1" | bc)
    y=$(echo "$y * -1" | bc)
    x0=$(echo "$x0 * -1" | bc)
    y0=$(echo "$y0 * -1" | bc)
    tckgen -act 5tt_coreg.mif -backtrack -seed_sphere $x,$y,$z,$r -select 1k wmfod_norm.mif "tracks_1k_${x0}_${y0}_${z0}.tck" -nthreads $cores -force
done

# Close the file descriptors
exec 3<&-
exec 4<&-
