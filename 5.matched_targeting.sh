#!/bin/bash

# This requires ${subject}_ses-00_$ref_type.nii(.gz) in the ses-00/anat folder for the selected participant. 

# Target tracks seeds are located in target_seeds_${seeds_tag}.csv in the main FUS folder. 
# An example of target_seeds_${seeds_tag}.csv file. Also note it should be in LPS (-R,-A,S) coordinate system.
# x,    y,  z,	r
# 8,   -2,	0,	2.5
# 8.5, -2,	0,	2.5

# It will create a directory called targeting_tracks_${seeds_tag} and sample all the possible combinations of coordinates from csv file stored in main FUS folder in order to review for targeting meetings. 

cores=15

# Absolute or relative path of the data folder to where the script located
data_path=/home/peng/Work/fusOUD/FUS-RCT/
subject=sub-006-RCT
session=ses-1-00

num_tracks=10M
radius=3

# Using target_seeds_${seeds_tag}.csv
seeds_tag="ACPC"

# Choose one of the following reference type which corresponds to different nifti volumes.
# AC or PC: LPS coordinates in the AC-PC-Midline coordinate system where AC or PC is (0,0,0)
ref_type='AC'


YELLOW='\033[0;33m'
GREEN='\033[0;32m'
NC='\033[0m'

cd $(dirname %0)

# Absolute data path
data_path_abs=$(readlink -f "$data_path")

if [ -d "$data_path/$subject/$session/dwi/mrtrix3" ]; then
    cd $data_path/$subject/$session/dwi/mrtrix3
else
    echo -e "${YELLOW}$subject/$session run mrtrix first.$NC"
    exit 1
fi

printf "\n${GREEN}Entering $subject/$session/dwi/mrtrix3/...$NC\n"

# Always use the reference volume and seed files from ses-00
if [ -d "../../../ses-1-00/anat" ]; then
    REF_nii=$(find ../../../ses-1-00/anat \( -name "${subject}_ses-00_$ref_type.nii" -o -name "${subject}_ses-1-00_$ref_type.nii.gz" \) | head -n 1)
fi

if [ -z "$REF_nii" ]; then
    echo -e "${YELLOW}${subject}_ses-1-00_$ref_type.nii(.gz) not found.$NC"
    exit 1
fi

if ! [ -f "$data_path_abs/target_seeds_${seeds_tag}.csv" ]; then
    echo -e "${YELLOW}$data_path_abs/target_seeds_${seeds_tag}.csv not found.$NC"
    exit 1
fi

# T1 volume in the dwi space
if ! [ -f "T1_FS_coreg.nii.gz" ]; then
    echo -e "${YELLOW}T1_FS_coreg.nii.gz not found.$NC"
    exit 1
fi

if ! [ -d "matched_targeting_tracks_${seeds_tag}" ]; then
    mkdir matched_targeting_tracks_${seeds_tag}
fi

cd matched_targeting_tracks_${seeds_tag}

# Register the ref volume to T1_FS_coreg.nii.gz
# ${ref_type}2dwi_0GenericAffine.mat and ${ref_type}2dwi_Warped.nii.gz will be produced
if ! [ -f "${ref_type}2dwi_0GenericAffine.mat" ]; then
    antsRegistrationSyNQuick.sh -d 3 -t r -f ../T1_FS_coreg.nii.gz -m "../$REF_nii" -o ${ref_type}2dwi_ -n $cores
    
    # Check registration
    mrview ../T1_FS_coreg.nii.gz -overlay.load AC2dwi_Warped.nii.gz -mode 2 &
fi 

# Use the track_file generated from BATMAN script or regenerate with matched parameters
# track_file="../tracks_10M.tck"

track_file="matched_tracks_$num_tracks.tck"
if ! [ -f "matched_tracks_$num_tracks.tck" ]; then
    tckgen -act ../5tt_coreg.mif -backtrack -seed_gmwmi ../gmwmSeed_coreg.mif -select $num_tracks ../wmfod_norm.mif matched_tracks_$num_tracks.tck \
           -maxlength 250 -cutoff 0.1 -step 0.5 -nthreads $cores -force
    
    # Clean up and exit if user presses Ctrl->C
    if ! [ $? -eq 0 ]; then
        rm matched_tracks_$num_tracks.tck
        exit 1
    fi
fi

# Apply the registration transformation to the seed points in the csv file.
# Note that ANTs software assume the coordinates in LPS system.
antsApplyTransformsToPoints -d 3 -i "$data_path_abs/target_seeds_${seeds_tag}.csv" -o "${subject}_seeds_to_dwi.csv" -t [${ref_type}2dwi_0GenericAffine.mat, 1]

exec 3< <(tail -n +2 "$data_path_abs/target_seeds_${seeds_tag}.csv")
exec 4< <(tail -n +2 "${subject}_seeds_to_dwi.csv")

# x0, y0, z0 are LPS coordinates in the reference volume
# x, y, z are transformed LPS coordinates in the dwi space
# r0 and r are radius of the seed sphere
while IFS=',' read -r x0 y0 z0 r0 label comment <&3 && IFS=',' read -r x y z r label comment <&4; do
    # Skip the nan lines (corresponding to blank lines in target_seeds_ACPC.csv) in the csv files
    if [ "$x" = "nan" ]; then
        continue
    fi
    
    # tckgen work with RAS coordinates. The following is conversion from LPS to RAS
    x=$(echo "$x * -1" | bc)
    y=$(echo "$y * -1" | bc)
    x0=$(echo "$x0 * -1" | bc)
    y0=$(echo "$y0 * -1" | bc)
    
    # Use rounded number in the output filenames
    x_rounded=$(echo $x | awk '{printf "%.1f", $1}')
    y_rounded=$(echo $y | awk '{printf "%.1f", $1}')
    z_rounded=$(echo $z | awk '{printf "%.1f", $1}')
    
    # Pierre suggested command line
    # tckgen -act {data_dir}/5tt_coreg_MNI.nii -backtrack -seed_dynamic {data_dir}/wmfod_norm_MNI.nii -nthreads 30 -cutoff 0.06 -maxlength 250 -step 0.5 -select 100M {data_dir}/wmfod_norm_MNI.nii {data_dir}/tracks_MNI_100M.tck -crop_at_gmwmi

    tckedit -include $x,$y,$z,$radius $track_file tracks_${ref_type}_${x0}_${y0}_${z0}_RAS_${x_rounded}_${y_rounded}_${z_rounded}.tck -nthreads $cores -force
done

# Close the file descriptors
exec 3<&-
exec 4<&-

 



           