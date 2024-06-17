#!/bin/bash

#---------------------  User's variables to be modified ---------------------

cores=18 

num_tracks=500k

# Absolute or relative path of the data folder to where the script located
data_path=~/Work/fusOUD/FUS/sub-224-FUS
mask_path=~/Nextcloud/Study/fusOUD/Treatment_Masks

# Set to 0 to disable quality control popup
QC=0

#---------------------------------------------------------------------------

cd $(dirname %0)

# Find all the subfolders named dwi and saved the paths to sessions_dir as an array
mapfile -t sessions_dir < <(find $data_path -type d -name dwi)

basedir=$(pwd)

YELLOW='\033[0;33m'
GREEN='\033[0;32m'
NC='\033[0m'
BOLD='\033[1m'

for (( n=0; n<${#sessions_dir[@]}; n++ ))
do
    printf "\n${YELLOW}Entering ${sessions_dir[$n]} ...$NC\n"
    cd ${sessions_dir[$n]}
    echo ${sessions_dir[$n]}
    
# The following use the session path to construct subject name, e.g. /FUS/sub-212/ses-1/dwi -> sub-212_ses-1    
    IFS='/' read -ra parts <<< ${sessions_dir[$n]}
    N=${#parts[@]}
    sub_name="${parts[N-3]}"
    ses_name="${parts[N-2]}"
    
    if ! [ -f "mrtrix3/FS2dwi_0GenericAffine.mat" ]; then
        echo "Run BATMAN script first!"
        cd $basedir
        continue
    fi
    
    if ! [ -f "$mask_path/${sub_name}_ses-00_mask_Left_T1.nii.gz" ]; then
        echo "No mask files found!"
        cd $basedir
        continue
    fi
    
    if ! [ -d connectome_from_mask ]; then
        mkdir connectome_from_mask
    fi
    cd connectome_from_mask
    
    if ! [ -f "mask_right_coreg.nii.gz" ]; then
        antsRegistrationSyNQuick.sh -d 3 -t r -f ../mrtrix4/T1_FS_coreg.nii.gz -m $mask_path/${sub_name}_ses-00_T1w.nii.gz -o T1w2dwi_ -n $cores
        
        # Use matlab method to apply image header transformation, avoiding interpolation of image data
        # Requires apply_rigid_transform.m in the script folder
        matlab -batch "addpath('$basedir'); apply_rigid_transform('$mask_path/${sub_name}_ses-00_mask_Left_T1.nii.gz', 'mask_left_coreg', 'T1w2dwi_0GenericAffine.mat')"
        matlab -batch "addpath('$basedir'); apply_rigid_transform('$mask_path/${sub_name}_ses-00_mask_Right_T1.nii.gz', 'mask_right_coreg', 'T1w2dwi_0GenericAffine.mat')"
        # mrconvert T1_FS_coreg.nii.gz T1_FS_coreg.mif -force
    fi

    # Check registration
    if ! [ $QC -eq 0 ]; then
        mrview ../mrtrix3/mean_b0_preprocessed.mif -overlay.load T1w2dwi_Warped.nii.gz -overlay.load mask_left_coreg.nii.gz -overlay.load mask_right_coreg.nii.gz &
    fi

# Tracks generation from masks
    if ! [ -f "tracks_from_mask_${num_tracks}.tck" ]; then
        tckgen -act ../mrtrix4/5tt_coreg.mif -backtrack -seed_image mask_left_coreg.nii.gz -seed_image mask_right_coreg.nii.gz \
               -select $num_tracks ../mrtrix3/wmfod_norm.mif tracks_from_mask_${num_tracks}.tck -nthreads $cores
        
        # Clean up and exit if user presses Ctrl->C
        if ! [ $? -eq 0 ]; then
            rm tracks_from_mask_${num_tracks}.tck
            exit 1
        fi
    fi
    
    # tckedit tracks_$num_tracks.tck -number 200k smallerTracks_200k.tck -force

    # mrview ${sub_name}_den_preproc_unbiased.mif -tractography.load smallerTracks_200k.tck
    # if ! [ -f "sift_from_mask_1M.tck" ]; then
    #    tcksift -act ../mrtrix4/5tt_coreg.mif -term_number 1M tracks_from_mask_${num_tracks}.tck ../mrtrix3/wmfod_norm.mif sift_from_mask_1M.tck -nthreads $cores -force
    # fi
    

    if ! [ -f "sift_from_mask_${num_tracks}.txt" ]; then
        tcksift2 -act ../mrtrix4/5tt_coreg.mif -out_mu sift_mu_left.txt -out_coeffs sift_coeffs_left.txt tracks_from_mask_${num_tracks}.tck ../mrtrix3/wmfod_norm.mif sift_from_mask_${num_tracks}.txt -nthreads $cores -force
        tckedit tracks_from_mask_${num_tracks}.tck -number 200k smallerTracts_200k.tck
    fi
            
    if ! [ $QC -eq 0 ]; then
        mrview ../mrtrix4/T1_FS_coreg.nii.gz –tractography.load smallerTracts_200k.tck &
    fi
    
    echo -e "${GREEN}${sessions_dir[$n]} ACT done.$NC"


# ----------------- Connectome from freesurfer parcels -----------------------
 
                     
    # Connectome with individual freesurfer atlas regions
    if ! [ -f "${sub_name}_${ses_name}_${num_tracks}_connectome.csv" ]; then
        matlab -batch "addpath('$basedir'); apply_rigid_transform('$mask_path/${sub_name}_ses-00_parcels_with_mask.nii.gz', 'parcels_with_mask_coreg', 'T1w2dwi_0GenericAffine.mat')"
        
        tck2connectome -symmetric -zero_diagonal -scale_invnodevol \
                       -tck_weights_in sift_from_mask_${num_tracks}.txt tracks_from_mask_${num_tracks}.tck parcels_with_mask_coreg.nii.gz \
                       ${sub_name}_${ses_name}_mask_connectome.csv -force \
                       -out_assignment ${sub_name}_${ses_name}_${num_tracks}_connectome_assignments.csv    
    fi
  
    echo -e "${YELLOW}${BOLD}All done for ${sessions_dir[$n]}.$NC  $(date '+%Y-%m-%d %H:%M')" 

    cd $basedir
done


