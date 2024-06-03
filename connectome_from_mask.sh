#!/bin/bash

#---------------------  User's variables to be modified ---------------------

cores=18 

# Absolute or relative path of the data folder to where the script located
data_path=~/Work/fusOUD/FUS/sub-224-FUS

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
    sub_name="${parts[N-3]}_${parts[N-2]}"
    
    if ! [ -f "mrtrix3/FS2dwi_0GenericAffine.mat" ]; then
        echo "Run BATMAN script first!"
        cd $basedir
        continue
    fi
    
    if ! [ -d connectome_mask ]; then
        mkdir connectome_mask
    fi
    cd connectome_mask
    
    if ! [ -f "aparc+aseg_coreg.nii.gz" ]; then       
        # Use matlab method to apply image header transformation, avoiding interpolation of image data
        # Requires apply_rigid_transform.m in the script folder
        matlab -batch "addpath('$basedir'); apply_rigid_transform('../mrtrix3/aparc+aseg.nii.gz', 'aparc+aseg_coreg', '../mrtrix3/FS2dwi_0GenericAffine.mat')"
        # Also apply the transform to anatomic volume for checking purpose
        matlab -batch "addpath('$basedir'); apply_rigid_transform('../mrtrix3/T1_FS.nii.gz', 'T1_FS_coreg', '../mrtrix3/FS2dwi_0GenericAffine.mat')"
        mrconvert T1_FS_coreg.nii.gz T1_FS_coreg.mif -force
    fi

    # Check registration
    if ! [ $QC -eq 0 ]; then
        mrview ../mrtrix3/mean_b0_preprocessed.mif -overlay.load T1_FS_coreg.nii.gz -overlay.load aparc+aseg_coreg.nii.gz &
    fi

# Create 5tt registered T1 volume and gray matter/white matter boundary seed using freesurfer segmentation
    if ! [ -f "5tt_coreg.mif" ]; then
        # '5ttgen freesurfer' will invoke labelconvert by itself  
        5ttgen freesurfer aparc+aseg_coreg.nii.gz 5tt_coreg.mif -nthreads $cores
        
        # Continue script if 5ttgen failed
        if ! [ $? -eq 0 ]; then
            cd $basedir
            continue
        fi
    fi
    
    if ! [ $QC -eq 0 ]; then
        mrview T1_FS_coreg.nii.gz -overlay.load 5tt_coreg.mif &
    fi
    
    if ! [ -f "gmwmSeed_coreg.mif" ]; then
        5tt2gmwmi 5tt_coreg.mif gmwmSeed_coreg.mif -nthreads $cores -force
    fi
    
    chmod a+x *
    
    if ! [ -f "tracks_10M.tck" ]; then
        tckgen -act 5tt_coreg.mif -backtrack -seed_gmwmi gmwmSeed_coreg.mif -select 10000k ../mrtrix3/wmfod_norm.mif tracks_10M.tck -nthreads $cores
        
        # Clean up and exit if user presses Ctrl->C
        if ! [ $? -eq 0 ]; then
            rm tracks_10M.tck
            exit 1
        fi
    fi
    # tckedit tracks_10M.tck -number 200k smallerTracks_200k.tck -force

    # mrview ${sub_name}_den_preproc_unbiased.mif -tractography.load smallerTracks_200k.tck
    if ! [ -f "sift_1M.tck" ]; then
        tcksift -act 5tt_coreg.mif -term_number 1M tracks_10M.tck ../mrtrix3/wmfod_norm.mif sift_1M.tck -nthreads $cores -force
    fi

    if ! [ -f "sift_10M.txt" ]; then
        tcksift2 -act 5tt_coreg.mif -out_mu sift_mu.txt -out_coeffs sift_coeffs.txt tracks_10M.tck ../mrtrix3/wmfod_norm.mif sift_10M.txt -nthreads $cores -force
        tckedit tracks_10M.tck -number 200k smallerTracts_200k.tck
    fi
            
    if ! [ $QC -eq 0 ]; then
        mrview T1_FS_coreg.nii.gz –tractography.load smallerTracts_200k.tck &
    fi
    
    echo -e "${GREEN}${sessions_dir[$n]} ACT done.$NC"
    note="${sessions_dir[$n]} ACT done.$NC  $(date '+%Y-%m-%d %H:%M')"
    echo -e "$note" >> $basedir/script4BATMAN_log.txt


# ----------------- Connectome from freesurfer parcels -----------------------
 
                     
    # Connectome with individual freesurfer atlas regions
    if ! [ -f "${sub_name}_10M_connectome.csv" ]; then
        labelconvert aparc+aseg_coreg.nii.gz \
                     $FREESURFER_HOME/FreeSurferColorLUT.txt \
                     $(dirname $(which mrview))/../share/mrtrix3/labelconvert/fs_default.txt \
                     fs_parcels_coreg.nii.gz -force
                     
        tck2connectome -symmetric -zero_diagonal -scale_invnodevol \
                       -tck_weights_in sift_10M.txt tracks_10M.tck fs_parcels_coreg.nii.gz \
                       ${sub_name}_10M_connectome.csv \
                       -out_assignment ${sub_name}_10M_connectome_assignments.csv
    fi
    
    # Connectome scaled with mean FA
    if ! [ -f "${sub_name}_meanFA_10M_connectome.csv" ]; then 
        # Computing fractional anisotropy of full 10M track file
        dwi2tensor ${sub_name}_den_unr_preproc_unbiased.mif tensor.mif -force -nthreads $cores
        tensor2metric tensor.mif -fa FA.mif -force -nthreads $cores  

        # Computing the mean FA of tracks 
        tcksample tracks_10M.tck FA.mif tracks_meanFA_10M.csv -stat_tck mean -force -nthreads $cores 
   
        tck2connectome -symmetric -zero_diagonal \
                       -tck_weights_in sift_10M.txt tracks_10M.tck fs_parcels_coreg.nii.gz \
                       ${sub_name}_meanFA_10M_connectome.csv \
                       -scale_file tracks_meanFA_10M.csv -stat_edge mean -force 
    fi
    
    chmod a+x *
    
    note="${sessions_dir[$n]} connectome done.$NC  $(date '+%Y-%m-%d %H:%M')" 
    echo -e "$note" >> $basedir/script4BATMAN_log.txt
   
    echo -e "${YELLOW}${BOLD}All done for ${sessions_dir[$n]}.$NC  $(date '+%Y-%m-%d %H:%M')" 

    note="All done for ${sessions_dir[$n]}.$NC  $(date '+%Y-%m-%d %H:%M')"
    echo -e "$note" >> $basedir/script4BATMAN_log.txt

    cd $basedir
done


