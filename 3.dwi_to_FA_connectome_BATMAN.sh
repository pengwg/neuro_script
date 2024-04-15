#!/bin/bash


# Created by Drs Wang and Thompson-Lake 
# June 2023 WVU RNI
# Up to date version at https://github.com/pengwg/neuro_script

# Requires files in BIDS naming format
# Softwares required: mrtrix3, FSL, ANTs, Matlab
# Requires segmentation and fs_default.txt from freesurfer

# functions within the script: 
# Combines the PA files into single file
# Performs denoising and degibbs correction on the DWI data.
# Performs preprocessing steps including b0 extraction, topup, eddy, and bias correction using ANTs.
# Estimates the response function for spherical deconvolution.
# Generates fiber orientation distributions (FODs).
# Performs intensity normalization and saves the normalized FODs.
# Runs DTIFIT to generate diffusion tensor images (DTI).
# Finds and converts the anatomical T1 image to the MIF format.
# Registers T1 image to DWI space and generates 5TT (tissue type) segmentation and seed masks.
# Performs tractography using tckgen and sifts the tracks using tcksift.
# Converts the FreeSurfer parcellation to format and computes the connectome using tck2connectome.

# BATMAN4 edits: commented out SIFT1, Freesurfer folder location changed, added mrview popups for QC, mrconvert T1_FS_correg to a mif file for viewing

#---------------------  User's variables to be modified ---------------------

cores=18

# Define external drive mount point you need these nest lines in all your scripts!!
external_drive="/media/dgt00003/dgytl"

# Define folder path and list of sessions
#TO EDIT DEPENDING ON FOLDER
relative_folder_path="CPO"

# Construct the complete folder path
data_path="$external_drive/$relative_folder_path"

# Absolute or relative path of the data folder to where the script located
# data_path=FUS/

# Set to 0 to disable quality control popup
QC=0

# Freesurfer subject path
SUBJECTS_DIR=/media/dgt00003/dgytl/FS


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

    if ! [ -d mrtrix3 ]; then
        mkdir mrtrix3
    fi
    
    chmod a+x *
# Search NIFTIs by HIGH_RES and 2mm_PA in filesnames and convert them to mif
    sub_dwi_nii=$(ls *HIGH_RES*.nii* | head -n 1)
    if [ -z "$sub_dwi_nii" ]; then
        echo -e "${YELLOW}DWI files not found in ${sessions_dir[$n]}.$NC" | tee -a "$data_path/script4BATMAN_log.txt"
        cd $basedir
        continue
    fi
    sub_dwi=$(echo "$sub_dwi_nii" | sed 's/\.nii.*$//')
    
# The following use the session path to construct subject name, e.g. /FUS/sub-212/ses-1/dwi -> sub-212_ses-1    
    IFS='/' read -ra parts <<< ${sessions_dir[$n]}
    N=${#parts[@]}
    sub_name="${parts[N-3]}_${parts[N-2]}"
    sub_name_ses_00="${parts[N-3]}_ses-00"
        
    if ! [ -f "mrtrix3/$sub_name.mif" ]; then
        mrconvert $sub_dwi_nii mrtrix3/$sub_name.mif -fslgrad $sub_dwi.bvec $sub_dwi.bval    
    fi
    
    sub_PA_niis=$(ls *2mm_PA*.nii*)
    for sub_PA_nii in $sub_PA_niis; do
        sub_PA=$(echo "$sub_PA_nii" | sed 's/\.nii.*$//')
        if ! [ -f "mrtrix3/$sub_PA.mif" ]; then
            mrconvert $sub_PA_nii mrtrix3/$sub_PA.mif
        fi
    done
    
    cd mrtrix3
    chmod a+x *
    
# Denoise and degibbs
    if ! [ -f "${sub_name}_den.mif" ]; then
        dwidenoise $sub_name.mif ${sub_name}_den.mif -noise noise.mif -nthreads $cores -force
        
        # Clean up and exit if user presses Ctrl->C
        if ! [ $? -eq 0 ]; then
            rm ${sub_name}_den.mif noise.mif
            exit 1
        fi
        mrcalc $sub_name.mif  ${sub_name}_den.mif -subtract residual.mif -force
    fi
    
    chmod a+x *
    
    if ! [ -f "${sub_name}_den_unr.mif" ]; then
        mrdegibbs ${sub_name}_den.mif ${sub_name}_den_unr.mif -nthreads $cores -force
        
        # Clean up and exit if user presses Ctrl->C
        if ! [ $? -eq 0 ]; then
            rm ${sub_name}_den_unr.mif
            exit 1
        fi
    fi

# Compute b0 AP and PA
    if ! [ -f b0_pair.mif ]; then
        dwiextract ${sub_name}_den_unr.mif - -bzero | mrmath - mean mean_b0_AP.mif -axis 3 -force
        mrcat *2mm_PA*.mif -axis 3 - | mrmath - mean mean_b0_PA.mif -axis 3 -force
        mrcat mean_b0_AP.mif mean_b0_PA.mif -axis 3 b0_pair.mif
    fi    
    
    # Check 3 PAs alignment
    if ! [ $QC -eq 0 ]; then
        mapfile -t PA_files < <(find . -type f -name \*2mm_PA\*)
        mrview ${PA_files[0]} -overlay.load ${PA_files[1]} -overlay.load ${PA_files[2]} &
    fi       

# Wrapper for FSL's topup and eddy
    if ! [ -f "${sub_name}_den_unr_preproc.mif" ]; then
        dwifslpreproc ${sub_name}_den_unr.mif ${sub_name}_den_unr_preproc.mif -pe_dir AP -rpe_pair -se_epi b0_pair.mif -topup_options " --nthr="$cores -eddy_options " --slm=linear --data_is_shelled"
        
        # Clean up and exit if user presses Ctrl->C
        if ! [ $? -eq 0 ]; then
            rm ${sub_name}_den_unr_preproc.mif
            exit 1
        fi
    fi

# Bias correction with ANTs
    if ! [ -f "${sub_name}_den_unr_preproc_unbiased.mif" ]; then
        dwibiascorrect ants ${sub_name}_den_unr_preproc.mif ${sub_name}_den_unr_preproc_unbiased.mif -bias bias.mif -nthreads $cores -force
        if ! [ $? -eq 0 ]; then
            rm ${sub_name}_den_unr_preproc_unbiased.mif 
            exit 1
        fi
    fi
    
    if ! [ $QC -eq 0 ]; then
        mrview ${sub_name}_den_unr_preproc_unbiased.mif -overlay.load bias.mif &
    fi
       
# Estimate response function(s) for spherical deconvolution
    if ! [ -f "wm.txt" ]; then
        dwi2response dhollander ${sub_name}_den_unr_preproc_unbiased.mif wm.txt gm.txt csf.txt -nthreads $cores
        
        # Continue script if dwi2response failed
        if ! [ $? -eq 0 ]; then
            rm wm.txt gm.txt csf.txt 
            cd $basedir
            continue
        fi
    fi    
    
    if ! [ $QC -eq 0 ]; then
        shview wm.txt &
    fi
    
# Generate fibre orientation distributions
    if ! [ -f "${sub_name}_den_unr_preproc_unbiased_mask.nii.gz" ]; then
        #dwi2mask ${sub_name}_den_unr_preproc_unbiased.mif mask.mif -force -nthreads $cores
        
        mrconvert ${sub_name}_den_unr_preproc_unbiased.mif ${sub_name}_den_unr_preproc_unbiased.nii.gz -force
        bet ${sub_name}_den_unr_preproc_unbiased.nii.gz ${sub_name}_den_unr_preproc_unbiased -m -n -f 0.2  
        mrconvert ${sub_name}_den_unr_preproc_unbiased_mask.nii.gz mask.mif
    fi
    
    if ! [ -f "wmfod.mif" ]; then
        dwi2fod msmt_csd ${sub_name}_den_unr_preproc_unbiased.mif -mask mask.mif wm.txt wmfod.mif gm.txt gmfod.mif csf.txt csffod.mif -nthreads $cores
        
        # Clean up and exit if user presses Ctrl->C
        if ! [ $? -eq 0 ]; then
            rm wmfod.mif gmfod.mif csffod.mif 
            exit 1
        fi
    fi
    
# Create volume fraction for result inspection
    if ! [ -f "vf.mif" ]; then
        mrconvert -coord 3 0 wmfod.mif - | mrcat csffod.mif gmfod.mif - vf.mif    
    fi      
    
    if ! [ $QC -eq 0 ]; then
        mrview ${sub_name}_den_unr_preproc_unbiased.mif -overlay.load vf.mif &
        # display the white matter FOD on a map which shows the estimated volume fraction of each tissue
        mrview vf.mif -odf.load_sh wmfod.mif &
    fi
 
# Intensity normalization 
    if ! [ -f "wmfod_norm.mif" ]; then
        mtnormalise wmfod.mif wmfod_norm.mif gmfod.mif gmfod_norm.mif csffod.mif csffod_norm.mif -mask mask.mif -nthreads $cores
        
        # Clean up and exit if user presses Ctrl->C
        if ! [ $? -eq 0 ]; then
            rm wmfod_norm.mif gmfod_norm.mif csffod_norm.mif
            exit 1
        fi
    fi
    
    echo -e "${GREEN}${sessions_dir[$n]} FOD done.$NC"
    
    note="${sessions_dir[$n]} FOD done.$NC  $(date '+%Y-%m-%d %H:%M')"
    echo -e "$note" >> $data_path/script4BATMAN_log.txt

    chmod a+x *
    

# Run DTIFIT with weighted least squares
    if ! [ -f "${sub_name}_dti_V1.nii.gz" ]; then
        dtifit -k ${sub_name}_den_unr_preproc_unbiased.nii.gz \
               -o ${sub_name}_dti \
               -m ${sub_name}_den_unr_preproc_unbiased_mask.nii.gz \
               -r ../$sub_dwi.bvec \
               -b ../$sub_dwi.bval \
               -w
    fi
    
    echo -e "${GREEN}${sessions_dir[$n]} DTI done.$NC"
    note="${sessions_dir[$n]} DTI done.$NC $(date '+%Y-%m-%d %H:%M')"
    echo -e "$note" >> "$data_path/script4BATMAN_log.txt"
   


# ----------------- Anatomically Constrained Tractography ---------------------



    fs_subject="FS_$sub_name"
    
    if ! [ -f "$SUBJECTS_DIR/$fs_subject/mri/aparc+aseg.mgz" ]; then
        echo -e "${YELLOW}$SUBJECTS_DIR/$fs_subject/mri/aparc+aseg.mgz not found.$NC ******" >> $data_path/script4BATMAN_log.txt
        cd $basedir
        continue
    fi
    
# Register freesurfer segmentation to dwi mean_b0
    if ! [ -f "aparc+aseg_coreg.nii.gz" ]; then
        mri_convert $SUBJECTS_DIR/$fs_subject/mri/T1.mgz T1_FS.nii.gz
        mri_convert $SUBJECTS_DIR/$fs_subject/mri/aparc+aseg.mgz aparc+aseg.nii.gz
        
        dwiextract ${sub_name}_den_unr_preproc_unbiased.mif - -bzero | mrmath - mean mean_b0_preprocessed.mif -axis 3 -force
        mrconvert mean_b0_preprocessed.mif mean_b0_preprocessed.nii.gz -force
        
        antsRegistrationSyNQuick.sh -d 3 -t r -f mean_b0_preprocessed.nii.gz -m T1_FS.nii.gz -o FS2dwi_ -n $cores
        
        # Use matlab method to apply image header transformation, avoiding interpolation of image data
        # Requires apply_rigid_transform.m in the script folder
        matlab -batch "addpath('$basedir'); apply_rigid_transform('aparc+aseg.nii.gz', 'aparc+aseg_coreg', 'FS2dwi_0GenericAffine.mat')"
        # Also apply the transform to anatomic volume for checking purpose
        matlab -batch "addpath('$basedir'); apply_rigid_transform('T1_FS.nii.gz', 'T1_FS_coreg', 'FS2dwi_0GenericAffine.mat')"
        mrconvert T1_FS_coreg.nii.gz T1_FS_coreg.mif -force
    fi

    # Check registration
    if ! [ $QC -eq 0 ]; then
        mrview mean_b0_preprocessed.mif -overlay.load T1_FS_coreg.nii.gz -overlay.load aparc+aseg_coreg.nii.gz &
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
        tckgen -act 5tt_coreg.mif -backtrack -seed_gmwmi gmwmSeed_coreg.mif -select 10000k wmfod_norm.mif tracks_10M.tck -nthreads $cores
        
        # Clean up and exit if user presses Ctrl->C
        if ! [ $? -eq 0 ]; then
            rm tracks_10M.tck
            exit 1
        fi
    fi
    # tckedit tracks_10M.tck -number 200k smallerTracks_200k.tck -force

    # mrview ${sub_name}_den_preproc_unbiased.mif -tractography.load smallerTracks_200k.tck
    if ! [ -f "sift_1M.tck" ]; then
        tcksift -act 5tt_coreg.mif -term_number 1M tracks_10M.tck wmfod_norm.mif sift_1M.tck -nthreads $cores -force
    fi

    if ! [ -f "sift_10M.txt" ]; then
        tcksift2 -act 5tt_coreg.mif -out_mu sift_mu.txt -out_coeffs sift_coeffs.txt tracks_10M.tck wmfod_norm.mif sift_10M.txt -nthreads $cores -force
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


