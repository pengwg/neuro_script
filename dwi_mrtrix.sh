#!/bin/bash

cores=8

# Absolute or relative path of the data folder to where the script located
data_path=FUS/

# Freesurfer subject path
# SUBJECTS_DIR=/vulcan/Work/fusOUD/FUS/FS/

cd $(dirname %0)

# Find all the subfolders named dwi and saved the paths to sessions_dir as an array
mapfile -t sessions_dir < <(find $data_path -type d -name dwi)

basedir=$(pwd)

YELLOW='\033[0;33m'
GREEN='\033[0;32m'
NC='\033[0m'

for (( n=0; n<${#sessions_dir[@]}; n++ ))
do
    printf "\n${GREEN}Entering ${sessions_dir[$n]} ...$NC\n"
    cd ${sessions_dir[$n]}
    if ! [ -d mrtrix ]; then
        mkdir mrtrix
    fi
    
    chmod a+x *
# Search NIFTIs by HIGH_RES and 2mm_PA in filesnames and convert them to mif
    sub_dwi_nii=$(ls *HIGH_RES_??.nii* | head -n 1)
    if [ -z "$sub_dwi_nii" ]; then
        echo -e "${YELLOW}DWI files not found in ${sessions_dir[$n]}.$NC"
        cd $basedir
        continue
    fi
    sub_dwi=$(echo "$sub_dwi_nii" | sed 's/\.nii.*$//')
    
# The following use the session path to construct subject name, e.g. /FUS/sub-212/ses-1/dwi -> sub-212_ses-1    
    IFS='/' read -ra parts <<< ${sessions_dir[$n]}
    N=${#parts[@]}
    sub_name="${parts[N-3]}_${parts[N-2]}"
        
    if ! [ -f "mrtrix/$sub_name.mif" ]; then
        mrconvert $sub_dwi_nii mrtrix/$sub_name.mif -fslgrad $sub_dwi.bvec $sub_dwi.bval    
    fi
    
    sub_PA_niis=$(ls *2mm_PA*.nii*)
    for sub_PA_nii in $sub_PA_niis; do
        sub_PA=$(echo "$sub_PA_nii" | sed 's/\.nii.*$//')
        if ! [ -f "mrtrix/$sub_PA.mif" ]; then
            mrconvert $sub_PA_nii mrtrix/$sub_PA.mif
        fi
    done
    
    cd mrtrix
    chmod a+x *
    
# Denoise and degibbs
    if ! [ -f "${sub_name}_den.mif" ]; then
        dwidenoise $sub_name.mif ${sub_name}_den.mif -noise noise.mif -nthreads $cores
        if ! [ $? -eq 0 ]; then
            rm ${sub_name}_den.mif noise.mif
            exit 1
        fi
        mrcalc $sub_name.mif  ${sub_name}_den.mif -subtract residual.mif -force
    fi
    
    if ! [ -f "${sub_name}_den_unr.mif" ]; then
        mrdegibbs ${sub_name}_den.mif ${sub_name}_den_unr.mif -nthreads $cores
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

# Wrapper for FSL's topup and eddy
    if ! [ -f "${sub_name}_den_unr_preproc.mif" ]; then
        dwifslpreproc ${sub_name}_den_unr.mif ${sub_name}_den_unr_preproc.mif -pe_dir AP -rpe_pair -se_epi b0_pair.mif -eddy_options " --slm=linear --data_is_shelled"
        if ! [ $? -eq 0 ]; then
            rm ${sub_name}_den_unr_preproc.mif
            exit 1
        fi
    fi

# Bias correction with ANTs
    if ! [ -f "${sub_name}_den_unr_preproc_unbiased.mif" ]; then
        dwibiascorrect ants ${sub_name}_den_unr_preproc.mif ${sub_name}_den_unr_preproc_unbiased.mif -bias bias.mif -nthreads $cores
        if ! [ $? -eq 0 ]; then
            rm ${sub_name}_den_unr_preproc_unbiased.mif 
            exit 1
        fi
    fi
   
# Estimate response function(s) for spherical deconvolution
    if ! [ -f "wm.txt" ]; then
        dwi2response dhollander ${sub_name}_den_unr_preproc_unbiased.mif wm.txt gm.txt csf.txt -nthreads $cores
        if ! [ $? -eq 0 ]; then
            rm wm.txt gm.txt csf.txt 
            exit 1
        fi
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
        if ! [ $? -eq 0 ]; then
            rm wmfod.mif gmfod.mif csffod.mif 
            exit 1
        fi
    fi
    
# Create volume fraction for result inspection
    if ! [ -f "vf.mif" ]; then
        mrconvert -coord 3 0 wmfod.mif - | mrcat csffod.mif gmfod.mif - vf.mif    
    fi
    # mrview vf.mif -odf.load_sh wmfod.mif
 
# Intensity normalization 
    if ! [ -f "wmfod_norm.mif" ]; then
        mtnormalise wmfod.mif wmfod_norm.mif gmfod.mif gmfod_norm.mif csffod.mif csffod_norm.mif -mask mask.mif -nthreads $cores
        if ! [ $? -eq 0 ]; then
            rm wmfod_norm.mif gmfod_norm.mif csffod_norm.mif
            exit 1
        fi
    fi
    
    echo -e "${GREEN}${sessions_dir[$n]} FOD done.$NC"
    chmod a+x *
    
# Run DTIFIT
    if ! [ -f "${sub_name}_dti_V1.nii.gz" ]; then
        dtifit -k ${sub_name}_den_unr_preproc_unbiased.nii.gz \
               -o ${sub_name}_dti \
               -m ${sub_name}_den_unr_preproc_unbiased_mask.nii.gz \
               -r ../$sub_dwi.bvec \
               -b ../$sub_dwi.bval
    fi
    
    echo -e "${GREEN}${sessions_dir[$n]} DTI done.$NC"
    
# Find and convert anatomical T1 to mif
    if [ -d "../../anat" ]; then
        sub_T1_nii=$(find ../../anat -name *T1*.nii* | head -n 1)
    fi
    
    if [ -z "$sub_T1_nii" ]; then
        echo -e "${YELLOW}T1 images not found in ${sessions_dir[$n]}/../../anat..$NC"
        cd $basedir
        continue
    fi

# Create 5tt registered T1 volume and gray matter/white matter boundary seed
    if ! [ -f "5tt_nocoreg.mif" ]; then
        mrconvert $sub_T1_nii T1_raw.mif -force
        chmod a+x T1_raw.mif
        5ttgen fsl T1_raw.mif 5tt_nocoreg.mif -nthreads $cores
    fi
    
    if ! [ -f "gmwmSeed_coreg.mif" ]; then
        dwiextract ${sub_name}_den_unr_preproc_unbiased.mif - -bzero | mrmath - mean mean_b0_preprocessed.mif -axis 3 -force
        mrconvert mean_b0_preprocessed.mif mean_b0_preprocessed.nii.gz -force
        flirt -in mean_b0_preprocessed.nii.gz -ref "$sub_T1_nii" -dof 6 -omat diff2struct_fsl.mat -verbose 1
    
        transformconvert diff2struct_fsl.mat mean_b0_preprocessed.nii.gz T1_raw.mif flirt_import diff2struct_mrtrix.txt -force
        mrtransform T1_raw.mif -linear diff2struct_mrtrix.txt -inverse T1_coreg.mif -force
        mrtransform 5tt_nocoreg.mif -linear diff2struct_mrtrix.txt -inverse 5tt_coreg.mif -force
    
        5tt2gmwmi 5tt_coreg.mif gmwmSeed_coreg.mif -nthreads $cores
    fi

    if ! [ -f "tracks_10M.tck" ]; then
        tckgen -act 5tt_coreg.mif -backtrack -seed_gmwmi gmwmSeed_coreg.mif -select 10000k wmfod_norm.mif tracks_10M.tck -nthreads $cores
        if ! [ $? -eq 0 ]; then
            rm tracks_10M.tck
            exit 1
        fi
    fi
    # tckedit tracks_10M.tck -number 200k smallerTracks_200k.tck -force

    # mrview ${sub_name}_den_preproc_unbiased.mif -tractography.load smallerTracks_200k.tck
    if ! [ -f "sift_1M.tck" ]; then
        tcksift -act 5tt_coreg.mif -term_number 1000k tracks_10M.tck wmfod_norm.mif sift_1M.tck -nthreads $cores
    fi

    if ! [ -f "sift_1M.txt" ]; then
        tcksift2 -act 5tt_coreg.mif -out_mu sift_mu.txt -out_coeffs sift_coeffs.txt sift_1M.tck wmfod_norm.mif sift_1M.txt -nthreads $cores
    fi
    
    echo -e "${GREEN}${sessions_dir[$n]} ACT done.$NC"
    
    fs_subject="FS_$sub_name"
    
    if ! [ -f "$SUBJECTS_DIR/$fs_subject/mri/aparc+aseg.mgz" ]; then
        echo -e "${YELLOW}$SUBJECTS_DIR/$fs_subject/mri/aparc+aseg.mgz not found.$NC"
        cd $basedir
        continue
    fi
 
 # Connectome with individual freesurfer atlas to get the regions
    if ! [ -f "fs_parcels.csv" ]; then      
        labelconvert $SUBJECTS_DIR/$fs_subject/mri/aparc+aseg.mgz \
                     $FREESURFER_HOME/FreeSurferColorLUT.txt \
                     $(dirname $(which mrview))/../share/mrtrix3/labelconvert/fs_default.txt \
                     fs_parcels.mif -force
           
        tck2connectome -symmetric -zero_diagonal -scale_invnodevol \
                       -tck_weights_in sift_1M.txt tracks_10M.tck fs_parcels.mif fs_parcels.csv \
                       -out_assignment fs_assignments_parcels.csv
    fi
    
    echo -e "${GREEN}${sessions_dir[$n]} connectome done.$NC"
    cd $basedir
done


