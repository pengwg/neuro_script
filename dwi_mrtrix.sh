#!/bin/bash

cores=6

# Find all the subfolders named dwi and saved the paths to sessions_dir as an array
mapfile -t sessions_dir < <(find FUS/ -type d -name dwi)

basedir=$(dirname $0)

for (( n=0; n<${#sessions_dir[@]}; n++ ))
do
    cd ${sessions_dir[$n]}
    if ! [ -d mrtrick ]; then
        mkdir mrtrick
    fi

# Search NIFTIs by HIGH_RES and 2mm_PA in filesnames and convert them to mif
    sub_dwi_nii=$(find . -name *HIGH_RES*.nii* | head -n 1)
    if [ -z "$sub_dwi_nii" ]; then
        echo "No dwi files in ${sessions_dir[$n]}."
        continue
    fi
    
    sub_dwi=$(echo "$sub_dwi_nii" | sed 's/\.nii.*$//')
    if ! [ -f "mrtrick/$sub_dwi.mif" ]; then
        mrconvert $sub_dwi_nii mrtrick/$sub_dwi.mif -fslgrad $sub_dwi.bvec $sub_dwi.bval    
    fi
    
    sub_PA_niis=$(ls *2mm_PA*.nii*)
    for sub_PA_nii in $sub_PA_niis; do
        sub_PA=$(echo "$sub_PA_nii" | sed 's/\.nii.*$//')
        if ! [ -f "mrtrick/$sub_PA.mif" ]; then
            mrconvert $sub_PA_nii mrtrick/$sub_PA.mif
        fi
    done
    
    cd mrtrick
    
# Denoise and degibbs
    if ! [ -f "${sub_dwi}_den.mif" ]; then
        dwidenoise $sub_dwi.mif ${sub_dwi}_den.mif -noise noise.mif -nthreads $cores
        if ! [ $? -eq 0 ]; then
            rm ${sub_dwi}_den.mif noise.mif
            exit 1
        fi
        mrcalc $sub_dwi.mif  ${sub_dwi}_den.mif -subtract residual.mif -force
    fi
    
    if ! [ -f "${sub_dwi}_den_unr.mif" ]; then
        mrdegibbs ${sub_dwi}_den.mif ${sub_dwi}_den_unr.mif -nthreads $cores
        if ! [ $? -eq 0 ]; then
            rm ${sub_dwi}_den_unr.mif
            exit 1
        fi
    fi

# Compute b0 AP and PA
    if ! [ -f b0_pair.mif ]; then
        dwiextract ${sub_dwi}_den_unr.mif - -bzero | mrmath - mean mean_b0_AP.mif -axis 3 -force
        mrcat *2mm_PA*.mif -axis 3 - | mrmath - mean mean_b0_PA.mif -axis 3 -force
        mrcat mean_b0_AP.mif mean_b0_PA.mif -axis 3 b0_pair.mif
    fi

# Wrapper for FSL's topup and eddy
    if ! [ -f "${sub_dwi}_den_unr_preproc.mif" ]; then
        dwifslpreproc ${sub_dwi}_den_unr.mif ${sub_dwi}_den_unr_preproc.mif -pe_dir AP -rpe_pair -se_epi b0_pair.mif -topup_options " --nthr="$cores -eddy_options " --slm=linear --data_is_shelled"
        if ! [ $? -eq 0 ]; then
            rm ${sub_dwi}_den_unr_preproc.mif
            exit 1
        fi
    fi

# Bias correction with ANTs
    if ! [ -f "${sub_dwi}_den_unr_preproc_unbiased.mif" ]; then
        dwibiascorrect ants ${sub_dwi}_den_unr_preproc.mif ${sub_dwi}_den_unr_preproc_unbiased.mif -bias bias.mif -nthreads $cores
        if ! [ $? -eq 0 ]; then
            rm ${sub_dwi}_den_unr_preproc_unbiased.mif 
            exit 1
        fi
    fi
   
# Estimate response function(s) for spherical deconvolution
    if ! [ -f "wm.txt" ]; then
        dwi2response dhollander ${sub_dwi}_den_unr_preproc_unbiased.mif wm.txt gm.txt csf.txt -nthreads $cores
        if ! [ $? -eq 0 ]; then
            rm wm.txt gm.txt csf.txt 
            exit 1
        fi
    fi
    
# Generate fibre orientation distributions
    if ! [ -f "mask.mif" ]; then
        dwi2mask ${sub_dwi}_den_unr_preproc_unbiased.mif mask.mif -force -nthreads $cores
        if ! [ $? -eq 0 ]; then
            rm mask.mif
            exit 1
        fi
    fi
    
    if ! [ -f "wmfod.mif" ]; then
        dwi2fod msmt_csd ${sub_dwi}_den_unr_preproc_unbiased.mif -mask mask.mif wm.txt wmfod.mif gm.txt gmfod.mif csf.txt csffod.mif -nthreads $cores
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

# Find and convert anatomical T1 to mif
    if ! [ -d "../../anat" ]; then
        sub_T1_nii=$(find ../../anat -name *T1*.nii* | head -n 1)
    fi
    
    if [ -z "$sub_T1_nii" ]; then
        echo "No T1 images found in ${sessions_dir[$n]}../../anat."
        continue
    fi
    
    if ! [ -f "5tt_nocoreg.mif" ]; then
        mrconvert $sub_T1_nii T1_raw.mif -force
        5ttgen fsl T1_raw.mif 5tt_nocoreg.mif
    fi
    
    dwiextract ${sub_dwi}_den_unr_preproc_unbiased.mif - -bzero | mrmath – mean mean_b0_preprocessed.mif –axis 3 -force
    mrconvert mean_b0_preprocessed.mif mean_b0_preprocessed.nii.gz -force
    flirt –in mean_b0_preprocessed.nii.gz –ref $sub_T1_nii –dof 6 –omat diff2struct_fsl.mat
    
    transformconvert diff2struct_fsl.mat mean_b0_preprocessed.nii.gz T1_raw.mif flirt_import diff2struct_mrtrix.txt -force
    mrtransform T1_raw.mif –linear diff2struct_mrtrix.txt –inverse T1_coreg.mif -force
    mrtransform 5tt_nocoreg.mif –linear diff2struct_mrtrix.txt –inverse 5tt_coreg.mif
    
    5tt2gmwmi 5tt_coreg.mif gmwmSeed_coreg.mif
    tckgen –act 5tt_coreg.mif –backtrack –seed_gmwmi gmwmSeed_coreg.mif –select 10000000 wmfod_norm.mif tracks_10mio.tck
    tckedit tracks_10mio.tck –number 200k smallerTracks_200k.tck
    tcksift –act 5tt_coreg.mif –term_number 1000000 tracks_10mio.tck wmfod_norm.mif sift_1mio.tck
    
    echo "${sessions_dir[$n]} processing finished."
    cd $basedir
done


