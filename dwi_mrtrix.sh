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
    sub_name_nii=$(ls *HIGH_RES*.nii* | head -n 1)
    sub_name=$(echo "$sub_name_nii" | sed 's/\.nii.*$//')
    if ! [ -f "mrtrick/$sub_name.mif" ]; then
        mrconvert $sub_name_nii mrtrick/$sub_name.mif -fslgrad $sub_name.bvec $sub_name.bval    
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
        dwifslpreproc ${sub_name}_den_unr.mif ${sub_name}_den_unr_preproc.mif -pe_dir AP -rpe_pair -se_epi b0_pair.mif -topup_options " --nthr="$cores -eddy_options " --slm=linear --data_is_shelled"
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
    if ! [ -f "mask.mif" ]; then
        dwi2mask ${sub_name}_den_unr_preproc_unbiased.mif mask.mif -force -nthreads $cores
        if ! [ $? -eq 0 ]; then
            rm mask.mif
            exit 1
        fi
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
    
    echo "${sub_name} processing finished."
    cd $basedir
done


