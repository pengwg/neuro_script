#!/bin/bash

cores=6

# Find all the subfolders named dwi and saved the paths to sessions_dir as an array
mapfile -t sessions_dir < <(find FUS/sub-210-FUS/sub-210-FUS_Day7/ -type d -name dwi)

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
        mrcalc $sub_name.mif  ${sub_name}_den.mif -subtract residual.mif -force
    fi
    
    if ! [ -f "${sub_name}_den_unr.mif" ]; then
        mrdegibbs ${sub_name}_den.mif ${sub_name}_den_unr.mif -nthreads $cores
    fi

# Compute b0 AP and PA
    if ! [ -f b0_pair.mif ]; then
        dwiextract ${sub_name}_den_unr.mif - -bzero | mrmath - mean mean_b0_AP.mif -axis 3 -force
        mrcat *2mm_PA*.mif -axis 3 - | mrmath - mean mean_b0_PA.mif -axis 3 -force
        mrcat mean_b0_AP.mif mean_b0_PA.mif -axis 3 b0_pair.mif
    fi

# Wrapper for FSL's topup and eddy
    if ! [ -f "${sub_name}_preproc.mif" ]; then
        dwifslpreproc ${sub_name}_den_unr.mif ${sub_name}_preproc.mif -pe_dir AP -rpe_pair -se_epi b0_pair.mif -topup_options " --nthr="$cores -eddy_options " --slm=linear --data_is_shelled"
    fi

# Bias correction with ANTs
    if ! [ -f "${sub_name}_preproc_unbiased.mif" ]; then
        dwibiascorrect ants ${sub_name}_preproc.mif ${sub_name}_preproc_unbiased.mif -bias bias.mif -nthreads $cores
    fi
   
# Estimate response function(s) for spherical deconvolution
    if ! [ -f "wm.txt" ]; then
        dwi2response dhollander ${sub_name}_preproc_unbiased.mif wm.txt gm.txt csf.txt -nthreads $cores
    fi
    
# Generate fibre orientation distributions
    if ! [ -f "wmfod.mif" ]; then
        dwi2mask ${sub_name}_preproc_unbiased.mif mask.mif -force -nthreads $cores
        dwi2fod msmt_csd ${sub_name}_preproc_unbiased.mif -mask mask.mif wm.txt wmfod.mif gm.txt gmfod.mif csf.txt csffod.mif -nthreads $cores
    fi
    
    # mrconvert -coord 3 0 wmfod.mif - | mrcat csffod.mif gmfod.mif - vf.mif    
    # mtnormalise wmfod.mif wmfod_norm.mif gmfod.mif gmfod_norm.mif csffod.mif csffod_norm.mif -mask mask.mif
    
    echo "$sub_name processed"
    cd $basedir
done


