#!/usr/bin/python

import nibabel as nib
import numpy as np
from scipy.ndimage import affine_transform

def resample_to_target(source_path, target_path, output_path):
    # Load the source and target NIfTI images
    source_img = nib.load(source_path)
    target_img = nib.load(target_path)

    # Get the data and affine matrices
    source_data = source_img.get_fdata()
    target_data = target_img.get_fdata()
    source_affine = source_img.affine
    target_affine = target_img.affine

    # Compute the transformation matrix
    source_to_target_affine = np.linalg.inv(source_affine).dot(target_affine)

    # Calculate the new shape of the source data to match the target resolution
    target_shape = target_img.shape

    # Resample the source data to the target resolution
    resampled_data = affine_transform(
        source_data,
        source_to_target_affine[:3, :3],
        offset=source_to_target_affine[:3, 3],
        output_shape=target_shape,
        order=0 
    )

    # Convert resampled data to int32
    resampled_data = resampled_data.astype(np.uint32)
    
    # Save the resampled image
    resampled_img = nib.Nifti1Image(resampled_data, target_affine)
    nib.save(resampled_img, output_path)

def main():
    resample_to_target('mask.nii', 'fs_parcels.nii', 'resampled_mask.nii')

if __name__ == '__main__':
    main()

