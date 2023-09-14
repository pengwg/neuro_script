function apply_rigid_transform(input_image, output_image, ants_transform)

% Apply the rigid transform matrix from ANTs to the header directional
% matrix of the input image

% input_image = 'FUS/sub-219-FUS/ses-00/anat/sub-219-FUS_ses-00_T1w.nii';
% output_image = 'T1_coreg_matlab';
% ants_transform = 'FUS/sub-219-FUS/ses-00/dwi/mrtrix/T1todwi_0GenericAffine.mat';

info = niftiinfo(input_image);
image_data = niftiread(info);
% Access the directional matrix from the input image's header
directional_matrix = info.Transform.T;

ants_transform_data = load(ants_transform);

% Define the rigid transformation matrix
local_transform = [reshape(ants_transform_data.AffineTransform_double_3_3, [3 4]);
                   0 0 0 1];
center_transform = diag([1 1 1 1]);
center_transform(1:3, 4) = ants_transform_data.fixed;
center_transform_reverse = diag([1 1 1 1]);
center_transform_reverse(1:3, 4) = -ants_transform_data.fixed;

rigid_transform = center_transform_reverse * local_transform * center_transform;

% Conversion from LPS to RAS system
lps2ras = diag([-1, -1, 1, -1]);
rigid_transform_ras = lps2ras * rigid_transform / lps2ras;

% Apply the rigid transformation by modifying the directional matrix
new_directional_matrix = directional_matrix * rigid_transform_ras';

% Update the directional matrix in the input image's header
info.Transform.T(1:12) = new_directional_matrix(1:12);

niftiwrite(image_data, output_image, info, 'Compressed', true);
disp(['Matlab: rigid transformation applied to ' input_image]);
% fprintf('\x1b[32m%s\x1b[0m\n', ['Rigid transformation applied to ' input_image]);

end