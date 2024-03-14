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

rigid_transform_ras = ea_antsmat2mat(ants_transform_data.AffineTransform_double_3_3, ...
                                     ants_transform_data.fixed);

% Apply the rigid transformation by modifying the directional matrix
new_directional_matrix = directional_matrix * rigid_transform_ras';

% Update the directional matrix in the input image's header
info.Transform.T(1:12) = new_directional_matrix(1:12);

niftiwrite(image_data, output_image, info, 'Compressed', true);
disp(['Matlab: rigid transformation applied to ' input_image]);
% fprintf('\x1b[32m%s\x1b[0m\n', ['Rigid transformation applied to ' input_image]);

end

% https://github.com/ANTsX/ANTs/wiki/ITK-affine-transform-conversion
% https://github.com/InsightSoftwareConsortium/ITK/blob/master/Modules/Core/Transform/include/itkMatrixOffsetTransformBase.hxx
% https://github.com/netstim/leaddbs/blob/master/helpers/ea_antsmat2mat.m

function mat = ea_antsmat2mat(afftransform, m_Center)
% followed instructions from
% https://www.neuro.polymtl.ca/tips_and_tricks/how_to_use_ants and ITK code
% in ComputeOffset() itkMatrixOffsetTransformBase.hxx

mat=[reshape(afftransform(1:9),[3,3])',afftransform(10:12)];

m_Translation=mat(:,4);
mat=[mat;[0,0,0,1]];

for i=1:3
    m_Offset(i) = m_Translation(i) + m_Center(i);
    for j=1:3
       m_Offset(i) = m_Offset(i)-(mat(i,j) * m_Center(j));  % (i,j) should remain the same since in C indexing rows before cols, too.
    end
end

mat(1:3,4)=m_Offset;
mat=inv(mat);

% convert RAS to LPS (ITK uses LPS)
mat=mat.*...
    [1  1 -1 -1
     1  1 -1 -1
    -1 -1  1  1
     1  1  1  1];
end