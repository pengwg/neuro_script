clear

mask_path = '~/Nextcloud/Study/fusOUD/Treatment_Masks';
fs_path = '~/Work/fusOUD/FS';

subjects_list = dir([mask_path '/*T1w.nii.gz']);
% subjects_list = {'sub-224-FUS', 'sub-218-FUS', 'sub-219-FUS'};

LOI = labels_of_interest();

for n = 9 : 9 %length(subjects_list)
    sub_name = subjects_list(n).name(1 : 11);
    fs_path_subject = [fs_path '/FS_' sub_name '_ses-00/mri/aparc+aseg.mgz'];
    system(['mri_convert ' fs_path_subject ' aparc+aseg.nii']);
    system(['labelconvert aparc+aseg.nii $FREESURFER_HOME/FreeSurferColorLUT.txt ' ...
            '$(dirname $(which mrview))/../share/mrtrix3/labelconvert/fs_default.txt ' ...
            'fs_parcels.nii -quiet -force']);

    mask_left_info = niftiinfo([mask_path '/' sub_name '_ses-00_mask_Left_T1.nii.gz']);
    mask_left_vol = niftiread(mask_left_info);

    mask_right_info = niftiinfo([mask_path '/' sub_name '_ses-00_mask_Right_T1.nii.gz']);
    mask_right_vol = niftiread(mask_right_info);

    mask_vol = 42 * mask_left_vol + 49 * mask_right_vol;
    niftiwrite(mask_vol, 'mask.nii', mask_left_info);

    parcels_info = niftiinfo('fs_parcels.nii');
    parcels_vol = niftiread(parcels_info);
    parcels_vol(ismember(parcels_vol, [42 49])) = 0;

    system('python ref_resample.py');
    resampled_mask_vol = niftiread('resampled_mask.nii');

    parcels_vol(resampled_mask_vol > 0) = 0;
    parcels_vol = parcels_vol + resampled_mask_vol;

    niftiwrite(parcels_vol, [mask_path '/' sub_name '_ses-00_parcels_with_mask'], parcels_info, 'Compressed', true);
    disp([mask_path '/' sub_name '_ses-00_parcels_with_mask.nii.gz created!'])
end

delete aparc+aseg.nii fs_parcels.nii resampled_mask.nii mask.nii


%%
function LOI = labels_of_interest()

brainRegions={'caudalanteriorcingulate', 'rostralanteriorcingulate',  'posteriorcingulate' ,  'lateralorbitofrontal', ...
    'medialorbitofrontal' , 'caudalmiddlefrontal'  ,  'rostralmiddlefrontal' ,'frontalpole' , 'insula' , 'parsopercularis', 'parstriangularis', 'parsorbitalis', 'Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Amygdala'}; %, 'Accumbens'};

[labels, names, ~] = xlsread('FS_default_labels.xlsx');

matching = false(size(labels));

for n = 1 : length(brainRegions)
    matching = matching | contains(lower(string(names(:, 2))), lower(brainRegions{n}));
end

LOI = labels(matching);

end
