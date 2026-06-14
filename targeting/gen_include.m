clear

data_path = '/mnt/evo/FUS-NV';
subject = 'sub-021-NAV';
session = 'ses-00';

mrtrix_path = fullfile(data_path, subject, session, 'dwi', 'mrtrix3');
out_path = fullfile(data_path, subject, session, 'anat');

if ~isfolder(mrtrix_path)
    fprintf([data_path, '/', subject, '/', session, ' run mrtrix first.\n']);
    return;
end

LOI = labels_of_interest();
parcels_info = niftiinfo([mrtrix_path '/fs_parcels_coreg.nii.gz']);
parcels_vol = niftiread(parcels_info);
parcels_vol(~ismember(parcels_vol, LOI)) = 0;
parcels_vol(ismember(parcels_vol, LOI)) = 1;
niftiwrite(parcels_vol, [out_path '/include_mofc'], parcels_info, 'Compressed', true);

return
%%
function LOI = labels_of_interest()

brainRegions={'medialorbitofrontal'};
[labels, names, ~] = xlsread('FS_default_labels.xlsx');

matching = false(size(labels));

for n = 1 : length(brainRegions)
    matching = matching | contains(lower(string(names(:, 2))), lower(brainRegions{n}));
end

LOI = labels(matching);

end
