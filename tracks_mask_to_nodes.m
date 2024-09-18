clear
close all

cores = 10;
num_seeds = '10k';

data_path = '~/Work/fusOUD/FUS/';
mask_path = '~/Nextcloud/Study/fusOUD/Treatment_Masks';
subject = 'sub-224-FUS';
session = 'ses-00';


mrtrix_path = fullfile(data_path, subject, session, 'dwi', 'mrtrix3');
targeting_path = fullfile(mrtrix_path, 'mask_to_ROIs');

if ~isfolder(mrtrix_path)
    fprintf([subject, '/', session, ' run mrtrix first.\n']);
    return;
end

ref_nii_file = fullfile(mask_path, [subject '_ses-00_T1w.nii.gz']);
mask_nac_left = fullfile(mask_path, [subject '_ses-00_mask_Left_T1.nii.gz']);
mask_nac_right = fullfile(mask_path, [subject '_ses-00_mask_Right_T1.nii.gz']);

% T1 volume in the dwi space
if ~isfile(fullfile(mrtrix_path, 'T1_FS_coreg.nii.gz'))
    fprintf('T1_FS_coreg.nii.gz not found.\n');
    return;
end

if ~isfolder(targeting_path)
    mkdir(targeting_path);
end

% Register the ref volume to T1_FS_coreg.nii.gz
% anat2dwi_0GenericAffine.mat and anat2dwi_Warped.nii.gz will be produced
if ~isfile([targeting_path '/anat2dwi_0GenericAffine.mat'])
    system(['antsRegistrationSyNQuick.sh -d 3 -t r -f ' mrtrix_path '/T1_FS_coreg.nii.gz -m ', ref_nii_file, ' -o ' targeting_path '/anat2dwi_ -n ', num2str(cores)]);
end

apply_rigid_transform(mask_nac_left, [targeting_path '/mask_NAc_left_dwi'], [targeting_path '/anat2dwi_0GenericAffine.mat'])
apply_rigid_transform(mask_nac_right, [targeting_path '/mask_NAc_right_dwi'], [targeting_path '/anat2dwi_0GenericAffine.mat'])

[Labels, Names] = labels_of_interest();

l_mean = zeros(size(Labels));
l_median = zeros(size(Labels));
l_std = zeros(size(Labels));
l_min = zeros(size(Labels));
l_max = zeros(size(Labels));
l_count = zeros(size(Labels));

r_mean = zeros(size(Labels));
r_median = zeros(size(Labels));
r_std = zeros(size(Labels));
r_min = zeros(size(Labels));
r_max = zeros(size(Labels));
r_count = zeros(size(Labels));

for k = 1 : length(Labels)
    include_region_nii = [targeting_path '/include_' Names{k} '.nii.gz'];
    if ~isfile(include_region_nii)
        parcels_info = niftiinfo([mrtrix_path '/fs_parcels_coreg.nii.gz']);
        parcels_vol = niftiread(parcels_info);
        parcels_vol(parcels_vol ~= Labels(k)) = 0;
        parcels_vol(parcels_vol == Labels(k)) = 1;
        niftiwrite(parcels_vol, [targeting_path '/include_' Names{k}], parcels_info, 'Compressed', true);
    end

    tck_file_left = [targeting_path '/seeds_' num_seeds '_left_NAc_to_' Names{k} '.tck'];
    % disp(['Compute ', tck_file_left, '...'])
    cmd = sprintf('tckgen -act %s/5tt_coreg.mif -backtrack -include %s -seed_image %s/mask_NAc_left_dwi.nii.gz %s/wmfod_norm.mif %s -nthreads %d -force -cutoff 0.06 -maxlength 250 -step 0.5 -crop_at_gmwmi -seeds %s -quiet', ...
        mrtrix_path, include_region_nii, targeting_path, mrtrix_path, tck_file_left, cores, num_seeds);
    system(cmd);

    [mean_length, median_length, std_dev, min_length, max_length, streamline_count] = stats_tracks(tck_file_left);
    fprintf('Left NAc to %s\nstats = %f %f %f %f %f %f\n\n', Names{k}, mean_length, median_length, std_dev, min_length, max_length, streamline_count)
    l_mean(k) = mean_length;
    l_median(k) = median_length;
    l_std(k) = std_dev;
    l_min(k) = min_length;
    l_max(k) = max_length;
    l_count(k) = streamline_count;

    tck_file_right = [targeting_path '/seeds_' num_seeds '_right_NAc_to_' Names{k} '.tck'];
    % disp(['Compute ', tck_file_right, '...'])
    cmd = sprintf('tckgen -act %s/5tt_coreg.mif -backtrack -include %s -seed_image %s/mask_NAc_right_dwi.nii.gz %s/wmfod_norm.mif %s -nthreads %d -force -cutoff 0.06 -maxlength 250 -step 0.5 -crop_at_gmwmi -seeds %s -quiet', ...
        mrtrix_path, include_region_nii, targeting_path, mrtrix_path, tck_file_right, cores, num_seeds);
    system(cmd);

    [mean_length, median_length, std_dev, min_length, max_length, streamline_count] = stats_tracks(tck_file_right);
    fprintf('Right NAc to %s\nstats = %f %f %f %f %f %f\n\n', Names{k}, mean_length, median_length, std_dev, min_length, max_length, streamline_count)
    r_mean(k) = mean_length;
    r_median(k) = median_length;
    r_std(k) = std_dev;
    r_min(k) = min_length;
    r_max(k) = max_length;
    r_count(k) = streamline_count;
end

T = table(Names, l_mean, l_median, l_std, l_min, l_max, l_count, r_mean, r_median, r_std, r_min, r_max, r_count); 
writetable(T, [subject '_' session '_streamline_stats.xlsx'])


%% 
function [mean_length, median_length, std_dev, min_length, max_length, streamline_count] = stats_tracks(tck_file)

% Run the command and capture the output
[status, cmdout] = system(['tckstats -quiet ' tck_file]);

% Check if the command was successful
if status == 0
     % Split the output into lines
    cmdout_lines = strsplit(cmdout, '\n');
    
    % The first line is the header, the second line contains the data
    data_line = cmdout_lines{2};
    
    % Convert the data line into numeric values
    data = sscanf(data_line, '%f');
    
    % Assign the values to respective variables
    if ~isempty(data)
        mean_length = data(1);
        median_length = data(2);
        std_dev = data(3);
        min_length = data(4);
        max_length = data(5);
        streamline_count = data(6);
    else
        mean_length = [];
        median_length = [];
        std_dev = [];
        min_length = [];
        max_length = [];
        streamline_count = [];
    end

end

end


%%
function [Labels, Names] = labels_of_interest()

brainRegions={'caudalanteriorcingulate', 'rostralanteriorcingulate',  'posteriorcingulate' ,  'lateralorbitofrontal', ...
    'medialorbitofrontal' , 'caudalmiddlefrontal'  ,  'rostralmiddlefrontal' ,'frontalpole' , 'insula' , 'Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Amygdala'}; %, 'Accumbens'};

[labels, names, ~] = xlsread('util/FS_default_labels.xlsx');

matching = false(size(labels));

for n = 1 : length(brainRegions)
    matching = matching | contains(lower(string(names(:, 2))), lower(brainRegions{n}));
end

Labels = labels(matching);
Names = names(matching, 2);

end