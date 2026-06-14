clear
close all

cores = 12;

data_path = "/mnt/evo/FUS-NV";
subject = "sub-019-NAV";
session = "ses-00";

target_AC = [9.2, 2.5, 0;
            -9.5, 3.0, -0.5];

session_path = fullfile(data_path, subject, session);
mrtrix_path = fullfile(session_path, "dwi", "mrtrix3");
targeting_path = fullfile(session_path, "targeting");

if ~isfolder(mrtrix_path)
    fprintf(session_path + " run mrtrix first.\n");
    return;
end

if ~isfolder(targeting_path)
    mkdir(targeting_path)
end

if ~isfile(fullfile(mrtrix_path, "fs_parcels_aligned.nii.gz"))
    disp("fs_parcels_aligned.nii.gz not found!");
    return;
end

if ~isfolder(targeting_path)
    mkdir(targeting_path);
end

LOI = labels_of_interest();
parcels_info = niftiinfo(fullfile(mrtrix_path, "fs_parcels_aligned.nii.gz"));
parcels_vol = niftiread(parcels_info);
parcels_vol(~ismember(parcels_vol, LOI)) = 0;
parcels_vol(ismember(parcels_vol, LOI)) = 1;
niftiwrite(parcels_vol, fullfile(targeting_path, "include_spec_aligned"), parcels_info, "Compressed", true);

%% Create streamlines from targets
r = 2.5;
step = -0.2 : 0.2 : 0.2;

[X, Y, Z] = meshgrid(step, step, step);
targets_AC = repmat(target_AC, length(X(:)), 1);
X = [X(:)'; X(:)'];
Y = [Y(:)'; Y(:)'];
Z = [Z(:)'; Z(:)'];
targets_AC(:, 1:3) = targets_AC(:, 1:3) + [X(:), Y(:), Z(:)];

for k = 1 : size(targets_AC, 1)
    R = targets_AC(k, 1);
    A = targets_AC(k, 2);
    S = targets_AC(k, 3);

    if R < 0
        track_file = sprintf('Left_r%2.1f_a%2.1f_s%2.1f_r%2.1f.tck', R, A, S, r);
    else
        track_file = sprintf('Right_r%2.1f_a%2.1f_s%2.1f_r%2.1f.tck', R, A, S, r);
    end
    disp(['Compute ', track_file, ' ...'])

    tckedit_cmd = "tckedit " + fullfile(mrtrix_path, "tracks_10M_aligned.tck") + ' ' + ...
        fullfile(targeting_path, track_file) + ' ' + ...
        "-include " + R + ',' + A + ',' + S + ',' + r + ' ' + ...
        "-include " + fullfile(targeting_path, "include_spec_aligned.nii.gz") + ' ' + ...
        "-nthreads 8 -force";

    [~, ~] = system(tckedit_cmd);
end


%%
function LOI = labels_of_interest()

brainRegions={'caudalanteriorcingulate', 'rostralanteriorcingulate',  'posteriorcingulate' ,  'lateralorbitofrontal', ...
     'medialorbitofrontal' , 'caudalmiddlefrontal'  ,  'rostralmiddlefrontal' ,'frontalpole' , 'insula' , 'Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Amygdala'}; %, 'Accumbens'};
% brainRegions={'insula'};

[labels, names, ~] = xlsread('util/FS_default_labels.xlsx');

matching = false(size(labels));

for n = 1 : length(brainRegions)
    matching = matching | contains(lower(string(names(:, 2))), lower(brainRegions{n}));
end

LOI = labels(matching);

end

