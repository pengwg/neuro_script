clear
close all

cores = 12;

data_path = '/home/pw0032/Data/FUS-RCT';
subject = 'sub-018-RCT';
session = 'ses-2-00';
session_ref = 'ses-2-00';

target_AC = [ 9.5, 4, -1, 1;
              -10, 3, -2, 1];

mrtrix_path = fullfile(data_path, subject, session, 'dwi', 'mrtrix3');
targeting_path = fullfile(mrtrix_path, 'targeting_tracks_AC');
tckmap_path = fullfile(mrtrix_path, 'tckmap');

if ~isfolder(mrtrix_path)
    fprintf([data_path, '/', subject, '/', session, ' run mrtrix first.\n']);
    return;
end

% Always use the reference volume and seed files from ses-1-00
if isfolder(fullfile(data_path, subject, session_ref, 'anat'))
    ref_nii = dir(fullfile(data_path, subject, session_ref, 'anat', [subject, '_', session_ref, '_AC.nii*']));
end

if isempty(ref_nii)
    fprintf([subject, '_ses-00_AC.nii(.gz) not found.\n']);
    return;
end

ref_nii_file = fullfile(ref_nii.folder, ref_nii.name);

% T1 volume in the dwi space
if ~isfile(fullfile(mrtrix_path, 'T1_FS_coreg.nii.gz'))
    fprintf('T1_FS_coreg.nii.gz not found.\n');
    return;
end

if ~isfolder(tckmap_path)
    mkdir(tckmap_path);
end

if ~isfolder(targeting_path)
    mkdir(targeting_path);
end

% Register the ref volume to T1_FS_coreg.nii.gz
% AC2dwi_0GenericAffine.mat and AC2dwi_Warped.nii.gz will be produced
if ~isfile([targeting_path '/AC2dwi_0GenericAffine.mat'])
    system(['antsRegistrationSyNQuick.sh -d 3 -t r -f ' mrtrix_path '/T1_FS_coreg.nii.gz -m ', ref_nii_file, ' -o ' targeting_path '/AC2dwi_ -n ', num2str(cores)]);
end

AC2dwi_transform_data = load([targeting_path '/AC2dwi_0GenericAffine.mat']);

AC2dwi_transform = ea_antsmat2mat(AC2dwi_transform_data.AffineTransform_double_3_3, ...
                                     AC2dwi_transform_data.fixed);

target_dwi = target_AC * AC2dwi_transform';

LOI = labels_of_interest();
parcels_info = niftiinfo([mrtrix_path '/fs_parcels_coreg.nii.gz']);
parcels_vol = niftiread(parcels_info);
parcels_vol(~ismember(parcels_vol, LOI)) = 0;
parcels_vol(ismember(parcels_vol, LOI)) = 1;
niftiwrite(parcels_vol, [tckmap_path '/include_spec'], parcels_info);

%% Create streamlines from large regions covering targets
x = target_dwi(1:2, 1);
y = target_dwi(1:2, 2);
z = target_dwi(1:2, 3);

r = 2.5;
select = '100k';
tck_file = sprintf('%s/sphere_%d_select_%s_NAc.tck', tckmap_path, r, select);
disp(['Compute ', tck_file, '...'])
cmd = sprintf(['tckgen -act %s/5tt_coreg_hsvs.mif -backtrack -include %s/include_spec.nii ' ...
    '-seed_sphere %f,%f,%f,%f -seed_sphere %f,%f,%f,%f ' ...
    '%s/wmfod_norm.mif %s -select %s -cutoff 0.06 -maxlength 250 ' ...
    '-step 0.5 -crop_at_gmwmi -nthreads %d -quiet -force'], ...
    mrtrix_path, tckmap_path, x(1), y(1), z(1), r, x(2), y(2), z(2), r, mrtrix_path, tck_file, select, cores);

system(cmd);

% return

%% Create streamlines from targets
r = 2.5;
seeds = '3k';

[X, Y, Z] = meshgrid(-1:1:1, -1:1:1, -1:1:1);
targets_AC = repmat(target_AC, length(X(:)), 1);
X = [X(:)'; X(:)'];
Y = [Y(:)'; Y(:)'];
Z = [Z(:)'; Z(:)'];
targets_AC(:, 1:3) = targets_AC(:, 1:3) + [X(:), Y(:), Z(:)];

targets_dwi = targets_AC * AC2dwi_transform';

for k = 1 : size(targets_dwi, 1)
    x0 = targets_AC(k, 1);
    y0 = targets_AC(k, 2);
    z0 = targets_AC(k, 3);

    x = targets_dwi(k, 1);
    y = targets_dwi(k, 2);
    z = targets_dwi(k, 3);

    tck_file = sprintf('%s/seeds_%s_AC_%2.1f_%2.1f_%2.1f_RAS_%2.1f_%2.1f_%2.1f.tck', targeting_path, seeds, x0, y0, z0, round(x), round(y), round(z));
    disp(['Compute ', tck_file, '...'])
    cmd = sprintf(['tckgen -act %s/5tt_coreg_hsvs.mif -backtrack -include %s/include_spec.nii ' ...
        '-seed_sphere %f,%f,%f,%f %s/wmfod_norm.mif %s -seeds %s -cutoff 0.04 -maxlength 250 ' ...
        '-step 0.5 -crop_at_gmwmi -nthreads %d -quiet -force'], mrtrix_path, tckmap_path, x, y, z, r, mrtrix_path, tck_file, seeds, cores);

    % cmd = sprintf(['tckgen -act %s/5tt_coreg_hsvs.mif -backtrack ' ...
    %     '-seed_sphere %f,%f,%f,%f %s/wmfod_norm.mif %s -select %s -cutoff 0.04 -maxlength 250 ' ...
    %     '-step 0.5 -crop_at_gmwmi -nthreads %d -quiet -force'], mrtrix_path, x, y, z, r, mrtrix_path, tck_file, seeds, cores);
    system(cmd);
end


%%
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

