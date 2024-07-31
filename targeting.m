clear
close all

cores = 10;

data_path = '/home/pw0032/Work/fusOUD/FUS-RCT/';
subject = 'sub-006-RCT';
session = 'ses-1-00';

num_tracks = '1k';

target_AC0 = [ 11, 4, 0, 1];
             % -9, 3, 0, 1];
r = 0.5;

% Using target_seeds_${seeds_tag}.csv
seeds_tag = 'ACPC';

mrtrix_path = fullfile(data_path, subject, session, 'dwi', 'mrtrix3');
targeting_path = fullfile(mrtrix_path, 'targeting_tracks_ACPC');

if ~isfolder(mrtrix_path)
    fprintf([subject, '/', session, ' run mrtrix first.\n']);
    return;
end

% Always use the reference volume and seed files from ses-1-00
if isfolder(fullfile(data_path, subject, 'ses-1-00', 'anat'))
    ref_nii = dir(fullfile(data_path, subject, 'ses-1-00', 'anat', [subject, '_ses-1-00_AC.nii*']));
end

if isempty(ref_nii)
    fprintf([subject, '_ses-1-00_AC.nii(.gz) not found.\n']);
    return;
end

ref_nii_file = fullfile(ref_nii.folder, ref_nii.name);

% T1 volume in the dwi space
if ~isfile(fullfile(mrtrix_path, 'T1_FS_coreg.nii.gz'))
    fprintf('T1_FS_coreg.nii.gz not found.\n');
    return;
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

[dx, dy] = meshgrid(-3:1:3, 3:-1:-3);
count = zeros(size(dx));

for k = 1 : length(dx(:))
    target_AC = target_AC0 + [dx(k), dy(k), 0, 0];

    target_dwi = target_AC * AC2dwi_transform';

    for i = 1 : size(target_dwi, 1)
        x0 = target_AC(i, 1);
        y0 = target_AC(i, 2);
        z0 = target_AC(i, 3);

        x = target_dwi(i, 1);
        y = target_dwi(i, 2);
        z = target_dwi(i, 3);

        tck_file = sprintf('%s/seeds_2000_AC_%d_%d_%d_RAS_%d_%d_%d.tck', targeting_path, x0, y0, z0, round(x), round(y), round(z));
        cmd = sprintf('tckgen -act %s/5tt_coreg.mif -backtrack -seed_sphere %f,%f,%f,%f %s/wmfod_norm.mif %s -nthreads %d -force -cutoff 0.06 -maxlength 250 -step 0.5 -crop_at_gmwmi -seeds 1000', ...
            mrtrix_path, x, y, z, r, mrtrix_path, tck_file, cores);
        system(cmd);
        count(k) = count_tracks(tck_file);
    end

end

figure
imagesc(count)
colormap jet
colorbar

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
function count = count_tracks(tck_file)

% Run the command and capture the output
[status, cmdout] = system(['tckinfo ' tck_file]);

% Check if the command was successful
if status == 0
    % Find the line containing 'count:'
    lines = strsplit(cmdout, '\n');
    count_line = '';
    for i = 1:length(lines)
        if contains(lines{i}, 'count:') && ~contains(lines{i}, 'total_count:')
            count_line = lines{i};
            break;
        end
    end
    
    % Extract the count value
    if ~isempty(count_line)
        count_value = regexp(count_line, '\d+', 'match');
        count = str2double(count_value{1});
    else
        error('count not found in tckinfo output.');
    end
else
    error('Failed to execute tckinfo command.');
end

end
