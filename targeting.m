clear
close all

cores = 10;

data_path = '/home/peng/Work/fusOUD/FUS';
subject = 'sub-221-FUS';
session = 'ses-07';
session_ref = 'ses-00';
output_path = '~/Nextcloud/Study/fusOUD/OUD221/';

[X, Y, Z] = meshgrid(-25:1:25, -10:1:15, -15:1:15);
target_AC0 = [ 10, 4, -1, 1;
              -10, 4, -1, 1];
r = 0.7;

D = 4;


mrtrix_path = fullfile(data_path, subject, session, 'dwi', 'mrtrix3');
targeting_path = fullfile(mrtrix_path, 'targeting_tracks_ACPC');

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

SUM = single(zeros(size(X)));
SUM_i = 1 : length(SUM(:));

NAc1 = abs(X(:) - target_AC0(1, 1)) <= D & abs(Y(:) - target_AC0(1, 2)) <= D & abs(Z(:) - target_AC0(1, 3)) <= D;
NAc2 = abs(X(:) - target_AC0(2, 1)) <= D & abs(Y(:) - target_AC0(2, 2)) <= D & abs(Z(:) - target_AC0(2, 3)) <= D;
NAc = NAc1 | NAc2;
SUM_i = SUM_i(NAc);

target_AC = [X(NAc), Y(NAc), Z(NAc), ones(size(X(NAc)))];
target_dwi = target_AC * AC2dwi_transform';

LOI = labels_of_interest();
parcels_info = niftiinfo([mrtrix_path '/fs_parcels_coreg.nii.gz']);
parcels_vol = niftiread(parcels_info);
parcels_vol(~ismember(parcels_vol, LOI)) = 0;
parcels_vol(ismember(parcels_vol, LOI)) = 1;
niftiwrite(parcels_vol, [targeting_path '/include_spec'], parcels_info);

for k = 1 : length(X(NAc))
    x0 = target_AC(k, 1);
    y0 = target_AC(k, 2);
    z0 = target_AC(k, 3);

    x = target_dwi(k, 1);
    y = target_dwi(k, 2);
    z = target_dwi(k, 3);

    % tck_file = sprintf('%s/seeds_1000_AC_%d_%d_%d_RAS_%d_%d_%d.tck', targeting_path, x0, y0, z0, round(x), round(y), round(z));
    tck_file = sprintf('%s/select_100_AC_%d_%d_%d_RAS_%d_%d_%d.tck', targeting_path, x0, y0, z0, round(x), round(y), round(z));

    disp(['Compute ', tck_file, '...'])
    if ~isfile(tck_file)
        % cmd = sprintf('tckgen -act %s/5tt_coreg.mif -backtrack -include %s/include_spec.nii -seed_sphere %f,%f,%f,%f %s/wmfod_norm.mif %s -nthreads %d -force -cutoff 0.06 -maxlength 250 -step 0.5 -crop_at_gmwmi -seeds 1000 -quiet', ...
        %     mrtrix_path, targeting_path, x, y, z, r, mrtrix_path, tck_file, cores);
        cmd = sprintf('tckgen -act %s/5tt_coreg.mif -backtrack -include %s/include_spec.nii -seed_sphere %f,%f,%f,%f %s/wmfod_norm.mif %s -nthreads %d -force -cutoff 0.06 -maxlength 250 -step 0.5 -crop_at_gmwmi -select 100 -quiet', ...
            mrtrix_path, targeting_path, x, y, z, r, mrtrix_path, tck_file, cores);

        system(cmd);
    end
    
    num_tracks = count_tracks(tck_file);
    fprintf('num_tracks = %d\n', num_tracks);

    sift_file = sprintf('%s/sift2_select_100_AC_%d_%d_%d_RAS_%d_%d_%d.txt', targeting_path, x0, y0, z0, round(x), round(y), round(z));
    disp(['Compute ', sift_file, '...'])

    if num_tracks > 1
        if ~isfile(sift_file)
            cmd = sprintf('tcksift2 -act %s/5tt_coreg.mif %s %s/wmfod_norm.mif %s -quiet -force -nthreads %d', mrtrix_path, tck_file, mrtrix_path, sift_file, cores);
            system(cmd);
        end        
        sift_sum = sift2_weights_sum(sift_file);
    else
        sift_sum = 0;
    end
 
    fprintf('sift_sum = %d\n\n', sift_sum);
    SUM(SUM_i(k)) = sift_sum;
end
SUM = permute(SUM, [2, 1, 3]);

ref_info = niftiinfo(ref_nii_file);
ref_info.Datatype = 'single';
ref_info.PixelDimensions = [1, 1, 1];
ref_info.Transform.T = eye(4);
ref_info.ImageSize = size(SUM);
ref_info.Transform.T(4, 1:3) = [X(1) Y(1) Z(1)];

if ~isfolder(output_path)
    mkdir(output_path);
end

niftiwrite(SUM, [output_path, 'tck_sift_07'], ref_info);

figure
imagesc(rot90(SUM(:,:,7)', 2))
axis equal
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


%%
function LOI = labels_of_interest()

% brainRegions={'caudalanteriorcingulate', 'rostralanteriorcingulate',  'posteriorcingulate' ,  'lateralorbitofrontal', ...
%     'medialorbitofrontal' , 'caudalmiddlefrontal'  ,  'rostralmiddlefrontal' ,'frontalpole' , 'insula' , 'Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Amygdala'}; %, 'Accumbens'};
brainRegions={'insula'};

[labels, names, ~] = xlsread('util/FS_default_labels.xlsx');

matching = false(size(labels));

for n = 1 : length(brainRegions)
    matching = matching | contains(lower(string(names(:, 2))), lower(brainRegions{n}));
end

LOI = labels(matching);

end

%%
function sift_sum = sift2_weights_sum(sift_file)


fileID = fopen(sift_file, 'r');

% Skip the first line (header line)
fgetl(fileID);  
secondLine = fgetl(fileID);

fclose(fileID);

values = str2double(strsplit(secondLine));

sift_sum = sum(values);

end
