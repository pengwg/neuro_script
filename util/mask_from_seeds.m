clear

data_path = 'FUS/';
subject = 'sub-214-FUS';

% Each seed will be represented as a square cuboid with the size DxDxH mm in the mask volume.
D = 2.5;
H = 7;

seed_file = [data_path subject '/ses-00/anat/' subject '_seeds_treatment.csv'];
seed_data = readmatrix(seed_file);
seed_RAS = [-seed_data(:, 1) -seed_data(:, 2) seed_data(:, 3)];

imgfile = [data_path subject '/ses-00/anat/' subject '_ses-00_treatment.nii.gz'];
info = niftiinfo(imgfile);
mask_volume = cast(zeros(info.ImageSize), info.Datatype);

origin = info.Transform.T(4, 1:3);
pixel_size = diag(info.Transform.T)';

sonic_size = [D / pixel_size(1) D / pixel_size(2) H / pixel_size(3)];
seed_kernel = ones(round(sonic_size));

seed_size = size(seed_kernel);

N = size(seed_RAS, 1);

for n = 1 : N
    ras = seed_RAS(n, :);
    spot_location = round((ras - origin) ./ pixel_size(1:3) - (seed_size - 1)/2);

    mask_volume(spot_location(1) + 1 : spot_location(1) + seed_size(1), ...
        spot_location(2) + 1 : spot_location(2) + seed_size(2), ...
        spot_location(3) + 1 : spot_location(3) + seed_size(3)) = seed_kernel;
end

niftiwrite(mask_volume, [data_path subject '/ses-00/anat/' subject '_ses-00_mask'], info, 'Compressed', true)
