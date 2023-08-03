clear

data_path = 'FUS/';
subject = 'sub-218-FUS';
session = 'ses-00';

D = 3;
H = 7;

subspot_file = [data_path subject '/' session '/anat/sub-218-FUS_seeds_treatment.csv'];
subspot_data = readmatrix(subspot_file);
subspot_RAS = [-subspot_data(:, 1) -subspot_data(:, 2) subspot_data(:, 3)];

imgfile = [data_path subject '/' session '/anat/' subject '_' session '_treatment.nii.gz'];
info = niftiinfo(imgfile);
mask_volume = single(zeros(info.ImageSize));

origin = info.Transform.T(4, 1:3);
pixle_size = info.Transform.T(1);

sonic_D = D / pixle_size;
sonic_H = H / pixle_size;
subspot_kernel = ones(round(sonic_D), round(sonic_D), round(sonic_H));

subspot_size = size(subspot_kernel);

N = length(subspot_RAS);

for n = 1:N
    ras = subspot_RAS(N, :);
    spot_location = round((ras - origin) / pixle_size - (subspot_size - 1)/2);

    mask_volume(spot_location(1) + 1 : spot_location(1) + subspot_size(1), ...
        spot_location(2) + 1 : spot_location(2) + subspot_size(2), ...
        spot_location(3) + 1 : spot_location(3) + subspot_size(3)) = subspot_kernel;
end

niftiwrite(mask_volume, 'Treat 1 Dose Spot.nii', info)
