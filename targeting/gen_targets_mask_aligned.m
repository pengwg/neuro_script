function gen_targets_mask_aligned(path, r)

disp('Matlab generate targets mask...')
targets_csv = [path '/targets_aligned_RAS.csv'];

data = readmatrix(targets_csv);
RAS = data(:, 1:4);
RAS(:, 4) = 1;

ref_file = [path '/include_spec_aligned.nii.gz'];
info = niftiinfo(ref_file);

IJK = ras2ijk(info.Transform.T, RAS);

mask_vol = single(zeros(info.ImageSize));

sx = size(mask_vol, 1);
sy = size(mask_vol, 2);

d = r;
[dx,dy,dz] = ndgrid(-d:d, -d:d, -d:d);
offsets = dx + dy*sx + dz*sx*sy;

IJK = int32(round(IJK));
idx0 = sub2ind(size(mask_vol), IJK(:, 1), IJK(:, 2), IJK(:, 3));
idx = bsxfun(@plus, idx0, offsets(:)'); 

mask_vol(idx(:)) = 1;
niftiwrite(mask_vol, [path '/targets_mask_aligned'], info, 'Compressed', true)
end

%%
function ijk = ras2ijk(T, ras)
    ijk_h = ras / T;
    ijk = ijk_h(:, 1 : 3) + 1;
end
