

function [Y, coords, mask] = superimpose_simple(objimg, mask, topLeft, im)

disp('-- Superimposing object...')
coords = [];
% crop img,mask to be tight on the object
[r,c] = find(mask>0);
left = min(c); right = max(c);
up = min(r); down = max(r);
objimg = objimg(up:down, left:right,:); 
mask = mask(up:down, left:right);

%figure; imagesc(objimg);
%figure; imagesc(mask);

% get only the single large connected component from the mask
%mask = eliminate_small_components(mask);

% use resized img and mask to segment the object
seg = segment_object(objimg, mask);
%figure; imagesc(uint8(seg));
%imwrite(uint8(seg), ['gmu_kitch_empty/mask_examples/mask_11.png']);

pw = size(seg,1); ph = size(seg,2);

bottom=topLeft(1)+pw-1; right=topLeft(2)+ph-1;
if bottom>size(im,1) || right>size(im,2), Y=0; return; end

patch = im(topLeft(1):bottom, topLeft(2):right, :);
%figure; imagesc(patch);

inds = find(mask>0);
patch_r = patch(:,:,1); patch_g = patch(:,:,2); patch_b = patch(:,:,3);
patch_r(inds)=0; patch_g(inds)=0; patch_b(inds)=0;
patch(:,:,1)=patch_r; patch(:,:,2)=patch_g; patch(:,:,3)=patch_b;
patch = patch + uint8(seg);
%figure; imagesc(patch);

im(topLeft(1):bottom, topLeft(2):right, :) = patch;
% create the world mask for this object
%imMask = zeros(size(im,1), size(im,2), 1);
%imMask(topLeft(1):bottom, topLeft(2):right) = mask;

Y = im;
coords = [topLeft(1) topLeft(2) bottom right];
