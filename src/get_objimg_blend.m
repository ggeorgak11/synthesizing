

function [objimg_to_blend, mask_to_blend, top, left, objparams] = get_objimg_blend(objects, objroot, root_dir, segH, segW, depth, umap, per)

% init outputs
objimg_to_blend=[];
mask_to_blend=[];
top=0; left=0;
objparams=[];


index_seq = 0:3:357;

% choose class, pose randomly
rClass = randi(length(objects),1,1); % random object label
objpath = [objroot, objects{rClass}, '/'];
% randomly select object image, implicitly selects azimuth   
rIndex = index_seq(randi(length(index_seq), 1, 1));
% choose random cam(1-3), exclude cam 4,5, implicitly selects elevation
rCam = randi(3,1,1);
% read img and mask, read depth also
objimg = imread([objpath, '/NP', num2str(rCam), '_', num2str(rIndex), '.jpg']);
hdepthinfo = hdf5info([objpath, '/NP', num2str(rCam), '_', num2str(rIndex),'.h5']);
objdepth = hdf5read(hdepthinfo.GroupHierarchy.Datasets(1));
objdepth = imrotate(objdepth, -90); objdepth = flipdim(objdepth,2); 
objdepth = double(objdepth); %depth = imresize(depth, [size(objimg,1) size(objimg,2)]);
%figure; imagesc(objimg)
%figure; imagesc(depth);
%keyboard 

% choose the masks from graphcut if refined_masks_flag==1
%if refined_masks_flag
    mask_path = [root_dir, 'object_data/our_objects_masks/', objects{rClass}, '/NP', num2str(rCam), '_', num2str(rIndex), '_mask_refined.png'];
    if ~exist(mask_path, 'file'), return; end;
    mask = imread(mask_path);
%else
%    mask = imread([objpath, '/masks/NP', num2str(rCam), '_', num2str(rIndex), '_mask.pbm']);
%    mask=abs(1-mask);
%end     

% scale of object should be chosen by depth ratio
% sample depth from the objimg. To sample this properly and make
% sense when finding the ratio in the image, we have to resize
% both depth and mask to the segmentation size
mask_tmp = imresize(mask, [segH segW]);
mask_tmp(mask_tmp<0.5*256)=0; mask_tmp(mask_tmp>=0.5*256)=1;
depth_tmp = imresize(objdepth, [segH segW]);
objimg_tmp = imresize(objimg, [segH segW]);

[r,c]=find(mask_tmp);
depth_sample = depth_tmp(min(r):max(r), min(c):max(c)); %figure; imagesc(depth_sample);
depth_sample = reshape(depth_sample, [size(depth_sample,1)*size(depth_sample,2) 1]);
depth_sample = depth_sample(find(depth_sample>0)); % remove zeros
smp = median(depth_sample)/10; % sample is in mm  
height_tmp = max(r)-min(r); width_tmp = max(c)-min(c); % keep mask_tmp size in pixels

% sample a position for object center in umap
% keep sampling until the position is valid
% valid position is one that contains more than per of the other props
valid=0; iter=0;
while valid==0
    %disp(iter);
    if iter>100, break; end
    iter=iter+1;
    
    [r,c] = find(umap);
    ind = randi([1 length(r)], 1, 1);
    c_y=r(ind); c_x=c(ind);
    % get a depth sample from the image
    %test_smp = depth(c_y, c_x);
    test_smp = sample_depth(depth, [c_x c_y]);
    if test_smp==0, continue; end % if test_smp is 0 then invalid

    ratio = smp/double(test_smp); 
    height_on_img = round(height_tmp*ratio); 
    width_on_img = round(width_tmp*ratio);
    % find the top left coordinate
    smp_top = round(c_y-height_on_img/2);
    smp_left = round(c_x-width_on_img/2);
    bottom = c_y+height_on_img;
    right = c_x+width_on_img;
    if smp_top<1 || smp_left<1 || bottom>size(umap,1) || right>size(umap,2), continue; end

    umap_patch = umap(smp_top:bottom, smp_left:right);
    %lbls=unique(umap_patch);
    other_props_ratio = length(find(umap_patch))/(size(umap_patch,1)*size(umap_patch,2));
    if other_props_ratio >= per, valid=1; end

    %figure; imagesc(segmentation); hold on;
    %rectangle('Position',[left top width_on_img height_on_img],'EdgeColor','r', 'LineWidth',3);
    %hold off;
end

if valid==0, return; end % never chosen a valid position

% need to rescale the mask and the objimg based on the ratio found
objimg_to_blend = imresize(objimg_tmp, ratio);
mask_to_blend = imresize(mask_tmp, ratio);
mask_to_blend(mask_to_blend<0.5)=0; mask_to_blend(mask_to_blend>=0.5)=1;
top = smp_top;
left = smp_left;

objparams.rClass = rClass;
objparams.rCam = rCam;
objparams.rIndex = rIndex;


