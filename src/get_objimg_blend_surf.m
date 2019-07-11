

function [objimg_to_blend, mask_to_blend, top, left, objparams] = get_objimg_blend_surf(seg_size, depth, umap, params, ind)

% choose the bottom right position and then find top left, instead of
% choosing top left directly. This is to ensure that the objects are on top
% of the support surface

% init outputs
objimg_to_blend=[];
mask_to_blend=[];
top=0; left=0;
objparams=[];

% choose class, pose randomly
rClass = randi(length(params.object_list),1,1); % random object label
objpath = [params.objroot, params.object_list{rClass}];
% randomly select object image, implicitly selects azimuth   
rIndex = params.index_seq(randi(length(params.index_seq), 1, 1));
% choose random cam(1-3), exclude cam 4,5, implicitly selects elevation
rCam = randi(3,1,1);
% read img and mask, read depth also
objimg = imread([objpath, '/NP', num2str(rCam), '_', num2str(rIndex), '.jpg']);
hdepthinfo = hdf5info([objpath, '/NP', num2str(rCam), '_', num2str(rIndex),'.h5']);
objdepth = hdf5read(hdepthinfo.GroupHierarchy.Datasets(1));
objdepth = imrotate(objdepth, -90); objdepth = flipdim(objdepth,2); 
objdepth = double(objdepth); %depth = imresize(depth, [size(objimg,1) size(objimg,2)]);
% figure; imagesc(objimg);
% figure; imagesc(objdepth);
% figure; imagesc(depth);
% keyboard 

% choose the masks from graphcut if refined_masks_flag==1
%if refined_masks_flag
    mask_path = [params.dataroot, 'object_data/our_objects_masks/', params.object_list{rClass}, '/NP', num2str(rCam), '_', num2str(rIndex), '_mask_refined.png'];
    if ~exist(mask_path, 'file'), return; end;
    mask = imread(mask_path);
%else
%    mask = imread([objpath, '/masks/NP', num2str(rCam), '_', num2str(rIndex), '_mask.pbm']);
%    mask=abs(1-mask);
%end     

% save the object files
index = num2str(ind);
if params.save_imgs
    imwrite(objimg, [params.save_path, params.img_id, '_obj_img_', index,'.png']);
    %imwrite(uint16(objdepth), [params.save_path, params.img_id, '_obj_depth_', index,'.png']);
    f=figure; imagesc(objdepth); set(gcf,'Visible', 'off');
    set(gca,'LooseInset',get(gca,'TightInset')); axis off;
    saveas(f, [params.save_path, params.img_id, '_obj_depth', index,'.png']);
    
    imwrite(mask, [params.save_path, params.img_id, '_obj_mask_', index,'.png']);
end


if params.scaling_mode==0
    % randomly choose a scale for the object
    scale = params.obj_scaling_list(randi(length(params.obj_scaling_list), 1, 1));
    % resize the img and mask
    objimg = imresize(objimg, scale);
    mask = imresize(mask, scale);
    mask(mask<0.5)=0; mask(mask>=0.5)=1;
    
    % get mask dimensions so to choose bottom right that fits the object
    [r,c]=find(mask);
    height=max(r)-min(r); width=max(c)-min(c);

    % randomly choose the topLeft coordinate
    %maxh = seg_size(1)-height-1; % maximum available height to choose from (valid range 1...im size)
    %maxW = seg_size(2)-width-1; % maximum available width to choose from
    minh = height+1; % choose from (minh:seg_size(1)) so that object fits in the image
    minw = width+1;
    
    % need to choose the position based on the umap that is given
    %[r,c]=find(umap);
    % if umaps maximum positions are below the mins then we cannot put object
    %if max(r)<minh || max(c)<minw, return; end
    % choose a position in umap that is larger than the mins
    % remove the invalid positions
    umap(1:minh,:)=0; umap(:,1:minw)=0;
    if sum(sum(umap)) == 0, return; end
    [rv, cv]=find(umap);
    %try
        ind = randi([1 length(rv)], 1, 1);
    %catch
    %    keyboard
    %end
    bottom=rv(ind); right=cv(ind);
    %rv=r(find(r>minh)); cv=c(find(c>minw));
    %ind = randi([1 length(rv)], 1, 1);
    %try
    %    bottom=rv(ind); right=cv(ind);
    %catch
    %    keyboard
    %end
        
    top = bottom - height;
    left = right - width;
    objimg_to_blend = objimg;
    mask_to_blend = mask;
else
    % scale of object should be chosen by depth ratio
    % sample depth from the objimg. To sample this properly and make
    % sense when finding the ratio in the image, we have to resize
    % both depth and mask to the segmentation size
    mask_tmp = imresize(mask, seg_size);
    mask_tmp(mask_tmp<0.5*256)=0; mask_tmp(mask_tmp>=0.5*256)=1;
    depth_tmp = imresize(objdepth, seg_size);
    objimg_tmp = imresize(objimg, seg_size);

    [r,c]=find(mask_tmp);
    depth_sample = depth_tmp(min(r):max(r), min(c):max(c)); %figure; imagesc(depth_sample);
    depth_sample = reshape(depth_sample, [size(depth_sample,1)*size(depth_sample,2) 1]);
    depth_sample = depth_sample(find(depth_sample>0)); % remove zeros
    smp = median(depth_sample)/10; % scene sample will be in mm  
    height_tmp = max(r)-min(r); width_tmp = max(c)-min(c); % keep mask_tmp size in pixels

    if isnan(smp), return; end % problem with the hersheys bar!, no depth values 
    
    % sample a position for object bottom right in umap
    % keep sampling until the position is valid
    % valid position is one that contains more than per of the other props
    valid=0; iter=0;
    while valid==0
        %disp(iter);
        if iter>100, break; end
        iter=iter+1;

        [r,c] = find(umap);
        ind = randi([1 length(r)], 1, 1);
        %c_y=r(ind); c_x=c(ind);
        bottom=r(ind); right=c(ind);
        % get the center

        % get a depth sample from the image
        %test_smp = depth(c_y, c_x);
        %test_smp = sample_depth(depth, [c_x c_y]);
        test_smp = sample_depth(depth, [right bottom]);
        if test_smp==0, continue; end % if test_smp is 0 then invalid

        ratio = smp/double(test_smp);
        
        % do small scaling pertubations, randomly select to slightly (0.8-1.2) distort the scale
        if params.scaling_mode==2
            pertubation = 0.8 + (1.2-0.8)*rand;
            ratio = ratio*pertubation;
        end
        
        height_on_img = round(height_tmp*ratio); 
        width_on_img = round(width_tmp*ratio);
        % find the top left coordinate
        %smp_top = round(c_y-height_on_img/2);
        %smp_left = round(c_x-width_on_img/2);
        %bottom = c_y+height_on_img;
        %right = c_x+width_on_img;
        smp_top = bottom - height_on_img;
        smp_left = right - width_on_img;
        if smp_top<1 || smp_left<1 || bottom>size(umap,1) || right>size(umap,2) || isnan(smp_top) || isnan(smp_left)
            continue; 
        end
        
        % problem with the hersheys bar!, no depth values 
        try
            umap_patch = umap(smp_top:bottom, smp_left:right);
        catch
           keyboard 
        end
            %lbls=unique(umap_patch);
        other_props_ratio = length(find(umap_patch))/(size(umap_patch,1)*size(umap_patch,2));
        if other_props_ratio >= params.per, valid=1; end
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
end
    
objparams.rClass = rClass;
objparams.rCam = rCam;
objparams.rIndex = rIndex;


