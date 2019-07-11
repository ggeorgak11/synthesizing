close all; clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used to generate the synthetic sets for the paper:
% Synthesizing Training Data for Object Detection in Indoor Scenes.
% Set-up to work for kitchen scenes and generates both a blended
% and a simple superimposed dataset.
% Georgios Georgakis 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize the parameters, information inside the script
init_params;

% create the params for the planes segmentation
plane_params = get_plane_params;

scene_dirs = dir([params.back_dir, 'kitchen_*']);
a=[scene_dirs.isdir];
scene_dirs = scene_dirs(a); % to get only the scene folders

for i=1:size(scene_dirs,1)
    scene = scene_dirs(i).name;
    img_files = dir([params.back_dir, scene, '/*im.png']);
    
    save_scene_blended_dir = [params.save_dir_blended, scene, '/'];
    save_scene_imposed_dir = [params.save_dir_imposed, scene, '/'];
    for j=1:size(img_files,1)
        path_dir = [params.back_dir, scene, '/'];
        img_name = img_files(j).name;
        params.img_id = img_name(1:end-4);
        %if exist([save_scene_blended_dir, params.imgs_dir, img_files(j).name], 'file'), continue; end
        
        depth_name = strrep(img_name, 'im.png', 'depth.mat');
        seg_name = strrep(img_name, 'im.png', 'seg.mat');
        
        im = imread([path_dir, img_name]);
        load([path_dir, depth_name]); % loads depth_filled
        load([path_dir, seg_name]); % loads segmentation
        
        seg_size = size(segmentation);
        im=imresize(im, seg_size); depth=imresize(depth_filled, seg_size);
        depth=depth*1e3;
       
        %show_seg(im, depth, segmentation, className); % uncomment to visualize the segmentation 
        
        sampling_area = get_sampling_area(im, depth, segmentation, seg_size, params, plane_params);       
        %figure; imagesc(sampling_area);
        if sum(sum(sampling_area)) < 1, continue; end % check if there is somewhere to sample
        
        % create k blendings of the same image
        for m=1:params.k
            im_name = img_files(j).name;
            im_name = [im_name(1:end-4), '_', num2str(m), im_name(end-3:end)];
            
            %%% blending process
            n=randi(params.min_obj,1,1)+(params.max_obj-params.min_obj); % from min_obj to max_obj objects
            [im_blended, im_imposed, info] = blend_with_seg(im, depth, n, seg_size, sampling_area, params);
            %%%
            if isempty(info), continue; end % if no objects were blended then dont save the image           
            %%% save the annotation struct of the img
            if ~exist([save_scene_blended_dir, params.imgs_dir], 'dir'), mkdir([save_scene_blended_dir, params.imgs_dir]); end
            if ~exist([save_scene_imposed_dir, params.imgs_dir], 'dir'), mkdir([save_scene_imposed_dir, params.imgs_dir]); end
            if ~exist([save_scene_blended_dir, params.annot_dir], 'dir'), mkdir([save_scene_blended_dir, params.annot_dir]); end
            if ~exist([save_scene_imposed_dir, params.annot_dir], 'dir'), mkdir([save_scene_imposed_dir, params.annot_dir]); end

            disp(['===> Writing ', im_name, ' ...']);

            % write the synthetic image
            imwrite(uint8(im_blended), [save_scene_blended_dir, params.imgs_dir, im_name]);
            imwrite(uint8(im_imposed), [save_scene_imposed_dir, params.imgs_dir, im_name]);
            % write the annotations mat file
            save([save_scene_blended_dir, params.annot_dir, strrep(im_name, '.png', '.mat')], 'info');
            save([save_scene_imposed_dir, params.annot_dir, strrep(im_name, '.png', '.mat')], 'info');
            % create the xml for this synthetic image
            create_xml(info, [save_scene_blended_dir, params.annot_dir], strrep(im_name,'.png',''), scene, seg_size);
            create_xml(info, [save_scene_imposed_dir, params.annot_dir], strrep(im_name,'.png',''), scene, seg_size); 
        end
        
    end
    
end

% save the params for the current exp
save([params.save_dir_blended, 'params.mat'], 'params');
save([params.save_dir_imposed, 'params.mat'], 'params');

