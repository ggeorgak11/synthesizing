close all; clear all; clc;

% Process the raw data as downloaded from the NYU v2 webpage.
% The script synchronizes and aligns the raw data.
% Need to do this before generating a synthetic set.

% path for the nyu v2 toolbox
% save the aligned video frames 
dataroot = '/project_root/';
addpath([dataroot, 'toolbox_nyu_depth_v2/']);

scene_group = 'living_rooms_part4'; % kitchens
raw_dir = ['raw_scenes/', scene_group,'/']; % dir where the raw scenes are
scene_dirs = dir([raw_dir, 'living_room*']); % names for the video scenes
a=[scene_dirs.isdir];
scene_dirs = scene_dirs(a); % to get only the scene folders
step = 10; % sample step of the videos

for sc=1:size(scene_dirs,1)
    scene = scene_dirs(sc).name
    out_dir = [dataroot, 'background_scenes/', scene_group, '/', scene];
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    stmp = strsplit(scene, '_');
    %scene_id = stmp{end};
    
    list = get_synched_frames([raw_dir, scene]);

    for i=1:step:size(list,2)
        % filename to be saved later
        %filename = [scene_id, '_', num2str(i), '_']; % scene_img_
        filename = [scene, '_', num2str(i), '_']; % scene_img_
        fprintf('%s\n', filename);
        %if exist([out_dir, '/', filename, 'im.png'], 'file'), continue; end
        if exist([out_dir, '/', filename, 'depth.mat'], 'file'), continue; end
        
        dfile = list(i).rawDepthFilename;
        imfile = list(i).rawRgbFilename;
        try
            depth = imread([raw_dir, scene, '/', dfile]);
            im = imread([raw_dir, scene, '/', imfile]);
        catch
            continue; 
        end
        %figure; imagesc(im);
        %figure; imagesc(depth);
        depth = swapbytes(depth);
        [depthOut, rgbOut] = project_depth_map(depth, im);
        %figure; imagesc(rgbOut);
        %figure; imagesc(depthOut);
        
        % crop the images a bit
        [imh, imw, ~] = size(im);
        im = rgbOut(40:imh-10, 40:imw-40,:);
        depth = depthOut(40:imh-10, 40:imw-40,:);
        %figure; imagesc(im);
        %figure; imagesc(depth);

        % fill in the missing depth
        depth_filled = fill_depth_colorization(double(im), depth);
        %figure; imagesc(depth);
        
        
        imwrite(im, [out_dir, '/', filename, 'im.png']);
        save([out_dir, '/', filename, 'depth.mat'], 'depth_filled');
        %imwrite(depth, [scene, '/', filename, 'depth.pgm']);
    end


end
