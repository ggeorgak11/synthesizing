

% init the params needed

% root path
params.dataroot = '/media/ggeorgak/Elements/synthesizing_project_georgakis/';
% background scenes path
params.back_dir = [params.dataroot, 'background_scenes/kitchens/'];
% objects path
params.objroot = [params.dataroot, 'object_data/objects/'];
% load the object list
params.object_list = textread([params.dataroot, 'object_data/our_dataset_objects'], '%s');

% save directories
params.save_dir_blended = [params.dataroot, 'Synthetic_Sets/synth_tmp_bl/'];
params.save_dir_imposed = [params.dataroot, 'Synthetic_Sets/synth_tmp_si/'];
params.imgs_dir = 'imgs/';
params.annot_dir = 'annotations/';

% Add dependencies
addpath([params.dataroot, 'ModifiedPoisson/']);
addpath(genpath([params.dataroot, 'support_surface/']));
addpath([params.dataroot, 'toolbox_nyu_depth_v2/']);
addpath('utils/');

% decide the sampling area to position the objects: 0 is random anywhere in
% image (random positioning), 1 is only with seg, 2 is seg+surf (selective positioning)
params.sampling_mode = 2; 
% define the useful labels from seg, given the sampling method
if params.sampling_mode==1, params.usf_lbls = [6, 11, 13, 28, 39]; % include also props
elseif params.sampling_mode==2, params.usf_lbls = [6, 11, 13];
else params.usf_lbls=[];
end
% define the overlap needed with the seg for a patch to be valid, given the sampling method
if params.sampling_mode==1, params.per = 0.6;
elseif params.sampling_mode==2, params.per = 0.2;
else params.per = 0;
end
% choose the scaling mode of the object. 0 is random (random scaling), 1 use the depth (selective scaling)  ...
% 2 use the depth but with small scaling pertubations
params.scaling_mode=1;
if params.scaling_mode==0, params.obj_scaling_list = 0.2:0.1:1;
else params.obj_scaling_list=[];
end
params.k = 4; % number of synthetic images created per image 
params.min_obj = 3; % minimum number of objects in an image
params.max_obj = 5; % maximum number of objects in an image
% set the threshold for maximum overlap between blended objects
params.max_bl_overlap = 0.4;
% set of azimuth poses to be sampled (the full set is 0:3:357)
params.index_seq = 0:3:357; % [0:3:90 270:3:357];
% Save example images from the procedure
params.save_imgs = 0;
params.save_path = [params.dataroot, 'Synthetic_Sets/synth_tmp_bl_examples/'];
if params.save_imgs
    if ~exist(params.save_path, 'dir'), mkdir(params.save_path); end
end
