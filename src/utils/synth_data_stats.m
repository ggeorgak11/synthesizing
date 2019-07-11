close all; clear all; clc;

% Collect data stats from the annotation
% for now count the occurences of the object classes

set_dir = 'NYU_synth_sets/NYU_synth4_bl/'; % assuming running from create_synth_set
scene_dirs = dir([set_dir, 'kitchen_*']);
a=[scene_dirs.isdir];
scene_dirs = scene_dirs(a); % to get only the scene folders 

class_freq = zeros(11,1); % keep class frequencies

for sc=1:size(scene_dirs,1)
    path = [set_dir, scene_dirs(sc).name, '/annotations/']
    files = dir([path, '*.mat']);
    
    for i=1:size(files,1)
        load([path, files(i).name]);
        for j=1:size(info,2), class_freq(info(j).label)=class_freq(info(j).label)+1; end        
    end
end

class_freq

%save([set_dir, 'class_occurences.mat'], 'class_freq');
