
% Finds all large surfaces in the image including the support surfaces
% Georgios Georgakis, Jana Kosecka, 2016

function [support_surfaces, surface_eq, large_planes, merged_plane_labels, plane_eq] = get_support_surface(RGB, depth, params, par)

%% Interpretation parameters

% Nominal focal length for the Kinect's RGB camera
focal_length            = 525; % 570.3

% min_points - minimum number of points in a planar segment
min_points              = 20; % 20

% ransac_trials - number of ransac trials to use
ransac_trials           = 80; %40

% inlier_threshold - threshold used to decide if a point is an inlier to a hypothesis
inlier_threshold        = 0.01; % 0.01

% outlier_ratio - fraction of leftover points that triggers a recursive refit
outlier_ratio           = 0.25; % 0.25

% inlier_ratio - ratio used in plane merging to decide whether to merge planes
inlier_ratio            = 0.90; % 0.90

% dotpthreshold_m - dot product threshold - used in merging to decide if planes are sufficiently similar
dotpthreshold_m         = 0.1; % default 0.10

% max_planes - maximum number of merged planes
max_planes              = 100; % 50

%% Disparity image

% W = 1000 ./ double(d);  % 1/depth where depth is in meters should be between 1/0.1 and 1/10
% W(d == 0) = NaN;

dz = depth/1e3;
W = 1./double(dz);     % 1/depth where depth is in meters should be between 1/0.1 and 1/10
W(dz == 0) = NaN;
W(dz > 7) = NaN;  % JK ignore far away data - not accurate

[nrows,ncols,~] = size(RGB);
SegmentationScriptV2;
    
[merged_plane_labels, ~] = mxFindPlanes (W, uint32(labels2), params.fx_rgb, min_points, ransac_trials, inlier_threshold, ...
        outlier_ratio, inlier_ratio, dotpthreshold_m, max_planes);

% JK remove small connected components
% for each plane label computed how many connected components
% it has and remove the small ones - they typically result
% from spurius plane fits and far away from real planes

% JK remove the planes  which have very small support
remove_small_planes=0;
if remove_small_planes
    np = max(merged_plane_labels(:));
    l = 0;
    for k=1:np
        ss = find(merged_plane_labels == k);
        if length(ss) < min_wall_points
            merged_plane_labels(ss) = 0;
        else
            l = l+1;
            merged_plane_labels(ss) = l;
        end
    end
end

nplanes = max(merged_plane_labels(:));

pcloud = depthToCloud(depth, params.fx_rgb, params.fy_rgb, params.center);
dx = pcloud(:,:,1);
dy = pcloud(:,:,2);
dz = pcloud(:,:,3);

plane_eq = [];
% This could be done a bit more efficiently by sorting the labeled points
% Excluding label 0 which denotes unknown areas
for label = 1:nplanes    
    t = find(merged_plane_labels == label);  
    dzc = dz(t);  
    dxc = dx(t);
    dyc = dy(t);
    [ normSpxl ] = fitPlaneAffine( dxc, dyc, dzc );

    plane_eq = [plane_eq; normSpxl];
end  

% detect the largest support surfaces in order to suppress detections 
% that are found on these surfaces
% Correct object hypos are usually not on large surfaces
for label = 1:nplanes
    pop(label) = length(find(merged_plane_labels == label));
end
[pop_sort lab_ind] = sort(pop, 'descend');
largest=find(pop_sort > params.nPixelsSupport); % large surfaces must have more than 27000 pixels in a 640x480 image
largest=lab_ind(1:length(largest));
large_planes = zeros(size(merged_plane_labels,1), size(merged_plane_labels,2));
for i=1:length(largest)
    large_label = largest(i);
    ids = merged_plane_labels==large_label;
    large_planes(ids) = 1;
end

% The supporting surface has to be sufficiently aligned with gravity
% direction and have a large number of pixels
gravity = [0 1 0];
for label = 1:nplanes
    n = plane_eq(label, 1:3);
    %n = plane_normal(merged_planes(label))';
    sim(label) = abs(dot(n, gravity));
end

[v,inds] = sort(sim, 'descend'); % sort the planes given their alignment
% filter out the very small ones
candid_labels=[];
while isempty(candid_labels)
    ctr=0;
    for i=1:length(inds)
        if length(find(merged_plane_labels == inds(i))) > params.plane_size_thresh && v(i) > params.gravity_align_thresh % decide on a thresh
            ctr=ctr+1;
            candid_labels(ctr) = inds(i);
        end
    end
    params.plane_size_thresh = params.plane_size_thresh - 1000;
    params.gravity_align_thresh = params.gravity_align_thresh - 0.05;
end


% keep the candid_labels all to be candidate support surfaces. Then,
% instead of refuting segments that are too far, keep those that are close
% enough to any of the supporting surfaces
% Refine the large_surfaces as well.

%f2=figure; if ~params.figs_visible, set(gcf,'Visible', 'off'); end 
cols = length(candid_labels);
for i=1:length(candid_labels)
    lbl_mask = zeros(size(merged_plane_labels,1), size(merged_plane_labels,2));
    lbl_mask(find(merged_plane_labels == candid_labels(i))) = 1;
    support_surfaces(:,:,i) = lbl_mask;
    %support_surfaces(:,:,i) = 1-lbl_mask;
    large_planes(find(lbl_mask>0)) = 1; 
    surface_eq(i,:) = plane_eq(candid_labels(i),:);
    
    %subplot(1,cols,i);
    %imagesc(support_surfaces(:,:,i));   
    %axis image;
    %title(['Label:', num2str(candid_labels(i)), ' Align:', num2str(sim(candid_labels(i)))]);
end

far = find(depth>params.far_thresh);
large_planes(far) = 1;

%figure; if ~params.figs_visible, set(gcf,'Visible', 'off'); end
%imagesc(large_planes); title('large surfaces');



% f2=figure; if ~figs_visible, set(gcf,'Visible', 'off'); end
% imagesc(surface); title('support surface');


if par.save_imgs
    %imwrite(double(merged_plane_labels), [par.save_path, par.img_id, '_planes.png']);
    f=figure; if ~params.figs_visible, set(gcf,'Visible', 'off'); end 
    imagesc(merged_plane_labels);
    axis off; set(gca,'visible','off');
    set(gca,'LooseInset',get(gca,'TightInset'));
    saveas(f, [par.save_path, par.img_id, '_planes.png']);
    
    %for i=1:size(support_surfaces,3)
    %    imwrite(support_surfaces(:,:,i), [par.save_path, par.img_id, '_surf_', num2str(i), '.png']);
    %end
    
end



