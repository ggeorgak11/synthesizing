%
% SegmentationScriptV2 : this script assumes that the variable RGB contains
% an image and nrows and ncols define its size.
%

%% Segmentation parameters

% Edge detection parameters
edge_thresh     = 0.04;
edge_sigma      = 2;

% sigma_b : parameter used to blur image prior to sampling color values
sigma_b         = 3;

% sigma_c : parameter used to compute color weights
sigma_c         = 3;

% merge_threshold : parameter that defines the final segmentation - larger
% values lead to under segmentation, smaller values lead to over
% segmentation.
merge_threshold = 0.4;


%% Edge detection

e = edge(rgb2gray(RGB), 'canny', edge_thresh, edge_sigma);

% imtool (e);

%% Extract edgel coordinates
[Y, X] = find (e);

%% Compute Delaunay Triangulation
% Note that the edgel coordinates are integral so they won't be in general
% position this means that the delaunay routines must implicitly 'joggle'
% the coordinates by adding random noise.

% Add in the four corners for completeness
corners = [1,1; ncols, 1; 1, nrows; ncols, nrows];

X = [corners(:,1); X];
Y = [corners(:,2); Y];


%% Compute Delaunay Triangulation

dt = DelaunayTri (X, Y);

%% Find the neighbors of each of the triangles in the triangulation

SN = dt.neighbors();

ntri = size(SN,1);

%% Compute circumcircle centers and radii

[CC, RCC] = dt.circumcenters();

%% Use rasterization to assign every pixel the label of its enclosing triangle
% There will be inevitable cases of triangle fighting here but it shouldn't
% cause any issues. Note that we use 0 for the bg_color so if a pixel is
% not in any triangle it will cause problems later.

labels = mxRasterizeTriangles (nrows, ncols, dt.X, int32(dt.Triangulation), 0);

%% Compute overlap ratios

[overlap_ratios, overlap_areas, edge_lengths] = ComputeCircumcircleOverlap (SN, CC, RCC);

%% Implement blurring

width = 2*ceil(sigma_b) + 1;

x = -width:width;

gauss  = exp (-(x/sigma_b).^2);

% normalize to unit sum
gauss = gauss / sum(gauss);

RGBb = imfilter (double(RGB),   gauss, 'replicate');
RGBb = imfilter (RGBb, gauss', 'replicate');

%% Compute Circumcircle statistics

[mean_color_rgb, pixel_count] = mxCircleStats (CC, RCC, RGBb);

%% Convert mean colors from RGB to Lab so that differences will hopefully be more perceptual
% Note that we could convert Ib to Lab prior to finding the mean color
% which would be another approach. You would need to be careful about lab
% values as uint8 vs doubles because of scaling issues

cform = makecform('srgb2lab');

% Note that we divide by 255 to scale the colors to the range 0-1
mean_color_lab = applycform (mean_color_rgb/255.0, cform);

% Note that the mean_color_lab should be double precision and distances in
% this space should be perceptually meaningful.


%% Compute color weights

dc2 = compute_squared_edge_distances (SN, mean_color_lab);
color_weights = exp(-(1/sigma_c^2) * dc2);

%% Merge Triangles based on edge costs.

mergeTimeId = tic;
[old_labels, new_labels, merge_costs, merge_lengths] = mxMergeTriangles (int32(SN), (1-color_weights), edge_lengths);
mergeTime = toc(mergeTimeId);

%% Threshold links and merge

t = merge_costs < merge_threshold;

% Integrate Merges
tri_labels = integrate_merges(ntri, old_labels(t), new_labels(t));

% Translate the triangle labels into pixel labels
labels2 = compactLabels(tri_labels(labels));

% segmentTime = toc;

%% Play with colorImageSegments

I2 = colorImageSegments (RGB, labels2, [1 0 0]);

%% Display image segmentation results
% imtool (I2);