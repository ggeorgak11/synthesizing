

function params = get_plane_params

camera_params;
params.fx_rgb = fx_rgb;
params.fy_rgb = fy_rgb;
params.center = [cx_rgb cy_rgb];

params.far_thresh = 4000;
params.nPixelsSupport = 27000; % 30000
params.gravity_align_thresh = 0.7;
params.plane_size_thresh = 4000;

params.figs_visible = 0;