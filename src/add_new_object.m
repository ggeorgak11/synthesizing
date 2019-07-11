
function [im, info, glb_mask, mask_list] = add_new_object(mask_list, info, Y, objects, glb_mask, tmp_mask, objparams, imp_coords)

if isempty(mask_list), idx=1; else idx=length(mask_list)+1; end
mask_list(idx).map = tmp_mask;                
glb_mask = glb_mask + tmp_mask;
im=Y; % to superimpose in the same image

% create the object annotation
if isempty(info), idx=1; else idx=length(info)+1; end
info(idx).category = objects{objparams.rClass};
info(idx).label = objparams.rClass;
info(idx).top = imp_coords(1);
info(idx).left = imp_coords(2);
info(idx).bottom = imp_coords(3);
info(idx).right = imp_coords(4);
info(idx).mask = tmp_mask;
info(idx).elevation_cam = objparams.rCam;
info(idx).azimuth = objparams.rIndex;