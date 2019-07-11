

function [im, im2, info] = blend_with_seg(im, depth, n, seg_size, umap, params)

%refined_masks_flag=1;

% im is the image with the blended objects and im2 is the image with the
% simple superimposed objects
im2 = im;
mask_list=[];
glb_mask = zeros(seg_size); % keep track of the imposed objects
info=[]; % annotation struct

for i=1:n
    [objimg_to_blend, mask_to_blend, top, left, objparams] = get_objimg_blend_surf(seg_size, depth, umap, params, i);
    
    if isempty(objimg_to_blend) || top==0, continue; end % no object chosen    
    
    [Y, imp_coords, mask_cropped] = blend_object_poisson([top left], im, objimg_to_blend, mask_to_blend, 0.1, 2); % mode 2 is masked blending
    [Y2, ~, ~] = superimpose_simple(objimg_to_blend, mask_to_blend, [top left], im2);
    
    % no need to do the verification for both
    [add_obj, tmp_mask] = verify_blending(Y, imp_coords, mask_list, mask_cropped, seg_size, params);
    
    if add_obj
        [im, info, glb_mask, mask_list] = add_new_object(mask_list, info, Y, params.object_list, glb_mask, tmp_mask, objparams, imp_coords);
        % info, glb_mask, mask_list are the same for both blended and
        % superimposed. No need to update them
        im2 = Y2;
    end                
    

end

