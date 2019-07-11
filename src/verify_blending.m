

function [add_obj, tmp_mask] = verify_blending(Y, imp_coords, mask_list, mask_cropped, seg_size, params)

tmp_mask = zeros(seg_size);
if size(Y,1)>1 % if blending was carried out                   
    tmp_mask(imp_coords(1):imp_coords(3), imp_coords(2):imp_coords(4)) = mask_cropped;

    % check whether the bounding box makes sense
    if imp_coords(1)>=imp_coords(3) || imp_coords(2)>=imp_coords(4)
        add_obj=0;
        disp('rejected due to wrong bounding box coordinates');
        return;
    end
    
    % we need to check whether the added object does not
    % completely overlap other objects in the synthetic image
    overlap=[]; add_obj=1;
    for m=1:size(mask_list,2)
        inter = intersect(find(mask_list(m).map), find(tmp_mask));
        %overlap(m) = length(inter) / ( length(find(mask_list(m).map))+length(find(tmp_mask)) - length(inter) ); % IoU
        overlap(m) = length(inter) / ( min(length(find(mask_list(m).map)),length(find(tmp_mask))) ); % inter/min(a1,a2)
        if overlap(m) > params.max_bl_overlap, add_obj=0; disp('rejected due to high overlap'); end
    end
else
    add_obj=0;
    disp('Blending not carried out!');
end