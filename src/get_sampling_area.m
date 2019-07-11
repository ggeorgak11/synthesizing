

function sampling_area = get_sampling_area(im, depth, segmentation, seg_size, params, plane_params)

if params.sampling_mode==0
    % randomly put the object anywhere on the image
    sampling_area = ones(seg_size);
    
elseif params.sampling_mode==1
    % using only seg
    uinds=[];
    for l=1:length(params.usf_lbls)
        a=find(segmentation==params.usf_lbls(l));
        uinds = [uinds;a];
    end
    sampling_area = zeros(seg_size);
    sampling_area(uinds) = 1;
else
    % using both seg and surf
    % do the support surface estimation
    [support_surfaces, ~, ~, ~, ~] = get_support_surface(im, depth, plane_params, params);
    %close all;
    
    uinds=[]; usf_lbl_map=zeros(seg_size);
    for l=1:length(params.usf_lbls)
        a=find(segmentation==params.usf_lbls(l));
        uinds = [uinds;a];
        usf_lbl_map(a) = params.usf_lbls(l);
    end
    umap = zeros(seg_size);
    umap(uinds) = 1;
    
    % check whether the support surface has any overlap with the semantic labels that
    % correspond to surfaces (dont use props).
    % If it has, use their union to choose the area for the object.
    % Need to change the way the sampling is being done, figure out how
    % to put the objects on the plane
%     sampling_area = zeros(seg_size);
%     for sf=1:size(support_surfaces,3)
%         surf = support_surfaces(:,:,sf);
%         comb = surf+umap; % figure; imagesc(comb);
%         % check their overlap
%         if ~isempty(find(comb==2)), sampling_area(find(surf))=1; end
%     end

    % get the map for the floor label first
    floor_inds = find(segmentation==1);
    floor_map = zeros(seg_size);
    floor_map(floor_inds) = 1;
    
    sampling_area = zeros(seg_size);
    for sf=1:size(support_surfaces,3)
        surf = support_surfaces(:,:,sf);
        % if the support surface has overlap with the label floor then
        % ignore it.
        fs_map = surf+floor_map;
        if ~isempty(find(fs_map==2)), continue; end
        
        % otherwise check for an overlap with the useful labels
        comb = surf+umap; % figure; imagesc(comb);       
        % check their overlap
        if ~isempty(find(comb==2)), sampling_area(find(surf))=1; end
    end

    
    % save the visualizations
    if params.save_imgs
        % save the semantic segmentation
        f1=figure; set(gcf,'Visible', 'off');
        imagesc(segmentation); axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        saveas(f1, [params.save_path, params.img_id, '_seg.png']);
        
        % create and save the visualization for the support surfaces
        over_surf_im=im;
        for sf=1:size(support_surfaces,3)
            surf = support_surfaces(:,:,sf);
            a = find(surf);
            r=over_surf_im(:,:,1); g=over_surf_im(:,:,2); b=over_surf_im(:,:,3);
            r(a)=255; g(a)=0; b(a)=0;
            over_surf_im(:,:,1)=r; over_surf_im(:,:,2)=g; over_surf_im(:,:,3)=b;
        end
        imwrite(over_surf_im, [params.save_path, params.img_id, '_surf.png']);
        %figure; imagesc(over_surf_im);
        
        % create and save the useful labels visualization
        f2=figure; set(gcf,'Visible', 'off'); 
        imagesc(usf_lbl_map); axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        saveas(f2, [params.save_path, params.img_id, '_usf_lbls.png']); 
        
        % create and save the sampling area
        sample_area_im=im;
        a=find(sampling_area);
        r=sample_area_im(:,:,1); g=sample_area_im(:,:,2); b=sample_area_im(:,:,3);
        r(a)=0; g(a)=0; b(a)=255;
        sample_area_im(:,:,1)=r; sample_area_im(:,:,2)=g; sample_area_im(:,:,3)=b;
        %figure; imagesc(sample_area_im);
        imwrite(sample_area_im, [params.save_path, params.img_id, '_sample.png']);
        
        % save the depth image
        f3=figure; set(gcf,'Visible', 'off'); 
        imagesc(depth); axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        saveas(f3, [params.save_path, params.img_id, '_depth.png']); 
        
    end
    
end