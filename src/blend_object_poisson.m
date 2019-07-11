
% INPUT:
%   im: background scene image 
%   objimg: the object image as given initially by BigBird
%   mask: the object mask from BigBird
%   p: the extend percentage on the object crop
%   topLeft: top-left coordinate superimposition of object in im
%   mode: blending mode (0 is for rectangular patch, 1 is for masked patch)

function [Y, coords, mask_cropped] = blend_object_poisson(topLeft, im, objimg, mask, p, mode)

Lf = imGradFeature(double(im));
Gf = imGradFeature(double(objimg));
coords=[];

if mode==1 % blending with rectangular patch
    disp('-- Blending with rectangular patch...');
    % crop img,mask to be tight on the object
    [r,c] = find(mask>0);
    left = min(c); right = max(c);
    up = min(r); down = max(r);
    mask_width = right-left;
    mask_height = down-up;

    % crop the object
    width_ex = round(mask_width*p);
    height_ex = round(mask_height*p);
    w_marg = round(width_ex/2);
    h_marg = round(height_ex/2);
    % add the extra margin
    left = left-w_marg;
    right = right+w_marg;
    up = up-h_marg;
    down = down+h_marg;
    %objimg_ex = objimg(up:down, left:right,:); 
    %figure; imagesc(objimg_ex);

    pw = down-up+1; ph = right-left+1;
    imbottom=topLeft(1)+pw-1; imright=topLeft(2)+ph-1;
    if imbottom>size(im,1) || imright>size(im,2); Y=0; return; end;
    
    Lf(topLeft(1):imbottom, topLeft(2):imright, :, :) = Gf(up:down, left:right, :, :);

    X = Lf(:,:,:,1);
    param = buildModPoissonParam( size(Lf) );
    Y = modPoisson( Lf, param, 1E-8 );

    %figure; imagesc(uint8(X));
    %figure; imagesc(uint8(Y)); title(['blended with rectangular patch, p=', num2str(p)]);
    
else
    % try with using only the object inside the mask
    disp('-- Blending with masked patch...');
    [r,c] = find(mask>0);
    left = min(c); right = max(c);
    up = min(r); down = max(r);
    %objimg = objimg(up:down, left:right,:); 
    mask_cropped = mask(up:down, left:right);
    %seg = segment_object(objimg, mask_cropped);

    %topLeft = [523 665]; % topLeft coordinate in the scene img
    pw = size(mask_cropped,1); ph = size(mask_cropped,2);
    imDown=topLeft(1)+pw-1; imRight=topLeft(2)+ph-1;
    if imDown>size(im,1) || imRight>size(im,2); Y=0; return; end;

    inds_cropped = find(mask_cropped>0);
    non_inds = find(mask_cropped<=0);
    % replace the mask inds with the corresponding features
    for i=1:size(Lf,4)
        Lfd = Lf(:,:,:,i);
        Gfd = Gf(:,:,:,i);

        l_patch = Lfd(topLeft(1):imDown, topLeft(2):imRight, :);
        l_patch_r = l_patch(:,:,1); l_patch_g = l_patch(:,:,2); l_patch_b = l_patch(:,:,3);
        l_patch_r(inds_cropped)=0; l_patch_g(inds_cropped)=0; l_patch_b(inds_cropped)=0;
        l_patch(:,:,1)=l_patch_r; l_patch(:,:,2)=l_patch_g; l_patch(:,:,3)=l_patch_b;
        %figure; imagesc(uint8(l_patch))

        g_patch = Gfd(up:down, left:right, :);
        g_patch_r = g_patch(:,:,1); g_patch_g = g_patch(:,:,2); g_patch_b = g_patch(:,:,3);
        g_patch_r(non_inds)=0; g_patch_g(non_inds)=0; g_patch_b(non_inds)=0;
        g_patch(:,:,1)=g_patch_r; g_patch(:,:,2)=g_patch_g; g_patch(:,:,3)=g_patch_b;
        %figure; imagesc(uint8(g_patch));

        l_patch = l_patch + g_patch;
        Lfd(topLeft(1):imDown, topLeft(2):imRight, :) = l_patch;

        Lf(:,:,:,i) = Lfd;
    end

    %X = Lf(:,:,:,1);
    param = buildModPoissonParam( size(Lf) );
    Y = modPoisson( Lf, param, 1E-8 );

    % return the coordinates of the superimposition
    coords = [topLeft(1) topLeft(2) imDown imRight];
    
    %figure; imagesc(uint8(X));
    %figure; imagesc(uint8(Y)); title('blended with masked patch');
    
end
