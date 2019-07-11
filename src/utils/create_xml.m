


function create_xml(gt, path, name, sc, imsize)
imh = imsize(1); imw = imsize(2);
fp = fopen([path, name, '.xml'],'w');

%fprintf(fp, '<annotation>\n\t<folder>%s</folder>\n', ['scene_', num2str(sc)]);
fprintf(fp, '<annotation>\n\t<folder>%s</folder>\n', sc);
fprintf(fp, '\t<filename>%s</filename>\n', name);

fprintf(fp, '\t<source>\n\t\t<database>GMU Database</database>\n');
fprintf(fp, '\t\t<annotation>GMU Synthetic Scenes</annotation>\n');
fprintf(fp, '\t\t<image>kinectV2</image>\n');
fprintf(fp, '\t</source>\n');

fprintf(fp, '\t<size>\n');
fprintf(fp, '\t\t<width>%d</width>\n',imw);
fprintf(fp, '\t\t<height>%d</height>\n',imh);
fprintf(fp, '\t\t<depth>%d</depth>\n',3);
fprintf(fp, '\t</size>\n');
fprintf(fp, '\t<segmented>0</segmented>\n');

for i=1:size(gt,2)
    fprintf(fp, '\t<object>\n');
    fprintf(fp, '\t\t<name>%s</name>\n', gt(i).category);   
    fprintf(fp, '\t\t<pose>Unspecified</pose>\n');
    fprintf(fp, '\t\t<truncated>0</truncated>\n');
    fprintf(fp, '\t\t<difficult>0</difficult>\n');
    fprintf(fp, '\t\t<bndbox>\n');
    fprintf(fp, '\t\t\t<xmin>%d</xmin>\n', gt(i).left);
    fprintf(fp, '\t\t\t<ymin>%d</ymin>\n', gt(i).top);
    fprintf(fp, '\t\t\t<xmax>%d</xmax>\n', gt(i).right);
    fprintf(fp, '\t\t\t<ymax>%d</ymax>\n', gt(i).bottom);
    fprintf(fp, '\t\t</bndbox>\n');
    fprintf(fp, '\t</object>\n');
end

fprintf(fp, '</annotation>');
fclose(fp);
