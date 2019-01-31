function DV = ELBP_DV(img, mapping)
% mapping=getmapping(16,'riu2');
% Various combinations of radius and the number of sampling points
DV = elbp(img,2,16,mapping,'h');
end