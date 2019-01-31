function DV =LBP_DV(img)
mapping=getmapping(16,'riu2');
% Various combinations of radius and the number of sampling points
DV = lbp(img,2,16,mapping,'nh'); 
