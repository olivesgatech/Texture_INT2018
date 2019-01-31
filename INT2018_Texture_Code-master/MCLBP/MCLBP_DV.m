function DV = MCLBP_DV(img, mapping1, mapping2, mapping3)
CLBP_H11 = mclbp(img,1,8,mapping1,'h');
CLBP_H12 = mclbp(img,2,16,mapping2,'h');
CLBP_H13 = mclbp(img,3,24,mapping3,'h');

DV = [CLBP_H11, CLBP_H12,CLBP_H13];
end