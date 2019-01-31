function descriptor = LRI_DV(img)
% K chosen as 3, as it give good results with low complxity
K=3;
% Normalize the image
mi=min(img(:));
ma=max(img(:));
img=(img-mi)/(ma-mi);
% Compute LRI-A for Normalized Image
LRI_A_im = LRIA(img,K);
% Createte edge image for the input
img=edge(img,'canny');
% Calculate LRI-A for the edge image
LRI_A_im_E = LRIA(img,K);
% Combine the two LRI results (one for normalized, 2nd for edge images)
descriptor=[LRI_A_im LRI_A_im_E];
end
