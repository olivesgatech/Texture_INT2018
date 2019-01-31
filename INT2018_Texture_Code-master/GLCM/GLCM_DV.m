function DV = GLCM_DV(img)

% Normalising
mi=min(img(:));
ma=max(img(:));
img=(img-mi)/(ma-mi);

% GLCM Contrat, Homogeneity, Energy, Correlation
glcm=graycomatrix(img,'NumLevels',64,'GrayLimits',[],'offset',[1 0;1 1;0 1;2 0;2 2;2 -2;1 2;-2 0]);
stats = graycoprops(glcm,{'contrast','homogeneity','energy','correlation'});

% GLCM Mutual Information
glcm_sum=sum(glcm,3);
glcm_norm=glcm_sum./sum(glcm_sum(:));
Pi=sum(glcm_norm,2);
Pj=sum(glcm_norm,1);
Pij=Pi*Pj;
M_I=sum(glcm_norm(:).*(log2((glcm_norm(:)+eps)./(Pij(:)+eps))));


% GLCM Entropy
entropy1=(sum(sum(glcm(:,:,1).*log(glcm(:,:,1)+eps))));
entropy2=(sum(sum(glcm(:,:,2).*log(glcm(:,:,2)+eps))));
entropy3=(sum(sum(glcm(:,:,3).*log(glcm(:,:,3)+eps))));

DV=[stats.Contrast stats.Homogeneity stats.Energy stats.Correlation entropy1 entropy2 entropy3 M_I];
end