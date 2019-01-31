%% ---------------------------------------------------------------------
% This code was used to generate the results in Figure 10 of the paper: 
%
% Z. Long, Y. Alaudah, M. Qureshi, Y. Hu, Z. Wang, M. Alfarraj, 
% G. AlRegib, A. Amin, M. Deriche, S. Al-Dharrab, and H. Di, 
% "A comparative study of texture attributes for characterizing 
% subsurface structures in seismic volumes," Interpretation, vol. 6, 
% no. 4, pp. T1055-T1066, Nov. 2018.
%
% If you use this code, please consider citing the paper. 
%------------------------------------------------------------------------
function [PA, MA, MIU, FWIU] = seismicValidate(classifiedImage, groundTruth, numClasses)

if nargin<3
    numClasses = 3;
end

% vectorize:
classifiedImage = classifiedImage(:);
groundTruth = uint8(groundTruth(:));

% calculate confusion matrix:
[C,~] = confusionmat(groundTruth, classifiedImage);
C = C(2:4,2:4); % remove 0 which represents non-labeled areas
t = sum(C,2); % # of pixels in each class

% pixel accuracy
PA = trace(C)/sum(t);

% mean accuracy:
MA = (1/numClasses)*sum(diag(C)./t);

% Mean IU (intersection over union)
MIU = (1/numClasses)*sum(diag(C)./(t+sum(C',2)-diag(C) ));

% Frequency Weighted IU
FWIU = 1/(sum(t)) * sum( (t.*diag(C)) ./ (t+sum(C',2)-diag(C) ));

end