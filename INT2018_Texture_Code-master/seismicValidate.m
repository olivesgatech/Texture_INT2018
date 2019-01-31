function [PA, MA, MIU, FWIU] = seismicValidate(classifiedImage, groundTruth, numClasses)

if nargin<3
    numClasses = 4;
end

%vectorize:
classifiedImage = classifiedImage(:);
groundTruth = groundTruth(:);

% Calculate confusion matrix:
[C,~] = confusionmat(groundTruth, classifiedImage);
t = sum(C,2); % # of pixels in each class

% pixel accuracy
PA = trace(C)/sum(t);

% mean accuracy:
MA = (1/numClasses)*sum(diag(C)./t);

% Mean IU (intersection over union)
MIU = (1/numClasses)*sum(diag(C)./(t+ sum(C',2)-diag(C)));

% Frequency Weighted IU
FWIU = 1/(sum(t))   *   sum(  (t.*diag(C))   ./   ( t+ sum(C',2)-diag(C)    ));

end