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
%% --------------------------------------------------------------------- %
%                       Clear, and set Parameters:                       %
% ---------------------------------------------------------------------- %

clc; close all; clear;

%%%%% DEFAULT PARAMETERS %%%%%%
% -------------- SVM -----------------
outlierFraction = 0.018; % SVM slack 
numClasses = 3;
costMatrix =  1 - eye(numClasses,numClasses);
CodingVar = 'onevsall';
% ------------ Defaults --------------
wsize = 49;
sigma = 25; % Sigma for the Gaussian kernel
% --------------- SLIC  --------------
regionSize = 15;
regularizer = 0.3;

%% --------------------------------------------------------------------- %
%                       Load data, and prepare                           %
% ---------------------------------------------------------------------- %
load('.\\facies\\IL530_HST.mat');
DataHST=patch;
clear patch;
load('.\\facies\\IL530_LST.mat');
DataLST=patch;
clear patch;
load('.\\facies\\IL530_TST.mat');
DataTST=patch;
clear patch;

% number of images in each class:
HSTx = size(DataHST,1);
LSTx = size(DataLST,1);
TSTx = size(DataTST,1);
numImages = HSTx + LSTx + TSTx;

% make the images into sets of horizontal vectors
VectorsHST = reshape(DataHST, HSTx, wsize*wsize); 
VectorsLST = reshape(DataLST, LSTx, wsize*wsize); 
VectorsTST = reshape(DataTST, TSTx, wsize*wsize); 

% combine data into one single matrix
x = [VectorsHST; VectorsLST; VectorsTST]; 

% Construct labels vector
y = cell(numImages,1);
for i = 1:numImages
    if i <= HSTx
        y{i} = 'HST';
    elseif i<= HSTx+LSTx
        y{i} = 'LST';
    elseif i<= HSTx+LSTx+TSTx
        y{i} = 'TST';
    end
end
y = y';

% randomly shift the images (rows)
rng(1);
rndIndices = randperm(numImages);
X = x(rndIndices,:);
Y = y(rndIndices);

clear y;
load .\\facies\\F3_350_550_by_10.mat; % loading y: 462x951x21 double; original data
Crossline=y;
clear y;

%% --------------------------------------------------------------------- %
%                         Apply Attributes to Data                       %
% ---------------------------------------------------------------------- %
descriptorNames = {'GLCM','MCLBP','ELBP','LRI'};
numDescriptors = length(descriptorNames);

% make new dir to save results
mkdir('results');
cd('results');

% first index: 2 = 2 labeled sections
resultsPA = zeros(2,numDescriptors);
resultsMA = zeros(2,numDescriptors);
resultsMIU = zeros(2,numDescriptors);
resultsFWIU = zeros(2,numDescriptors);

% Make Gaussian Mask:
[Rows,Cols] = ndgrid(1:wsize, 1:wsize);
center = (wsize)/2;
exponent = ((Rows-center).^2 + (Cols-center).^2)./(2*sigma^2);
gMask   = 1*(exp(-exponent));

% for MCLBP:
mapping1 = getmapping(8,'u2');
mapping2 = getmapping(16,'u2');
mapping3 = getmapping(24,'u2');

for descriptor = 1:numDescriptors
    switch descriptor
        case 1 % GLCM
            feature_matrix = zeros(numImages, length(GLCM_DV(reshape(squeeze(X(1,:)),[wsize wsize]))));
            parfor i = 1:numImages
                old_img = reshape(squeeze(X(i,:)),[wsize wsize]);
                masked_img = gMask.*old_img;
                new_img = GLCM_DV(masked_img);
                feature_matrix(i,:) = new_img;
            end
        case 2 % MCLBP
            feature_matrix = zeros(numImages, length(MCLBP_DV(reshape(squeeze(X(1,:)),[wsize wsize]), mapping1, mapping2, mapping3)));
            parfor i = 1:numImages
                old_img = reshape(squeeze(X(i,:)),[wsize wsize]);
                masked_img = gMask.*old_img;
                new_img = MCLBP_DV(masked_img, mapping1, mapping2, mapping3);
                feature_matrix(i,:) = new_img;
            end
        case 3 % ELBP
            feature_matrix = zeros(numImages, length(ELBP_DV(reshape(squeeze(X(1,:)),[wsize wsize]),mapping2)));
            parfor i = 1:numImages
                old_img = reshape(squeeze(X(i,:)),[wsize wsize]);
                masked_img = gMask.*old_img;
                new_img = ELBP_DV(masked_img,mapping2);
                feature_matrix(i,:) = new_img;
            end
        case 4 % LRI
            feature_matrix = zeros(numImages, length(LRI_DV(reshape(squeeze(X(1,:)),[wsize wsize]))));
            parfor i = 1:numImages
                old_img = reshape(squeeze(X(i,:)),[wsize wsize]);
                masked_img = gMask.*old_img;
                new_img = LRI_DV(masked_img);
                feature_matrix(i,:) = new_img;
            end
    end
    
    XX = feature_matrix;

    % --------------------------------------------------------------------- %
    %                           Setup SVM and Train                         %
    % --------------------------------------------------------------------- %
    
    rng(1);
    t=templateSVM('CacheSize','maximal','OutlierFraction', outlierFraction, ...
        'Standardize',true, 'KernelFunction','linear','Verbose',0, ...
        'prior', 'uniform');
    Mdl = fitcecoc(XX,Y, 'Learners',t,'Options',statset('UseParallel',1),...
        'Cost',costMatrix,'FitPosterior',1 ,...
        'ClassNames',{'HST', 'LST', 'TST'},...
        'prior', 'uniform', 'Coding', 'binarycomplete');
    
    % compute the in-sample (training) error:
    isLoss = resubLoss(Mdl);
    disp(['Training accuracy is: %' num2str(100*(1-isLoss))]);
    
    % --------------------------------------------------------------------- %
    %                      Label Seismic Section (2D)                       %
    % --------------------------------------------------------------------- %
    
    crosslineNumber = [3,15]; % inline 370 & 490
    
    for crossline_idx = 1:2
        if crossline_idx == 1
            load ..\\facies\\IL370_GT.mat; % load gt (ground truth)
            gt = gt(105:end, 1:end); % crop TOP 104 rows (noise)
        else
            load ..\\facies\\IL490_GT.mat;
            gt = gt(105:end, 1:end); % crop TOP 104 rows (noise)
        end
        
        img = squeeze(Crossline(:,:,crosslineNumber(crossline_idx)));
        img = img(105:end, 1:end); % crop TOP 104 rows (noise)
        
        buffer = (wsize-1)/2;
        
        img_ext = padarray(img,[0,buffer],'symmetric', 'both');
        img_ext = padarray(img_ext,[buffer,0],'replicate','both');
        [imHeight, imWidth] = size(img_ext);
        
        classifiedImage =  uint8(zeros(size(img)));
        confidenceMap =  zeros(size(img));
        posteriorMap = zeros([size(img),numClasses]);
        
        [gx, gy] = vl_grad(img);
        img3D = cat(3, img, gx, gy);
        img_slic = vl_slic(single(img3D), regionSize, regularizer);
        numSeg = single(max(img_slic(:)));
        
        for i = 0:numSeg
            region = zeros(size(img_slic));
            xy2D = find(img_slic == i);
            if isempty(xy2D)
                continue;
            end
            
            region(xy2D) = 1;
            regionProperties = regionprops(region);
            centroid = round(regionProperties.Centroid);
            centroid = centroid(2:-1:1);
            
            Window = img_ext(centroid(1):centroid(1)+2*buffer, centroid(2):centroid(2)+2*buffer);
            
            % extract features and evaluate:
            Window = gMask.*Window;

            switch descriptor
                case 1 % GLCM
                    windowDV = GLCM_DV(Window);
                case 2 % MCLBP
                    windowDV = MCLBP_DV(Window,mapping1, mapping2, mapping3);
                case 3 % ELBP
                    windowDV = ELBP_DV(Window,mapping2);
                case 4 % LRI
                    windowDV = LRI_DV(Window);
            end
            
            [x,y] = ind2sub(size(img_slic),xy2D);
            
            [label] = predict(Mdl,windowDV);
            
            if  strcmp(label,'HST')
                classifiedImage(xy2D)  = 1;
            elseif strcmp(label,'LST')
                classifiedImage(xy2D)  = 2;
            elseif strcmp(label,'TST')
                classifiedImage(xy2D)  = 3;
            end
        end
        
        % colorize:
        coloredImage =  uint8(zeros([size(img),3]));
        [height, width] = size(classifiedImage);
        for i = 1:height
            for j = 1:width
                if gt(i,j)==0 % not labeled areas
                    classifiedImage(i,j) = 0;
                    coloredImage(i,j,:) = [255,255,255]; % white
                    continue;
                end
                if classifiedImage(i,j) == 1 % HST
                    coloredImage(i,j,:) = [255,0,0]; % red
                elseif classifiedImage(i,j) == 2 % LST
                    coloredImage(i,j,:) = [0,255,0]; % green
                elseif classifiedImage(i,j) == 3  % TST
                    coloredImage(i,j,:) = [0,0,255]; % blue
                end
            end
        end
        
        hh = figure;
        imshow(img,[]);
        name = strcat('C',num2str(crosslineNumber(crossline_idx)*10+340),descriptorNames(descriptor));
        title(name);
        hold on;
        h = imagesc(coloredImage);
        set(h, 'AlphaData', 0.35 );
        % save figure as an image with appropriate title: 
        saveas(hh, strcat(name{1},'.png'))

        % close figure
        close;
        
        % compute numerical results: 
        [PA, MA, MIU, FWIU] = seismicValidate(classifiedImage, gt);
        
        % All the results are stored here: 
        resultsPA(crossline_idx, descriptor) = PA;
        resultsMA(crossline_idx, descriptor) = MA;
        resultsMIU(crossline_idx, descriptor) = MIU;
        resultsFWIU(crossline_idx, descriptor) = FWIU;
        
    end
end

disp('Finally Done');
cd ..
