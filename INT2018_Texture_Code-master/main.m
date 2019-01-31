%% ---------------------------------------------------------------------
% This code was used to generate the results in the following paper: 
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
numClasses = 4;
costMatrix =  1 - eye(numClasses,numClasses);
CodingVar = 'onevsall';
% ------------ Defaults --------------
windowSize = 99;
sigma = 25; % Sigma for the Gaussian kernel
% --------------- SLIC  --------------
regionSize = 15; % used to be 35
regularizer = 0.3; % or 0.3 or 0.05

%% --------------------------------------------------------------------- %
%                       Load data, and prepare                           %
% ---------------------------------------------------------------------- %
DataCH = load(['data',filesep,'Patches_ch.mat']);
DataSM = load(['data',filesep,'Patches_sm.mat']);
DataHO = load(['data',filesep,'Patches_ho.mat']);
DataFA = load(['data',filesep,'Patches_fa.mat']);
DataS1 = load(['data',filesep,'Patches_sa1.mat']);
DataS2 = load(['data',filesep,'Patches_sa2.mat']);
DataS3 = load(['data',filesep,'Patches_sa3.mat']);

DataCH = DataCH.Patches_ch; % chaotic
DataSM = DataSM.Patches_sm; % smooth
DataHO = DataHO.Patches_ho; % horizon
DataOT = cat(1, DataSM, DataHO); %other class
DataFA = DataFA.Patches_fa; % faults
DataS1 = DataS1.Patches_sa1;
DataS2 = DataS2.Patches_sa2;
DataS3 = DataS3.Patches_sa3;
DataSA = cat(1, DataS1, DataS2, DataS3); % salt domes

% number of images in each class:
chx = size(DataCH,1);
otx = size(DataOT,1);
fax = size(DataFA,1);
sax = size(DataSA,1);

numImages = chx + otx + fax + sax;

% make the images into sets of horizontal vectors
VectorsCH = reshape(DataCH, chx, 99*99);
VectorsOT = reshape(DataOT, otx, 99*99);
VectorsFA = reshape(DataFA, fax, 99*99);
VectorsSA = reshape(DataSA, sax, 99*99);

% combine data into one single matrix
x = [VectorsCH; VectorsOT; VectorsFA; VectorsSA];

% Construct labels vector
y = cell(numImages,1);
for i = 1:numImages
    if i <= chx
        y{i} = 'Chaotic';
    elseif i<= chx+otx
        y{i} = 'Other';
    elseif i<= chx+otx+fax
        y{i} = 'Fault';
    elseif i<= chx+otx+fax+sax
        y{i} = 'Salt';
    end
end
y = y';

% randomly shift the images (rows)
rng(1);
rndIndices = randperm(numImages);
X = x(rndIndices,:);
Y = y(rndIndices);

load(['data',filesep,'crossline']);

% load Ground Truth:
load(['data',filesep,'GT61.mat' ]);
load(['data',filesep,'GT211.mat']);
load(['data',filesep,'GT231.mat']);
load(['data',filesep,'GT281.mat']);
%% --------------------------------------------------------------------- %
%                         Apply Attributes to Data                       %
% ---------------------------------------------------------------------- %
descriptorNames = {'GLCM','Semblance','LBP','CLBP','MCLBP','ELBP','CLDP','LRI'};
numDescriptors = length(descriptorNames);

% make new dir to save results
mkdir('results');
cd('results');

% first index: 4 = 4 labeled sections
resultsPA = zeros(4,numDescriptors);
resultsMA = zeros(4,numDescriptors);
resultsMIU = zeros(4,numDescriptors);
resultsFWIU = zeros(4,numDescriptors);

% Make Gaussian Mask:
[Rows,Cols] = ndgrid(1:windowSize, 1:windowSize);
center = (windowSize)/2;
exponent = ((Rows-center).^2 + (Cols-center).^2)./(2*sigma^2);
gMask   = 1*(exp(-exponent));
% imagesc(gMask) surf(gMask)

% for MCLBP:
mapping1 = getmapping(8,'riu2');
mapping2 = getmapping(16,'riu2');
mapping3 = getmapping(24,'riu2');

for descriptor = 1:numDescriptors
    switch descriptor
        case 1 % GLCM
            feature_matrix = zeros(numImages, length(GLCM_DV(reshape(squeeze(X(1,:)),[99 99]))));
            parfor i = 1:numImages
                old_img = reshape(squeeze(X(i,:)),[99 99]);
                masked_img = gMask.*old_img;
                new_img = GLCM_DV(masked_img);
                feature_matrix(i,:) = new_img;
            end
        case 2 % Semblance
            feature_matrix = zeros(numImages, length(SEM_DV(reshape(squeeze(X(1,:)),[99 99]))));
            parfor i = 1:numImages
                old_img = reshape(squeeze(X(i,:)),[99 99]);
                masked_img = gMask.*old_img;
                new_img = SEM_DV(masked_img);
                feature_matrix(i,:) = new_img;
            end
        case 3 % LBP
            feature_matrix = zeros(numImages, length(LBP_DV(reshape(squeeze(X(1,:)),[99 99]))));
            parfor i = 1:numImages
                old_img = reshape(squeeze(X(i,:)),[99 99]);
                masked_img = gMask.*old_img;
                new_img = LBP_DV(masked_img);
                feature_matrix(i,:) = new_img;
            end
        case 4 % CLBP
            feature_matrix = zeros(numImages, length(CLBP_DV(reshape(squeeze(X(1,:)),[99 99]),20)));
            parfor i = 1:numImages
                old_img = reshape(squeeze(X(i,:)),[99 99]);
                masked_img = gMask.*old_img;
                new_img = CLBP_DV(masked_img, 20);
                feature_matrix(i,:) = new_img;
            end
        case 5 % MCLBP
            feature_matrix = zeros(numImages, length(MCLBP_DV(reshape(squeeze(X(1,:)),[99 99]), mapping1, mapping2, mapping3)));
            parfor i = 1:numImages
                old_img = reshape(squeeze(X(i,:)),[99 99]);
                masked_img = gMask.*old_img;
                new_img = MCLBP_DV(masked_img, mapping1, mapping2, mapping3);
                feature_matrix(i,:) = new_img;
            end
        case 6 % ELBP
            feature_matrix = zeros(numImages, length(ELBP_DV(reshape(squeeze(X(1,:)),[99 99]),mapping2)));
            parfor i = 1:numImages
                old_img = reshape(squeeze(X(i,:)),[99 99]);
                masked_img = gMask.*old_img;
                new_img = ELBP_DV(masked_img,mapping2);
                feature_matrix(i,:) = new_img;
            end
        case 7 % CLDP
            feature_matrix = zeros(numImages, length(CLDP_DV(reshape(squeeze(X(1,:)),[99 99]),mapping2)));
            parfor i = 1:numImages
                old_img = reshape(squeeze(X(i,:)),[99 99]);
                masked_img = gMask.*old_img;
                new_img = CLDP_DV(masked_img,mapping2);
                feature_matrix(i,:) = new_img;
            end
        case 8 % LRI
            feature_matrix = zeros(numImages, length(LRI_DV(reshape(squeeze(X(1,:)),[99 99]))));
            parfor i = 1:numImages
                old_img = reshape(squeeze(X(i,:)),[99 99]);
                masked_img = gMask.*old_img;
                new_img = LRI_DV(masked_img);
                feature_matrix(i,:) = new_img;
            end
    end
    
    XX = feature_matrix;

    % --------------------------------------------------------------------- %
    %                           Setup SVM and Train                          %
    % ---------------------------------------------------------------------- %
    
    rng(1);
    t=templateSVM('CacheSize','maximal','OutlierFraction', outlierFraction, ...
        'Standardize',true, 'KernelFunction','linear','Verbose',0, ...
        'prior', 'uniform');
    Mdl = fitcecoc(XX,Y, 'Learners',t,'Options',statset('UseParallel',1),...
        'Cost',costMatrix,'FitPosterior',1 ,...
        'ClassNames',{'Chaotic', 'Other', 'Fault', 'Salt'},...
        'prior', 'uniform', 'Coding', 'binarycomplete');
    
    % CodingMat = Mdl.CodingMatrix;
    % Mdl.BinaryLearners{1}                % The first binary learner
    % Mdl.BinaryLearners{1}.SupportVectors % Support vector indices
    
    % compute the in-sample (training) error:
    isLoss = resubLoss(Mdl);
    disp(['Training error is: %' num2str(100*isLoss)]);
    
    % --------------------------------------------------------------------- %
    %                      Label Seismic Section (2D)                        %
    % ---------------------------------------------------------------------- %
    
    crosslineNumber = [61,211,231,281];
    
    for crossline_idx = 1:4
        
        img = squeeze(Crossline(crossline_idx,:,:));
        img = img(105:end, 1:end); % crop TOP 104 rows (noise)
        
        buffer = (windowSize-1)/2;
        
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
                case 2 % Semblance
                    windowDV = SEM_DV(Window);
                case 3 % LBP
                    windowDV = LBP_DV(Window);
                case 4 % CLBP
                    windowDV = CLBP_DV(Window, 20);
                case 5 % MCLBP
                    windowDV = MCLBP_DV(Window,mapping1, mapping2, mapping3);
                case 6 % ELBP
                    windowDV = ELBP_DV(Window,mapping2);
                case 7 % CLDP
                    windowDV = CLDP_DV(Window,mapping2);
                case 8 % LRI
                    windowDV = LRI_DV(Window);
            end
            
            [x,y] = ind2sub(size(img_slic),xy2D);
            
            [label] = predict(Mdl,windowDV);
            
            if  strcmp(label,'Chaotic')
                classifiedImage(xy2D)  = 2;
            elseif strcmp(label,'Other')
                classifiedImage(xy2D)  = 1; % these numbers correspond to those in ground truth. Don't change
            elseif strcmp(label,'Fault')
                classifiedImage(xy2D)  = 3;
            elseif strcmp(label,'Salt')
                classifiedImage(xy2D)  = 4;
            end
        end
        
        % colorize:
        coloredImage =  uint8(zeros([size(img),3]));
        [height, width] = size(classifiedImage);
        for i = 1:height
            for j = 1:width
                if classifiedImage(i,j) == 2 % CHaotic
                    coloredImage(i,j,:) = [0,0,255]; % blue
                elseif classifiedImage(i,j) == 1 % other:
                    coloredImage(i,j,:) = [184,184,184]; % Gray
                elseif classifiedImage(i,j) == 3  % Fault
                    coloredImage(i,j,:) = [0,255,0]; % Green
                elseif classifiedImage(i,j) == 4 % salt
                    coloredImage(i,j,:) = [255,0,0]; % red
                end
            end
        end
        
        hh = figure;
        imshow(img,[]);
        name = strcat('C',num2str(crosslineNumber(crossline_idx)),descriptorNames(descriptor));
        title(name);
        hold on;
        h = imagesc(coloredImage);
        set(h, 'AlphaData', 0.45 );
        % save figure as an image with appropriate title: 
        saveas(hh, strcat(name{1},'.png'))
        saveas(hh, strcat(name{1},'.bmp'))
        saveas(hh, strcat(name{1},'.fig'))
        saveas(hh, strcat(name{1},'.mfig'))

        % close figure
        close;
        
        switch crossline_idx
            case 1
                groundTruth = GT61;
            case 2
                groundTruth = GT211;
            case 3
                groundTruth = GT231;
            case 4
                groundTruth = GT281;
        end
        
        % compute results and print: 
        [PA, MA, MIU, FWIU] = seismicValidate(classifiedImage, groundTruth)
        
        % All the results are stored here: 
        resultsPA(crossline_idx, descriptor) = PA;
        resultsMA(crossline_idx, descriptor) = MA;
        resultsMIU(crossline_idx, descriptor) = MIU;
        resultsFWIU(crossline_idx, descriptor) = FWIU;
        
    end
end

disp('Finally Done');
cd ..




