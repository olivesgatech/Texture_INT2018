%CLBP returns the complete local binary pattern image or LBP histogram of an image.
%  [CLBP_S,CLBP_M,CLBP_C] = CLBP(I,R,N,MAPPING,MODE) returns either a local binary pattern
%  coded image or the local binary pattern histogram of an intensity
%  image I. The CLBP codes are computed using N sampling points on a 
%  circle of radius R and using mapping table defined by MAPPING. 
%  See the getmapping function for different mappings and use 0 for
%  no mapping. Possible values for MODE are
%       'h' or 'hist'  to get a histogram of LBP codes
%       'nh'           to get a normalized histogram
%  Otherwise an CLBP code image is returned.

%  [CLBP_S,CLBP_M,CLBP_C] = CLBP(I,SP,MAPPING,MODE) computes the CLBP codes using n sampling
%  points defined in (n * 2) matrix SP. The sampling points should be
%  defined around the origin (coordinates (0,0)).
%
%  Examples
%  --------
%       I=imread('rice.png');
%       mapping=getmapping(8,'u2'); 
%       [CLBP_SH,CLBP_MH]=CLBP(I,1,8,mapping,'h'); %CLBP histogram in (8,1) neighborhood
%                                  %using uniform patterns

% function [CLBP_S,CLBP_M,CLBP_D,CLBP_C] = clbp_ldp(varargin) % image,radius,neighbors,mapping,mode)
function CLBP_SMDCjoint = cldp(varargin) % image,radius,neighbors,mapping,mode)

% Version 0.2
% Authors: Zhenhua Guo, Lei Zhang, and David Zhang

% The implementation is based on lbp code from MVG, Oulu University, Finland
% http://www.ee.oulu.fi/mvg/page/lbp_matlab


% Check number of input arguments.
error(nargchk(1,5,nargin));

image=varargin{1};
d_image=double(image);

if nargin==1
%     spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
    spoints(:,:,1)=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];

    neighbors=8;
    mapping=0;
    mode='h';
end

if (nargin == 2) && (length(varargin{2}) == 1)
    error('Input arguments');
end

if (nargin > 2) && (length(varargin{2}) == 1)
    radius=varargin{2};
    neighbors=varargin{3};
    
    spoints=zeros(neighbors,2);

    % Angle step.
    a = 2*pi/neighbors;
    
    for i = 1:neighbors
        for j=1:2
            % j=1 => radius; j=2 => radis-1;
            spoints(i,1,j) = -(radius-j+1)*sin((i-1)*a);
            spoints(i,2,j) = (radius-j+1)*cos((i-1)*a);
        end
    end
    
    if(nargin >= 4)
        mapping=varargin{4};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 5)
        mode=varargin{5};
    else
        mode='h';
    end
end

if (nargin > 1) && (length(varargin{2}) > 1)
    spoints=varargin{2};
    neighbors=size(spoints,1);
    
    if(nargin >= 3)
        mapping=varargin{3};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 4)
        mode=varargin{4};
    else
        mode='h';
    end   
end

% Determine the dimensions of the input image.
[ysize xsize] = size(image);



miny=min(spoints(:,1,1));
maxy=max(spoints(:,1,1));
minx=min(spoints(:,2,1));
maxx=max(spoints(:,2,1));

% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;

% Coordinates of origin (0,0) in the block
origy=1-floor(min(miny,0));
origx=1-floor(min(minx,0));

% Minimum allowed size for the input image depends
% on the radius of the used LBP operator.
if(xsize < bsizex || ysize < bsizey)
  error('Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dx = xsize - bsizex;
dy = ysize - bsizey;

% Fill the center pixel matrix C.
C = image(origy:origy+dy,origx:origx+dx);
d_C = double(C);

bins = 2^neighbors;

% Initialize the result matrix with zeros.
CLBP_S=zeros(dy+1,dx+1);
CLBP_M=zeros(dy+1,dx+1);
% Similar to local derivative pattern
CLBP_D=zeros(dy+1,dx+1);
CLBP_C=zeros(dy+1,dx+1);

%Compute the LBP code image

for i = 1:neighbors
  y = spoints(i,1,1)+origy;
  x = spoints(i,2,1)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  % Check if interpolation is needed.
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
    N = d_image(ry:ry+dy,rx:rx+dx);
    D{i} = N >= d_C;   
    Diff{i} = abs(N-d_C);    
    MeanDiff(i) = mean(mean(Diff{i}));
  else
    % Interpolation needed, use double type images 
    ty = y - fy;
    tx = x - fx;

    % Calculate the interpolation weights.
    w1 = (1 - tx) * (1 - ty);
    w2 =      tx  * (1 - ty);
    w3 = (1 - tx) *      ty ;
    w4 =      tx  *      ty ;
    % Compute interpolated pixel values
    N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
        w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    D{i} = N >= d_C;     
    Diff{i} = abs(N-d_C);    
    MeanDiff(i) = mean(mean(Diff{i}));
  end  
end

% Encoding for smaller circle
for i = 1:neighbors
  y = spoints(i,1,2)+origy;
  x = spoints(i,2,2)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  % Check if interpolation is needed.
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
    N = d_image(ry:ry+dy,rx:rx+dx);
    D_smallercircle{i} = N >= d_C;   

  else
    % Interpolation needed, use double type images 
    ty = y - fy;
    tx = x - fx;

    % Calculate the interpolation weights.
    w1 = (1 - tx) * (1 - ty);
    w2 =      tx  * (1 - ty);
    w3 = (1 - tx) *      ty ;
    w4 =      tx  *      ty ;
    % Compute interpolated pixel values
    N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
        w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    D_smallercircle{i} = N >= d_C;     

  end  
end


% Difference threshold for CLBP_M
DiffThreshold = mean(MeanDiff);
% compute CLBP_S and CLBP_M
for i=1:neighbors
  % Update the result matrix.
  v = 2^(i-1);
  CLBP_S = CLBP_S + v*D{i};
  CLBP_M = CLBP_M + v*(Diff{i}>=DiffThreshold);
  % Interscale LDP
  CLBP_D = CLBP_D + v* xor(D{i},D_smallercircle{i});
end

% CLBP_C
CLBP_C = d_C>=mean(d_image(:));

%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    sizarray = size(CLBP_S);
    CLBP_S = CLBP_S(:);
    CLBP_M = CLBP_M(:);
    CLBP_D = CLBP_D(:);
    CLBP_S = mapping.table(CLBP_S+1);
    CLBP_M = mapping.table(CLBP_M+1);
    CLBP_D = mapping.table(CLBP_D+1);
    CLBP_S = reshape(CLBP_S,sizarray);
    CLBP_M = reshape(CLBP_M,sizarray);
    CLBP_D = reshape(CLBP_D,sizarray);

    % % another implementation method
%     for i = 1:size(CLBP_S,1)
%         for j = 1:size(CLBP_S,2)
%             CLBP_S(i,j) = mapping.table(CLBP_S(i,j)+1);
%             CLBP_M(i,j) = mapping.table(CLBP_M(i,j)+1);
%         end
%     end
end


CLBP_SMDCjoint=zeros(dy+1,dx+1);
CLBP_SMDCjoint=bins^2*CLBP_S+bins*CLBP_M+CLBP_D+bins^3*CLBP_C;

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
%     CLBP_S=hist(CLBP_S(:),0:(bins-1));
%     CLBP_M=hist(CLBP_M(:),0:(bins-1));
%     CLBP_D=hist(CLBP_D(:),0:(bins-1));
    CLBP_SMDCjoint=hist(CLBP_SMDCjoint(:),0:(2*bins^3-1));


    if (strcmp(mode,'nh'))
%         CLBP_S=CLBP_S/sum(CLBP_S);
%         CLBP_M=CLBP_M/sum(CLBP_M);
%         CLBP_D=CLBP_D/sum(CLBP_D);
        CLBP_SMDCjoint=CLBP_SMDCjoint/sum(CLBP_SMDCjoint);
    end
else
%     %Otherwise return a matrix of unsigned integers
%     if ((bins-1)<=intmax('uint8'))
%         CLBP_S=uint8(CLBP_S);
%         CLBP_M=uint8(CLBP_M);
%         CLBP_D=uint8(CLBP_D);
%     elseif ((bins-1)<=intmax('uint16'))
%         CLBP_S=uint16(CLBP_S);
%         CLBP_M=uint16(CLBP_M);
%         CLBP_D=uint16(CLBP_D);
%     else
%         CLBP_S=uint32(CLBP_S);
%         CLBP_M=uint32(CLBP_M);
%         CLBP_D=uint32(CLBP_D);

%     end
end






