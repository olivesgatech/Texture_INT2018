function S = dipSemb2D03(Dip, Data)
% DIPSEMB2D calculates the dip-guided semblance of 2D seismic data
%   S = dipSemb2D(Dip, Data) takes the 2D seismic array "Data" and its
%   corresponding dip array "Dip" as the input, "Data" and "Dip" have the
%   same size with the output "S". S is the dip-guided semblance based on
%   the algorithm from the book of Kurt J. Marfurt. In dipSemb2D, the
%   radius is set to 1 for better performance.

[rows, cols] = size(Data);        % Size of Data which is same with semblance
cropDip = zeros(rows,cols);   % Dip array will be cropped to the size of 1+2*r: rows-2*r, 1+r:cols-r
r = 1;                                     % Radius

% % model** is the array with the size of rows * cols * ((2*r+1).^2)
% % It stores the neighboring information of every pixel
% % Different dip will cause different neighboring
model0 = zeros(rows,cols,(2*r+1).^2);           % Degree 0
model90 = zeros(rows,cols,(2*r+1).^2);         % Degree 90
modelP45 = zeros(rows,cols,(2*r+1).^2);       % Degree +45
modelN45 = zeros(rows,cols,(2*r+1).^2);       % Degree -45

% Store the neighboring information into the array of "model**"
for i = 1:(2*r+1)
    for j = 1:(2*r+1)
        model0(1+2*r:rows-2*r,1+r:cols-r,3*(i-1)+j) = Data(i+1:rows-(4-i), j : cols+j-3);
        model90(1+2*r:rows-2*r,1+r:cols-r,3*(i-1)+j) = Data(1+j:rows-(4*r-j), i:cols-(2*r+1-i));
        modelP45(1+2*r:rows-2*r,1+r:cols-r,3*(i-1)+j) = Data((1*(i-1)+(4-j)) : (rows-(4*r+1-(1*(i-1)+(4-j)))), j:cols+j-3);
        modelN45(1+2*r:rows-2*r,1+r:cols-r,3*(i-1)+j) = Data((1*(i-1)+j):(rows-(4*r+1-(1*(i-1)+j))),j:cols+j-3);
    end
end

% Seperate the dips into four regions between -90 to 90
% Region1: 67.5<Dip<=90 or -90<Dip<=-67.5
% Region2: 22.5<Dip<=67.5
% Region3: -22.5<Dip<=22.5
% Region4: -67.5<Dip<=-22.5
cropDip(1+2*r:rows-2*r,1+r:cols-r) = Dip(1+2*r:rows-2*r,1+r:cols-r);
idx0 = find(cropDip>=-22.5 & cropDip<=22.5);
idx90 = find((cropDip>67.5 & cropDip<=90) | (cropDip>-90 & cropDip<=-67.5));
idxP45 = find(cropDip>22.5 & cropDip<=67.5);
idxN45 = find(cropDip<=-22.5 & cropDip>-67.5);

% Calculate the semblance based on the pixelwise dip information
% model** contains the dip guided neighboring information
% idx** contains the pixel indices in the corresponding model array
sem0 = semb03(model0, idx0, rows, cols, r);
sem90 = semb03(model90, idx90, rows, cols, r);
semP45 = semb03(modelP45, idxP45, rows, cols, r);
semN45 = semb03(modelN45,idxN45,rows,cols,r);

S = sem0 + sem90 + semP45 + semN45;

% Pixel-based Algorithm with Higher Time Cost
% Radius r = 1
% num = zeros(rows, cols);        % Numerator of semblance
% denom = zeros(rows, cols);    % Denominator of semblance
% for i = 1+2*r:rows-2*r              
%     for j = 1+r:cols-r
%         Classify Dip into four groups corresponding to different angel
%         intervals
%         [-22.5, 22.5]
%         if Dip(i,j)>=-22.5 && Dip(i,j)<=22.5
%             tmp = Data(i-r:i+r,j-r:j+r);
%         (22,5,67.5]
%         elseif Dip(i,j)>22.5 && Dip(i,j)<=67.5
%             tmp = [Data(i,j-1), Data(i-1,j), Data(i-2,j+1);
%                        Data(i+1,j-1), Data(i,j), Data(i-1,j+1);
%                        Data(i+2,j-1), Data(i+1,j), Data(i,j+1)];
%         (-22.5, -67.5]
%         elseif Dip(i,j)<-22.5 && Dip(i,j)>=-67.5
%             tmp = [Data(i-2,j-1), Data(i-1,j), Data(i,j+1);
%                        Data(i-1,j-1), Data(i,j), Data(i+1,j+1);
%                        Data(i,j-1), Data(i+1,j), Data(i+2,j+1)];
%         (67.5,90] & (-67.5,-90)
%         else
%             tmp0 = Data(i-r:i+r,j-r:j+r);
%             tmp = tmp0';
%         end
%         Calculate the dip-guided semblance with an approximated method
%         for k=1:2*r+1
%             tmprw = tmp(k,:);
%             num(i,j) = num(i,j) + (1/length(tmprw)*sum(real(tmprw(:)))).^2+(1/length(tmp)*sum(imag(tmprw(:)))).^2;
%             denom(i,j) = denom(i,j) + 1/length(tmprw)*(sum(real(tmprw(:)).^2)+sum(imag(tmprw(:)).^2));
%         end
%     end
% end
% 
% S = num./denom;

% Pixel-based Algorithm with Higher Time Cost
% radius r2 = 2
% Dip(i,j)>22.5 && Dip(i,j)<=67.5
%             tmp = [Data(i,j-2), Data(i-1,j-1), Data(i-2,j), Data(i-3,j+1), Data(i-4,j+2);
%                        Data(i+1,j-2), Data(i,j-1), Data(i-1,j), Data(i-2,j+1), Data(i-3,j+2);
%                        Data(i+2,j-2), Data(i+1,j-1), Data(i,j), Data(i-1,j+1), Data(i-2,j+2);
%                        Data(i+3,j-2), Data(i+2,j-1), Data(i+1,j), Data(i,j+1), Data(i-1,j+2);
%                        Data(i+4,j-2), Data(i+3,j-1), Data(i+2,j), Data(i+1,j+1), Data(i,j+2)];

% radius r2 = 2
% Dip(i,j)<-22.5 && Dip(i,j)>=-67.5
%             tmp = [Data(i-4,j-2), Data(i-3,j-1), Data(i-2,j), Data(i-1,j+1), Data(i,j+2);
%                        Data(i-3,j-2), Data(i-2,j-1), Data(i-1,j), Data(i,j+1), Data(i+1,j+2);
%                        Data(i-2,j-2), Data(i-1,j-1), Data(i,j), Data(i+1,j+1), Data(i+2,j+2);
%                        Data(i-1,j-2), Data(i,j-1), Data(i+1,j), Data(i+2,j+1), Data(i+3,j+2);
%                        Data(i,j-2), Data(i+1,j-1), Data(i+2,j), Data(i+3,j+1), Data(i+4,j+2)];
