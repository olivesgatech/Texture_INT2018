function S = dipSemb2D05(Dip, Data)
% DIPSEMB2D05 calculates the dip-guided semblance of 2D seismic data based
% on the 5*5 block.
%   S = dipSemb2D05(Dip, Data) takes the 2D seismic array "Data" and its
%   corresponding dip array "Dip" as the input, "Data" and "Dip" have the
%   same size with the output "S". S is the dip-guided semblance based on
%   the algorithm from the book of Kurt J. Marfurt. In dipSemb2D, the
%   radius is set to 2 for better performance.

[rows, cols] = size(Data);        % Size of Data which is same with semblance
cropDip = zeros(rows,cols);   % Dip array will be cropped to the size of 1+2*r: rows-2*r, 1+r:cols-r
r = 2;                                     % Radius

% model** is the array with the size of rows * cols * ((2*r+1).^2)
% It stores the neighboring information of every pixel
% Different dip will cause different neighboring
model0 = zeros(rows,cols,(2*r+1).^2);           % Degree 0
model90 = zeros(rows,cols,(2*r+1).^2);         % Degree 90
modelP45 = zeros(rows,cols,(2*r+1).^2);       % Degree +45
modelN45 = zeros(rows,cols,(2*r+1).^2);       % Degree -45

% Here for the degree of +/-22 and +/-68,  to have more accurate
% results, we choose the size of neighboring region as (2*r+1)*(2*r+2).
modelP22 = zeros(rows, cols, (2*r+1)*(2*r+2));      % Degree 22
modelP68 = zeros(rows, cols, (2*r+1)*(2*r+2));      % Degree 68
modelN22 = zeros(rows, cols, (2*r+1)*(2*r+2));      % Degree -22
modelN68 = zeros(rows, cols, (2*r+1)*(2*r+2));      % Degree -68


% Store the neighboring information into the array of "model**"
for i = 1:(2*r+1)
    for j = 1:(2*r+1)
        model0(1+2*r :rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data(i+r : rows-(7-i),1+ j : cols+j-6);
        model90(1+2*r :rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data(r+j:rows-7+j, i+1:cols-6+i);
        modelP45(1+2*r :rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data((1*(i-1)+(6-j)) : rows+i-j-4, j+1:cols+j-6);
        modelN45(1+2*r :rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data((1*(i-1)+j):rows-10+i+j,j+1:cols+j-6);
    end
end

% The size of neighboring region has been changed to (2*r+1)*(2*r+2)
for i = 1:(2*r+1)
    for j = 1:(2*r+2)
        if j<=r+1
%         modelP22(1+2*r:rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data(i+3: rows-6+i, 1+j : cols-6+j);  % ver1
    modelP22(1+2*r:rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data(i+2: rows-7+i, 1+j : cols-6+j);  % ver2
        modelN22(1+2*r:rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data(2+i: rows-7+i, 1+j: cols-6+j);
        modelP68(1+2*r:rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data(8-j:rows-j-1,i+1:cols-6+i);
        modelN68(1+2*r:rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data(j+2:rows-7+j, 1+i: cols-6+i);
        else
%         modelP22(1+2*r:rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data(i+2:rows-7+i, j: cols+j-7);  % ver1
    modelP22(1+2*r:rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j)=Data(i+1:rows-8+i, j: cols+j-7);    % ver2                     
        modelN22(1+2*r:rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data(3+i: rows-6+i, j:cols-7+j);
        modelP68(1+2*r:rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data(9-j:rows-j,i+2:cols-5+i);
        modelN68(1+2*r:rows-2*r,1+(r+1):cols-(r+1),(2*r+1)*(i-1)+j) = Data(j+1:rows-8+j, 2+i: cols-5+i);
        end
    end
end

% Seperate the dips into eight regions between -90 to 90
% Region1: ~=90:  -90<Dip<=-79 && 79<Dip<=90
% Region2: ~=-68: -79<Dip<=-56.5
% Region3: ~=-45: -56.5<Dip<=-33.5
% Region4: ~=-22: -33.5<Dip<=-11
% Region5: ~==0: -11<Dip<=11
% Region6: ~=22: 11<Dip<=33.5
% Region7: ~=45: 33.5<Dip<=56.5
% Region8: ~=68: 56.5<Dip<=79
cropDip(1+2*r:rows-2*r,2+r:cols-r-1) = Dip(1+2*r:rows-2*r,2+r:cols-r-1);
idx90 = find((cropDip>79 & cropDip<=90) | (cropDip>-90 & cropDip<=-79));
idxP68 = find(cropDip>56.5 & cropDip<=79);
idxP45 = find(cropDip>33.5 & cropDip<=56.5);
idxP22 = find(cropDip>11 & cropDip<=33.5);
idx0 = find(cropDip>-11 & cropDip<=11);
idxN22 = find(cropDip>-33.5 & cropDip<=-11);
idxN45 = find(cropDip>-56.5 & cropDip<=-33.5);
idxN68 = find(cropDip>-79 & cropDip<=-56.5);

% Calculate the semblance based on the pixelwise dip information
% model** contains the dip guided neighboring information
% idx** contains the pixel indices in the corresponding model array
sem0 = semb05(model0, idx0, rows, cols, r);
sem90 = semb05(model90, idx90, rows, cols, r);
semP45 = semb05(modelP45, idxP45, rows, cols, r);
semN45 = semb05(modelN45,idxN45,rows,cols,r);
semP22 = semb05(modelP22, idxP22, rows, cols, r);
semN22 = semb05(modelN22, idxN22, rows, cols, r);
semP68 = semb05(modelP68, idxP68, rows, cols, r);
semN68 = semb05(modelN68,idxN68,rows,cols,r);

S = sem0 + sem90 + semP45 + semN45+ semP22+semN22+semP68+semN68;



% % Pixel-based Algorithm with Higher Time Cost
% % Radius r = 1
% num = zeros(rows, cols);        % Numerator of semblance
% denom = zeros(rows, cols);    % Denominator of semblance
% for i = 1+2*r:rows-2*r              
%     for j = 1+r:cols-r
%         % Classify Dip into four groups corresponding to different angel
%         % intervals
%         % [-22.5, 22.5]
%         if Dip(i,j)>=-22.5 && Dip(i,j)<=22.5
%             tmp = Data(i-r:i+r,j-r:j+r);
%         % (22,5,67.5]
%         elseif Dip(i,j)>22.5 && Dip(i,j)<=67.5
%             tmp = [Data(i,j-1), Data(i-1,j), Data(i-2,j+1);
%                        Data(i+1,j-1), Data(i,j), Data(i-1,j+1);
%                        Data(i+2,j-1), Data(i+1,j), Data(i,j+1)];
%         % (-22.5, -67.5]
%         elseif Dip(i,j)<-22.5 && Dip(i,j)>=-67.5
%             tmp = [Data(i-2,j-1), Data(i-1,j), Data(i,j+1);
%                        Data(i-1,j-1), Data(i,j), Data(i+1,j+1);
%                        Data(i,j-1), Data(i+1,j), Data(i+2,j+1)];
%         % (67.5,90] & (-67.5,-90)
%         else
%             tmp0 = Data(i-r:i+r,j-r:j+r);
%             tmp = tmp0';
%         end
%         % Calculate the dip-guided semblance with an approximated method
%         for k=1:2*r+1
%             tmprw = tmp(k,:);
%             num(i,j) = num(i,j) + (1/length(tmprw)*sum(real(tmprw(:)))).^2+(1/length(tmp)*sum(imag(tmprw(:)))).^2;
%             denom(i,j) = denom(i,j) + 1/length(tmprw)*(sum(real(tmp(:)).^2)+sum(imag(tmprw(:)).^2));
%         end
%     end
% end
% 
% S = num./denom;


% Pixel-based Algorithm with Higher Time Cost
% This implementation is still not accurate because the radius increases from 1 to
% 2 without the increasing of the number of directions
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
