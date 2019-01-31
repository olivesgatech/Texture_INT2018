function aveDip = dipweight2D(Dip, Data, r)
% DIPWEIGHT2D Weights the dips based on energy
%   aveDip = dipweight2D(Dip, Data, radius) weights the array Dip based 
%   on the energy of elements in array Data with a neighboring of radius r.
%   Dip and Data should have the same size and r must be an integer
%   The Algorithm from Arthur E. Barnes et al, "weighted average seismic attributes"

[rows, cols] = size(Data);     % Size of Data
enve = abs(Data).^2;            % Energy of elements in Data
aveDip = zeros(rows,cols);  % Averaged dip with the same size of Data
weightedDip = zeros(rows-2,cols-2,(2*r+1)^2);
weight = zeros(rows-2,cols-2);
% for i = 1+r : rows-r
%     for j = 1+r : cols-r
%         tmp = Dip(i-r:i+r,j-r:j+r);
%         energy = enve(i-r:i+r,j-r:j+r);
%         aveDip(i,j) = 1/sum(energy(:))*(sum(tmp(:).*energy(:)));
%     end
% end
for i = 1:2*r+1
    for j =1:2*r+1
        weightedDip(:,:,3*(i-1)+j) = Dip(i:rows-(2*r+1-i),j:cols-(2*r+1-j)).*enve(i:rows-(2*r+1-i),j:cols-(2*r+1-j));
        weight = weight + enve(i:rows-(2*r+1-i),j:cols-(2*r+1-j));
    end
end
aveDip(1+r:rows-r,1+r:cols-r) = 1./weight.*sum(weightedDip, 3);