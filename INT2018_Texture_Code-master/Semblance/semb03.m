function Sem = semb03(Model, Idx, Rows, Cols, R)
% SEMB calculates semblance of every pixel based on its model and index.
% Sem = semb(Model, Idx, Rows, Cols, R) takes the 3D array "Model" and 2D
% array "Idx" as its inputs. "Model" contains the (2*R+1) neighboring information
% of every pixel with direction considered. The directions, like 0, +45,
% -45, 90, means the direction of faults. "Idx" here actually is the
% classification results based on direction regions. Variables "Rows",
% "Cols" are used to set the size of Sem. 2D array "Sem" is the output and
% represents the dicontinuity on the perpendicular direction of fault
% lines. 
Num = zeros(Rows, Cols);
Denom = zeros(Rows,Cols);
Sem = zeros(Rows,Cols);

for i = 1:2*R+1
    tmp1 = Model(:,:,3*(i-1)+1);
    tmp2 = Model(:,:,3*(i-1)+2);
    tmp3 = Model(:,:,3*(i-1)+3);
    % Semblance calculation algorithm
    % Page 51, in the book "Seismic Attributes for Prospect Identification and Reservoir Characterization"
    Num = Num + (1/(2*R+1)*real(tmp1+tmp2+tmp3)).^2+(1/(2*R+1)*imag(tmp1+tmp2+tmp3)).^2;
    Denom = Denom + (1/(2*R+1))*(real(tmp1).^2+imag(tmp1).^2+real(tmp2).^2+imag(tmp2).^2+real(tmp3).^2+imag(tmp3).^2);
end

Sem(Idx) = Num(Idx)./Denom(Idx);