function D = dipesti2D (Data)
% DIPESTI2D Estimate the dips in 2D seismic sections
%   DIPESTI2D takes 2D seismic section "Data" with a
%   complex data format as its input and return a matrix of 
%   the same size as "Data", where the dips are stored. 
%   The default unit of dips is radian. 
%   (Add: a parameter for choosing unit rad or deg)

s = real(Data);         % s is real part of seismic data
sH = imag(Data);     % sH is imaginary part of seismic data

% Gradient calculation
[s_x, s_t] = gradient(s);
[sH_x,sH_t] = gradient(sH);

% Instantaneous wavenumber calculation
% kx on the crossline direction
% kz on the time direction
kx = (s.*sH_x-sH.*s_x)./(s.^2+sH.^2);
kz = (s.*sH_t-sH.*s_t)./(s.^2+sH.^2);
D = atan(kx./kz);


