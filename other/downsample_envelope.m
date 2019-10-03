function [envelopedData,newFs] = downsample_envelope(data_in,windowLength,time)

% Check if osl_movavg is in your MATLAB path 
if exist('osl_movavg') ~= 2
    error('Add ...osl/osl-core/ to your MATLAB path');
end

% if exist('ROInets.rows') ~= 2
%     error('Add ...osl/MEG-ROI-nets/ to your MATLAB path');
% end


% remove 1st sample bias
data_in(:,1) = [];
% Get time index and remove 1st sample bias
time(:,1)  = [];

% Calculate sample frequency
Fs           = 1.0 / median(diff(time));
% Get no. of voxels
nVoxels      = ROInets.rows(data_in);

% calculate window size
winsize      = round(windowLength * Fs);

% Use osl_movavg to use a sliding window approach with windowLength = x
envelopedData = [];

%loop over parcels
for iVoxel = 1:1:nVoxels
    envelopedData(iVoxel,:) = osl_movavg(data_in(iVoxel,:), ...
                                         time, ...
                                         winsize);
end

DsFactor = length(time) / size(envelopedData, 2);
newFs    = Fs / DsFactor; clear DsFactor;

% Use logarithmic?
%envelopedData = 2 .* log(real(envelopedData) + eps(min(abs(envelopedData(:))))); % prevent log(0).

% Remove NaNs
envelopedData = envelopedData(:,any(~isnan(envelopedData)));  % for columns
end