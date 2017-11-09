%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script uses the MEG ROINets 
%
% ROInet: A Matlab package for performing leakage-robust network inference 
% between ROIs in MEG data

% The methodology used in this pipeline is set out in

% Colclough, G. L., Brookes, M., Smith, S. M. and Woolrich, M. W., 
% "A symmetric multivariate leakage correction for MEG connectomes," 
% NeuroImage 117, pp. 439-448 (2015)
%
% We are grateful for Colclough and colleagues at OHBA for making their
% code openly available for re-use. Please check copyright information in 
% OSL for more information.
%
% These scripts use the MEG ROINets code, but performs source analysis in
% Fieldtrip.
%
% Users will need to download OSL (https://github.com/OHBA-analysis).
% However please don't initiate OSL - this script will add the relevent
% bits of code for you.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Specify locations of data and directory to store data

subject = '2704'

rawfile = '/Users/44737483/Documents/mcq_data/2704/meg/run-rs/2704_AT_ME160_2017_09_27_rs.con';
data_mrk = '/Users/44737483/Documents/mcq_data/2704/meg/run-rs/2704_AT_ME160_2017_09_27_rs_PRE.mrk';
dir_name = '/Users/44737483/Documents/scripts_mcq/resting_state/2704'; 
mkdir(dir_name); cd(dir_name);

cfg = [];
cfg.headerfile = rawfile; 
cfg.datafile = rawfile;
cfg.trialdef.triallength = Inf;
cfg.trialdef.ntrials = 1;
cfg = ft_definetrial(cfg)

cfg.continuous = 'yes';
alldata = ft_preprocessing(cfg);

% Resample to speed up data analysis - check whether 200Hz is optimal?
cfg = [];
cfg.resamplefs = 200;
cfg.detrend = 'no';
alldata = ft_resampledata(cfg,alldata);

% Create layout file for later + save
cfg             = [];
cfg.grad        = alldata.grad;
lay             = ft_prepare_layout(cfg, alldata);
save lay lay

% Define 300s of rs-data from the trigger
cfg = [];
cfg.dataset                 = rawfile;
cfg.continuous              = 'yes';
cfg.trialdef.prestim        = 0;         % pre-stimulus interval  0s
cfg.trialdef.poststim       = 300;         % post-stimulus interval 300s (5min)
cfg.trialfun                = 'mytrialfun_MMN';
data_raw                    = ft_definetrial(cfg);

% Redefine Trial
cfg = [];
data = ft_redefinetrial(data_raw,alldata); %redefines the filtered data

% Band-pass filter between 0.5-250Hz
cfg = [];
cfg.continuous = 'yes';
cfg.bpfilter = 'yes';
cfg.bpfreq = [1 40];
data = ft_preprocessing(cfg,data);

% %Deal with 50Hz line noise
% cfg = [];
% cfg.bsfilter = 'yes';
% cfg.bsfreq = [49 51];
% data = ft_preprocessing(cfg,data);
% 
% %Deal with 100Hz line noise
% cfg = [];
% cfg.bsfilter = 'yes';
% cfg.bsfreq = [99 101];
% data = ft_preprocessing(cfg,data);

% Cut out spare channels
cfg = [];
cfg.channel = data.label(1:160);
data = ft_selectdata(cfg,data);

%% COREG

% Addpath to coreg scripts
addpath(genpath('/Users/44737483/Documents/scripts_mcq/coreg'));

% Specify location of data
mri_file = '/Users/44737483/Documents/mcq_data/2704/anat/2704.nii';
hspfile = '/Users/44737483/Documents/mcq_data/2704/meg/2704_AT_ME160_2017_09_27.hsp';
elpfile = '/Users/44737483/Documents/mcq_data/2704/meg/2704_AT_ME160_2017_09_27.elp'

% Do coreg - this takes 3-5 minutes
coreg_yokogawa_icp_adjust_weights(dir_name,rawfile,data_mrk,mri_file,hspfile,elpfile,100,0.04,'yes',10)

% Load the coregistered variables
load('grad_trans'); load('headmodel_singleshell'); load('trans_matrix'); load('lay');
load('mri_realigned'); clear alldata 

%% Create leadfields in subject's brain warped to MNI space
%Load template sourcemodel - currently 10mm to speed up computation time
%(could easily be 8mm)
load('/Users/44737483/Documents/fieldtrip-20170501/template/sourcemodel/standard_sourcemodel3d10mm.mat');
%load('/Users/44737483/Documents/scripts_mcq/resting_state/Yeo_JNeurophysiol11_MNI152/template_grid_10mm.mat');

template_grid = sourcemodel;
template_grid = ft_convert_units(template_grid,'cm');
clear sourcemodel;

% create the subject specific grid, using the template grid that has just been created
cfg                = [];
cfg.grid.warpmni   = 'yes';
cfg.grid.template  = template_grid;
cfg.grid.nonlinear = 'yes'; % use non-linear normalization
cfg.mri            = mri_realigned;
cfg.grid.unit      ='cm';
cfg.inwardshift = '1.5';
grid               = ft_prepare_sourcemodel(cfg);
figure;ft_plot_mesh((grid.pos(grid.inside,:)));

% Check the mesh is the right shape
grid.pos = ft_warp_apply(trans_matrix,grid.pos);

%% Create leadfield
cfg = [];
cfg.grid = grid;
cfg.headmodel = headmodel_singleshell;
cfg.grad = grad_trans;
cfg.normalize = 'yes'; % May not need this - for rs we probably do
lf = ft_prepare_leadfield(cfg,data);

% make a figure of the single subject headmodel, and grid positions
figure; hold on;
ft_plot_vol(headmodel_singleshell,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
ft_plot_mesh(lf.pos(lf.inside,:));
ft_plot_sens(grad_trans, 'style', 'r*')

% Here we are keeping all parts of the trial to compute the 
% covariance matrix --> common filter
cfg = [];
cfg.covariance = 'yes';
cfg.vartrllength = 2;
cfg.keeptrials = 'no';
cfg.covariancewindow = 'all';
avg = ft_timelockanalysis(cfg,data);

%% Do ya beamforming
% Source reconstruction for the whole 300s
cfg=[];
cfg.method='lcmv';
cfg.grid=lf;
cfg.headmodel=headmodel_singleshell;
cfg.keeptrials = 'no';
cfg.lcmv.keepfilter='yes';
cfg.lcmv.fixedori = 'yes';
sourceavg=ft_sourceanalysis(cfg, avg);

% Make sure your field positions match the template grid
sourceavg.pos=template_grid.pos; % right(?)

% Get AAL atlas
atlas = ft_read_atlas('/Users/44737483/Documents/fieldtrip-20170501/template/atlas/aal/ROI_MNI_V4.nii');
%atlas = ft_read_atlas('/Users/44737483/Documents/scripts_mcq/osl/parcellations/dk_cortical.nii');
%atlas = ft_convert_units(atlas,'cm');

% Interpolate the atlas onto 10mm grid
cfg = []; 
cfg.interpmethod = 'nearest'; 
cfg.parameter = 'tissue'; 
sourcemodel2 = ft_sourceinterpolate(cfg, atlas, template_grid); 

% Create fake VE structure
VE          = [];
VE.label    = sourcemodel2.tissuelabel(1:90)';
VE.time     = data.time;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Maybe switch to 78 regions? Or 40 region one as in the Colclough 
% paper? Need to understand issues of dimensionality...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For each label find vertices of interest, perform PCA and multply by
% spatial filter
for label = 1:90 % only take 90 ROIs (AAL-90)

vertices_of_interest = find(sourcemodel2.tissue == label);
F = cat(1,sourceavg.avg.filter{vertices_of_interest});
[u,s,v] = svd(F*avg.cov*F');
filter = u'*F;

    
for sub=1:(length(data.trial))
    % note that this is the non-filtered "raw" data
    VE.trial{sub}(label,:) = filter(1,:)*data.trial{sub}(:,:);
end

disp(sprintf('Label = %s is DONE', VE.label{label}));

end

save VE VE

% Define frequency bands
alpha_freq = [8 13]; beta_freq = [13 30]; theta_freq = [3 7];
VE_broadband = VE;

% BP filter the VE data
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = alpha_freq;
VE = ft_preprocessing(cfg,VE);

% Play around with correlation 
raw = corr(VE.trial{1,1}');
figure;imagesc(raw);

% Add path of ROInets toolbox
addpath(genpath('/Users/44737483/Documents/scripts_mcq/osl/MEG-ROI-nets'));

% Orthoganalise correlation matrix to remove source leakage
data_orthog = ROInets.remove_source_leakage(VE.trial{1,1}, 'closest');
corrected = corr(data_orthog');
figure;imagesc(corrected);

% Create new VE but replace data with ortoganalised data
VE_corrected = VE;
VE_corrected.trial{1,1} = data_orthog;

% Take Hilbert transform
cfg                         = [];
cfg.hilbert                 = 'abs';
alpha_envelope              = ft_preprocessing(cfg,VE);
alpha_envelope_corrected    = ft_preprocessing(cfg,VE_corrected);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Could improve by using sliding window (?) - YES this is now implemented
% and it seems to increase Z-values... Understand this in greater deal as
% it seems a little bit of a fudge?

% Addpath for osl_movavg
addpath('/Users/44737483/Documents/scripts_mcq/osl/osl-core/');

[envelopedData,newFs] = downsample_envelope(alpha_envelope_corrected.trial{1,1}...
    ,2,alpha_envelope_corrected.time{1});



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial Correlation (try this out later)
covariance              = cov((envelopedData'));
partialCorrSingleRun    = ROInets.convert_precision_to_pcorr(...
                                     pinv(covariance));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create correlation matrix figures to check stuffs 
figure
imagesc(corr(alpha_envelope.trial{1,1}')+diag(nan(90,1)))
axis square
colorbar
title('Envelope correlation before leakage correction');
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

figure
imagesc(corr(alpha_envelope_corrected.trial{1,1}')+diag(nan(90,1)))
axis square
colorbar
title('Envelope correlation after leakage correction')
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

figure
imagesc(corr(resample(alpha_envelope.trial{1,1}',1,200))+diag(nan(90,1)))
axis square
colorbar
title('Envelope correlation before leakage correction downsample to 1Hz')
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

figure
imagesc(corr(resample(alpha_envelope_corrected.trial{1,1}',1,200))+diag(nan(90,1)))
axis square
colorbar
title('Envelope correlation after leakage correction downsample to 1Hz')
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

figure
imagesc(partialCorrSingleRun+diag(nan(90,1)))
axis square
colorbar
title('Partial Correlation Envelope correlation after leakage correction downsample to 1Hz')
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

% Calculate matrix to carry forward (ie. the enveloped downsampled
% orthoganalised data)
matrix_to_use = (corr(envelopedData'));
[Z, z, p] = ROInets.Fisher_r_to_z(matrix_to_use, 90,0);

figure
imagesc(z+diag(nan(90,1)));
axis square
colorbar
title('Envelope correlation after leakage correction downsample to 1Hz')
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Could improve r-to-Z through non-parametric statistics?

% Number of iterations used for statistics
nIter = 100;

% Create AR model parameters from beamformed (but unfiltered) data
[ARmodel.coeffs,           ...
    ARmodel.varianceEstimate, ...
    ARmodel.partial_coeffs  ] = ROInets.estimate_AR_coeffs(VE_broadband.trial{1,1}, ...
    1);

% Create random data using the AR model parameters
clear filter; disp(' Estimating the null data');
randData = filter(1,              ...
                  ARmodel.coeffs, ...
                  sqrt(ARmodel.varianceEstimate) .* randn(length(VE_broadband.trial{1,1}), 90, nIter));

% Loop for each iteration - this takes a while.. any way to speed up?
for iIter = nIter:-1:1,
    disp(iIter);
    % BP filter randomised data
    bpfilt_data = ft_preproc_bandpassfilter(randData(:,:,iIter), 200, alpha_freq)';
    
    % Take hilbert envelope
    bpfilt_data = abs(hilbert(bpfilt_data));
    
    % Downsample envelope using sliding window approach
    [envelopedData_random,newFs] = downsample_envelope(bpfilt_data...
    ,2,VE.time{1});
    
    % Calculate correlation matrix
    rCov  = real(cov(envelopedData_random));
    rCorr = corrcov(rCov,1);
    
    % extract only unique values
    uniqueInd           = triu(true(size(rCorr)), 1);
    
    % Add this to rEmpirical
    rEmpirical(:,iIter) = rCorr(uniqueInd);
end

% Get sigma Z
sigma_z             = std(ROInets.Fisher_r_to_z(rEmpirical(:)));

% Compute r-to-Z correlations using the sigma Z
[Z, z, p] = ROInets.Fisher_r_to_z(matrix_to_use, 90,0,sigma_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
imagesc(z+diag(nan(90,1)))
axis square
hcb=colorbar
title(hcb,'Z-Value')
title('Envelope correlation after leakage correction downsample to 1Hz')
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

cd(dir_name);
edge = z>3.4;
dlmwrite('edge.edge',edge,'\t');

% Addpath to BrainNet Viewer
addpath(genpath('/Users/44737483/Documents/scripts_mcq/BrainNet-Viewer-master'));

% Addpath to CaptureFigVid
addpath(genpath('/Users/44737483/Documents/scripts_mcq/coreg/CaptureFigVid'));

% Make Figure and Save Video
BrainNet_MapCfg('/Users/44737483/Documents/scripts_mcq/BrainNet-Viewer-master/Data/SurfTemplate/BrainMesh_ICBM152.nv','/Users/44737483/Documents/scripts_mcq/BrainNet-Viewer-master/Data/ExampleFiles/AAL90/Node_AAL90.node',[dir_name '/edge.edge']);
OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
CaptureFigVid([0,0; 360,0], 'alpha_corr_connectivity_Z_greater_3_point4',OptionZ)



