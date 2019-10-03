%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ROInet: A Matlab package for performing leakage-robust network inference
% between ROIs in MEG data

% The methodology used in this pipeline is set out in:
%
% Colclough, G. L., Brookes, M., Smith, S. M. and Woolrich, M. W.,
% "A symmetric multivariate leakage correction for MEG connectomes,"
% NeuroImage 117, pp. 439-448 (2015)
%
% We are grateful for Colclough and colleagues at OHBA for making their
% code openly available for re-use. Please check copyright information in
% OSL for more information.
%
% Data is from an Elekta-Neuromag (306 channel) MEG system collected at 
% the Aston Brain Centre, Birmingham, UK. ASD/control participants are
% compared.
%
% The script uses VE in Fieldtrip format but various MEG-ROInets scripts
% to perform the key analysis steps
%
% Parcellation is based on 39 region parcellation scheme discussed in 
% Colclough, Giles L., Mark W. Woolrich, P. K. Tewarie, Matthew J. Brookes,
% Andrew J. Quinn, and Stephen M. Smith. "How reliable are MEG 
% resting-state connectivity metrics?." NeuroImage 138 (2016): 284-293.
%
% Users will need to download OSL (https://github.com/OHBA-analysis).
% However please don't initiate OSL - this script will add the relevent
% bits of code for you.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subject = {'0401','0402','0403','0404','0405','0406','0407','0408','0409',... %ASD
    '0411','0413','0414','0415','0416','0417','0418',...
    '1401','1402','1404','1405','1406','1407','1408','1410','1411','1412',... %Control
    '1413','1415','1416','1417'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subject = {'2660','2708','2735','2852','2857','2870','2880','2881'...   % ASD
    '2704','2721','2842','2850','2851','2877','2883','2890','2891'};    % Control

  % Define frequency bands
    frequency_bands = {[3 7],[8 13],[13 30],[31 48],[35 70] % Top row 
        'theta','alpha','beta','gamma','broadband_gamma'};
    
    
for i = 1:length(subject)
    % Close figures from the previous loop
    close all force
    
    % Go to directory and load VE 
    dir_name = ['/Users/44737483/Dropbox/MEG_data/rs_data/' subject{i}];
    cd(dir_name); load('VE.mat');
    % Get colourmap
    ft_hastoolbox('brewermap', 1);
    
    %% Rearrange VE to correct order

    LABELS = {'Left Visual Cortex'; ...
          'Right Visual Cortex'; ...
          'Left Occipital Lobe'; ...
          'Right Occipital Lobe'; ...
          'Left Somatosensory Cortex'; ...
          'Right Somatosensory Cortex'; ...
          'Left Temporal Lobe'; ...
          'Right Temporal Lobe'; ...
          'Right Parietal Lobe'; ...
          'Left Motor Cortex'; ...
          'Left Parietal Lobe'; ...
          'Right Parietal Lobe'; ...
          'Left Occipital Lobe'; ...
          'Right Occipital Lobe'; ...
          'Left Parietal Lobe'; ...
          'Left Temporal Lobe'; ...
          'Right Temporal Lobe'; ...
          'Left Motor Cortex'; ...
          'Right Motor Cortex'; ...
          'Right Parietal Lobe'; ...
          'Right Parietal Lobe'; ...
          'Left Parietal Lobe'; ...
          'Right Parietal Lobe'; ...
          'Left Frontal Lobe'; ...
          'Right Frontal Lobe'; ...
          'Left Visual Cortex'; ...
          'Right Visual Cortex'; ...
          'Left Frontal Lobe'; ...
          'Right Frontal Lobe'; ...
          'Right Frontal Lobe'; ...
          'Left Frontal Lobe'; ...
          'Left Frontal Lobe'; ...
          'Right Frontal Lobe'; ...
          'Left Temporal Lobe'; ...
          'Right Temporal Lobe'; ...
          'Left Frontal Lobe'; ...
          'Right Frontal Lobe'; ...
          'Posterior Cingulate'; ... % only with PCC included
          'Medial Frontal Regions'};
      
      % Number the labels and replace spaces with underscores
      for lab = 1:length(LABELS)
        LABELS{lab,1} = sprintf('%d%s',lab,strrep(LABELS{lab,1},' ','\_'));
      end
      
      NEW_ORDER = [32 36 31 24 28 5 10 18 11 15 22 1 26 3 13 34 7 16 ...
          17 8 35 14 4 27 2 20 23 21 9 12 19 6 29 30 25 33 37 38 39];
      
      VE.trial{1,1} = (VE.trial{1,1}(NEW_ORDER,:));
      
      % Reorder LABELS
      LABELS = LABELS(NEW_ORDER,:);
      VE.label = LABELS;

    %% Add paths to
    
    % Add path of ROInets toolbox
    addpath(genpath('/Users/44737483/Documents/scripts_mcq/osl/MEG-ROI-nets'));
    
    % Save broadband version of the VE for later
    VE_broadband = VE;
    
    % Create AR model parameters from beamformed (but unfiltered) data
    [ARmodel.coeffs,           ...
        ARmodel.varianceEstimate, ...
        ARmodel.partial_coeffs  ] = ROInets.estimate_AR_coeffs(VE_broadband.trial{1,1}, ...
        1);
    
    % Number of iterations used for r-to-z transform
    nIter = 20;
    
    % Create random data using the AR model parameters
    clear filter; fprintf(' Estimating the null data\n');
    randData = filter(1,              ...
        ARmodel.coeffs, ...
        sqrt(ARmodel.varianceEstimate) .* randn(length(VE.label),length(VE_broadband.trial{1,1}),nIter));
    
    %% Loop for frequency bands
    for freq_b = 1:length(frequency_bands(1,:))
    
        fprintf(' Processing frequency band -  %s  - for subject %s \n',frequency_bands{2,freq_b},subject{i});

    % BP filter the VE data
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = frequency_bands{1,freq_b};
    VE_filtered = ft_preprocessing(cfg,VE);
    
    % Orthoganalise data to remove source leakage
    fprintf(' Orthoganalising the data using "closest" algorithm". May take a while\n');
    data_orthog = (ROInets.remove_source_leakage(VE_filtered.trial{1,1}, 'closest'));
    
    % Create new VE but replace data with ortoganalised data
    VE_corrected = VE_filtered;
    VE_corrected.trial{1,1} = data_orthog;
    
    % Create Figure to show the effect of orthoganalisation
    figure; subplot(2,1,1);
    plot(VE_filtered.trial{1,1}(1,1:1000),'LineWidth',3); hold on;
    plot(VE_filtered.trial{1,1}(2,1:1000),'LineWidth',3);
    title(sprintf('Raw data %s band',frequency_bands{2,freq_b}),'Interpreter','none'); set(gca,'FontSize',20);
    subplot(2,1,2);
    plot(VE_corrected.trial{1,1}(1,1:1000),'LineWidth',3); hold on;
    plot(VE_corrected.trial{1,1}(2,1:1000),'LineWidth',3);
    title(sprintf('Corrected data %s band',frequency_bands{2,freq_b}),'Interpreter','none'); set(gca,'FontSize',20); drawnow;
    
    % Take Hilbert transform
    fprintf('Taking the Hilbert Transform\n');
    cfg                   = [];
    cfg.hilbert           = 'abs';
    envelope              = ft_preprocessing(cfg,VE_filtered);
    envelope_corrected    = ft_preprocessing(cfg,VE_corrected);
    

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Could improve by using sliding window (?) - YES this is now implemented
    % and it seems to increase Z-values... Understand this in greater deal as
    % it seems a little bit of a fudge?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Addpath for osl_movavg
    addpath(genpath('/Users/44737483/Documents/scripts_mcq/osl/osl-core/'));
    addpath('/Users/44737483/Documents/scripts_mcq/resting_state');
    
    % Envelope the data using a 2 second window
    [envelopedData,newFs] = downsample_envelope(envelope_corrected.trial{1,1}...
        ,2,envelope_corrected.time{1});
    
    fprintf('The new sampling frequency of the envelope is %d Hz\n',newFs);
    
    % Run Correlations
    Settings.Regularize.do            = false;
    
    mat = ROInets.run_correlation_analysis(data_orthog,     ...
        envelopedData,      ...
        Settings.Regularize);
    
%     % Make a figure to show corrected and uncorrected r-scores
%     figure; subplot(1,2,1);
%     raw = corr(VE_filtered.trial{1,1}');
%     imagesc(raw+diag(nan(length(VE.label),1)));
%     title('Raw Correlations');
%     axis square
%     hcb=colorbar;
%     colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
%     hold on;
%     subplot(1,2,2);
%     imagesc(mat.envCorrelation+diag(nan(length(VE.label),1)))
%     axis square
%     hcb=colorbar;
%     title('Corrected r-Value')
%     colormap(flipud(brewermap(64,'RdBu'))); % change the colormap
%     drawnow;
%     
    %% Convert r-to-Z Using Null Data
   
    % Loop for each iteration - this takes a while.. any way to speed up?
    for iIter = nIter:-1:1,
        disp(iIter);
        % BP filter randomised data
        bpfilt_data = ft_preproc_bandpassfilter(randData(:,:,iIter), VE.fsample, frequency_bands{1,freq_b});


        % Take hilbert envelope
        bpfilt_data = abs(hilbert(bpfilt_data));
                
        % Downsample envelope using sliding window approach
        [envelopedData_random,newFs] = downsample_envelope(bpfilt_data,...
        2,envelope_corrected.time{1});
        
        % Calculate correlation matrix
        rCov  = real(cov(envelopedData_random));
        rCorr = corrcov(rCov,1);
        
        % extract only unique values
        uniqueInd           = triu(true(size(rCorr)), 1);
        
        % Add this to rEmpirical
        rEmpirical(:,iIter) = rCorr(uniqueInd);
    end
    
    % Get sigma Z (std of null r values)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % N.B. Also need to apply the same for partial correlations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    sigma_z             = std(ROInets.Fisher_r_to_z(rEmpirical(:)));
    sigma_to_use = [];
    sigma_to_use.z = sigma_z;
    sigma_to_use.z_partial = sigma_z;
    sigma_to_use.partial_reg = sigma_z;
        
    % Get p-value and save
    [~,~,p]                 = ROInets.Fisher_r_to_z(mat.envCorrelation,   ...
        mat.nSamples,         ...
        0, ...
        sigma_z);
    
    p = p+diag(nan(length(VE.label),1));
    
    cd(dir_name);
    save(sprintf('p_values_%s',frequency_bands{2,freq_b}),'p');
    
    % Convert r-to-z based on sigma
    mat2 = ROInets.convert_correlations_to_normal_variables(mat, ...
        sigma_to_use,      ...
        false);
    
    %cd(dir_name);
    save(sprintf('corr_matrix_%s',frequency_bands{2,freq_b}),'mat2');
    
    %% Make a figure to show corrected Z-scores
    figure
    imagesc(mat2.envCorrelation_z+diag(nan(length(VE.label),1)))
    axis square; hcb=colorbar;
    lim = max(max(mat2.envCorrelation_z+diag(nan(length(VE.label),1))));
    zlim([-lim lim]); caxis([-lim lim])
    title(hcb,'Z-Value')
    colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
    title(sprintf('Subject %s',subject{i}));
    set(gca,'FontSize',20);
    
    % Plot labels for parcels
    yticks([1 6 9 12 16 19 22 26 31 33 38 39]);
    xticks([]);

    yticklabels(LABELS([1 7 9 12 16 19 24 26 31 33 38 39]));
    
    % Change yLabel FontSize
    yl = get(gca,'YLabel');
    ylFontSize = get(yl,'FontSize');
    yAX = get(gca,'YAxis');
    set(yAX,'FontSize', 15); set(yl, 'FontSize', 20);
    
    % Draw and save
    set(gcf, 'Position', [686 431 840 539])
    drawnow; fprintf(' Saving figure to file\n');
    print(sprintf('Z_subject_%s_%s',subject{i},...
        frequency_bands{2,freq_b}),'-dpng','-r300');
    
    %% Create Circular Graph
%     addpath(genpath('/Users/44737483/Documents/scripts_mcq/resting_state/circularGraph-master'));
%     
%     % Make diagonal = 0
%     z = mat2.envCorrelation_z + diag(nan(90,1));
%     
%     % Mask out everything but top 1% of connections
%     V = sort(z(:), 'descend');
%     indices = find(z< min(V(91:91+(length(V)/100))));
%     z(indices) = 0;
%         
%     atlas = ft_read_atlas('/Users/44737483/Documents/fieldtrip-20170501/template/atlas/aal/ROI_MNI_V4.nii');
%     
%     for f = 1:90
%         atlas.tissuelabel{1,f} = strrep(atlas.tissuelabel{1,f}, '_', ' ');
%     end
%    
%     figure; circularGraph(z,'Label',atlas.tissuelabel);
% % 
%     % 
%     cd(dir_name);
%     edge = z;
%     dlmwrite('edge.edge',edge,'\t');
%    
    
%% Make video of connections passing p<.2

% try
% z = zeros(size(mat2.envCorrelation));
% 
% z(find(p<.05)) = mat2.envCorrelation(find(p<.05));
% 
% z = z+diag(nan(39,1));
% 
% z(isnan(z)) = 1; %z_all = z_all+z;
% 
%     cd(dir_name);
%     edge = z;
%     dlmwrite('edge.edge',edge,'\t');
% 
%      % Addpath to BrainNet Viewer
%      addpath(genpath('/Users/44737483/Documents/scripts_mcq/BrainNet-Viewer-master'));
%     
%     % Addpath to CaptureFigVid
%     addpath(genpath('/Users/44737483/Documents/scripts_mcq/coreg/CaptureFigVid'));
%     
%     % Make Figure and Save Video
%     BrainNet_MapCfg('BrainMesh_Ch2_smoothed.nv','/Users/44737483/Dropbox/MEG_data/rs_data/atlas_39_binary_giles.node',[dir_name '/edge.edge']);
%     OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
%     CaptureFigVid([0,0; 360,0], 'connectivity_p_0_02',OptionZ)
%     
% catch
%     fprintf('Could not make video of connections\n');
% end
    
% %     
%     %% Create correlation matrix figures to check stuffs
%     figure;
%     subplot(2,2,1);imagesc(corr(envelope.trial{1,1}')+diag(nan(90,1)))
%     axis square
%     colorbar
%     title({'Envelope correlation before'; 'leakage correction'});
%     ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
%     colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
%     
%     subplot(2,2,2);
%     imagesc(corr(envelope_corrected.trial{1,1}')+diag(nan(90,1)))
%     axis square
%     colorbar
%     title({'Envelope correlation after'; 'leakage correction'})
%     colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
%     
%     subplot(2,2,3);
%     imagesc(corr(resample(envelope.trial{1,1}',1,200))+diag(nan(90,1)))
%     axis square
%     colorbar
%     title({'Envelope correlation before'; 'leakage correction downsample to 1Hz'})
%     colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
%     
%     subplot(2,2,4);
%     imagesc(corr(resample(envelope_corrected.trial{1,1}',1,200))+diag(nan(90,1)))
%     axis square
%     colorbar
%     title({'Envelope correlation after'; 'leakage correction downsample to 1Hz'})
%     colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
%     
    end
end










%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot connectivity matrices for control and ASD groups
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Addpath to BrainNet Viewer
addpath(genpath('/Users/44737483/Documents/scripts_mcq/BrainNet-Viewer-master'));

% Addpath to CaptureFigVid
addpath(genpath('/Users/44737483/Documents/scripts_mcq/coreg/CaptureFigVid'));

% Define frequency bands
frequency_bands = {[3 7],[8 13],[13 35],[31 48],[35 70] % Top row 
        'theta','alpha','beta','gamma','broadband gamma'};

    for freq_b = 1:length(frequency_bands(1,:))
    
%% ASD Group

% subject = {'0401','0402','0403','0404','0405','0406','0407','0408','0409',... %ASD
%     '0411','0413','0414','0415','0416','0417','0418','2660','2708','2735'...
%     ,'2852','2857','2870','2880','2881'};


subject = {'2660','2708','2735'...
    ,'2852','2857','2870','2880','2881'};

z_all_ASD = zeros(39); %ASD_num = [1:1:16];

for i = 1:length(subject)
    % Go to directory and load VE
    dir_name = ['/Users/44737483/Dropbox/MEG_data/rs_data/' subject{i}];
    cd(dir_name); load(sprintf('corr_matrix_%s.mat',frequency_bands{2,freq_b})); %load('p_values_alpha.mat');
    
    z = mat2.envCorrelation_z;
    % Make diagnonal = 1
    z = z+diag(nan(39,1));
    z(isnan(z)) = 1;
    z_all_ASD = z_all_ASD+z;
    
end

%% Do group statistics
z_all_ASD = z_all_ASD./length(subject);

sd_of_means = 1.0 ./ sqrt(length(subject));
z_all_ASD = z_all_ASD ./ sd_of_means;

[h, z_thresh, q] = ROInets.false_discovery_rate(z_all_ASD, 0.05,'dep');
z_all_ASD_thresh = z_all_ASD;
z_all_ASD_thresh(z_all_ASD_thresh<3.1) = 0;

fprintf('The threshold for Z (FDR corrected, p<.05) is %3d\n'.z_thresh);

cd('/Users/44737483/Dropbox/MEG_data/rs_data/Group');

% try
%     % Make video
%     dlmwrite(sprintf('edge_ASD_%s.edge',frequency_bands{2,freq_b}),z_all_ASD_thresh,'\t');
%     
%     % Make Figure and Save Video
%     BrainNet_MapCfg('BrainMesh_Ch2_smoothed.nv','/Users/44737483/Dropbox/MEG_data/rs_data/atlas_39_binary_giles.node',sprintf('edge_ASD_%s.edge',frequency_bands{2,freq_b}));
%     OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
%     CaptureFigVid([0,0; 360,0], sprintf('connectivity_ASD_%s',frequency_bands{2,freq_b}),OptionZ)
% catch
% end

%% Make ASD-Group Matrix Diagram
figure
imagesc(z_all_ASD+diag(nan(39,1)))
title(sprintf('ASD group Z-value %s band',frequency_bands{2,freq_b}))
axis square
hcb=colorbar
lim = max(max(z_all_ASD+diag(nan(39,1))));
zlim([-lim lim])
caxis([-lim lim])
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
set(gca,'FontSize',20);

% Plot labels for parcels
yticks([1 6 9 12 16 19 22 26 31 33 38 39]);
xticks([]);

yticklabels(LABELS([1 7 9 12 16 19 24 26 31 33 38 39]));

% Change yLabel FontSize
yl = get(gca,'YLabel');
ylFontSize = get(yl,'FontSize');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 15); set(yl, 'FontSize', 20);

% Draw and save
set(gcf, 'Position', [686 431 840 539])
drawnow; fprintf(' Saving figure to file\n');
print(sprintf('ASD Group Z-Value %s Band',frequency_bands{2,freq_b}),'-dpng','-r300');

%% Control
% subject = {'1401','1402','1404','1405','1406','1407','1408','1410','1411','1412',... %Control
%     '1413','1415','1416','1417','2704','2721','2842','2850','2851','2877',...
%     '2883','2890','2891'};

subject = {'2704','2721','2842','2850','2851','2877',...
    '2883','2890','2891'};

z_all_control = zeros(39);

for i = 1:length(subject)
    % Go to directory and load VE
    dir_name = ['/Users/44737483/Dropbox/MEG_data/rs_data/' subject{i}];
    cd(dir_name); load(sprintf('corr_matrix_%s.mat',frequency_bands{2,freq_b})); %load('p_values_broadband gamma.mat');
    
    z = mat2.envCorrelation_z;
    z = z+diag(nan(39,1));
    z(isnan(z)) = 1; %z_all = z_all+z;  
    z_all_control = z_all_control+z;
end

z_all_control = z_all_control./length(subject);

sd_of_means = 1.0 ./ sqrt(length(subject));
z_all_control = z_all_control ./ sd_of_means;

[h, z_thresh, q] = ROInets.false_discovery_rate(z_all_control, 0.05,'dep');

z_all_control_thresh = z_all_control;
z_all_control_thresh(z_all_control_thresh<3.1) = 0;

cd('/Users/44737483/Dropbox/MEG_data/rs_data/Group');

% try
%     % Make video
%     dlmwrite(sprintf('edge_control_%s.edge',frequency_bands{2,freq_b}),z_all_control_thresh,'\t');
%     
%     % Make Figure and Save Video
%     BrainNet_MapCfg('BrainMesh_Ch2_smoothed.nv','/Users/44737483/Dropbox/MEG_data/rs_data/atlas_39_binary_giles.node',[cd '/edge_control_' frequency_bands{2,freq_b} '.edge']);
%     OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
%     CaptureFigVid([0,0; 360,0], sprintf('connectivity_control_%s',frequency_bands{2,freq_b}),OptionZ)
% catch
% end

%% Make Matrix Diagram
cd('/Users/44737483/Dropbox/MEG_data/rs_data/Group');
figure
imagesc(z_all_control+diag(nan(39,1)))
title(sprintf('Control group Z-value %s band',frequency_bands{2,freq_b}))
axis square
hcb=colorbar
lim = max(max(z_all_control+diag(nan(39,1))));
zlim([-lim lim])
caxis([-lim lim])
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
set(gca,'FontSize',20);

% Plot labels for parcels
yticks([1 6 9 12 16 19 22 26 31 33 38 39]);
xticks([]);

yticklabels(LABELS([1 7 9 12 16 19 24 26 31 33 38 39]));

% Change yLabel FontSize
yl = get(gca,'YLabel');
ylFontSize = get(yl,'FontSize');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 15); set(yl, 'FontSize', 20);

% Draw and save
set(gcf, 'Position', [686 431 840 539])
drawnow; fprintf(' Saving figure to file\n');
print(sprintf('Control Group Z-Value %s Band',frequency_bands{2,freq_b}),'-dpng','-r300');

%% Difference between ASD/control

diff = z_all_ASD-z_all_control;

%[h, z_thresh, q] = ROInets.false_discovery_rate(diff, 0.05,'dep');

%diff(diff<1.96) = 0;

cd('/Users/44737483/Dropbox/MEG_data/rs_data/Group');
figure
imagesc(diff+diag(nan(39,1)))
axis square
hcb=colorbar
lim = max(max(diff+diag(nan(39,1))));
zlim([-lim lim])
caxis([-lim lim])
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
title(sprintf('%s band ASD>Control Z-Value',frequency_bands{2,freq_b}));
set(gca,'FontSize',20);
% Plot labels for parcels
yticks([1 6 9 12 16 19 22 26 31 33 38 39]);
xticks([]);

yticklabels(LABELS([1 7 9 12 16 19 24 26 31 33 38 39]));

% Change yLabel FontSize
yl = get(gca,'YLabel');
ylFontSize = get(yl,'FontSize');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 15); set(yl, 'FontSize', 20);

% Draw and save
set(gcf, 'Position', [686 431 840 539])
drawnow; 
print(sprintf('ASD versus Control Difference %s Band',frequency_bands{2,freq_b}),'-dpng','-r300');

end





%% Videos

% Make video
edge_ASD = z_all_ASD;
dlmwrite('edge_ASD.edge',edge_ASD,'\t');



% Make Figure and Save Video
BrainNet_MapCfg('BrainMesh_Ch2_smoothed.nv','/Users/44737483/Dropbox/MEG_data/rs_data/atlas_39_binary_giles.node',[cd '/edge_ASD.edge']);
OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
CaptureFigVid([0,0; 360,0], 'connectivity_ASD_gamma',OptionZ)


% Make video
edge = z_all_control;
dlmwrite('edge_control.edge',edge,'\t');

% Make Figure and Save Video
BrainNet_MapCfg('BrainMesh_Ch2_smoothed.nv','/Users/44737483/Dropbox/MEG_data/rs_data/atlas_39_binary_giles.node',[cd '/edge_control.edge']);
OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
CaptureFigVid([0,0; 360,0], 'connectivity_gamma_theta',OptionZ)















z_all_ASD = z_all_ASD./16;
z_all_control = z_all_control./14;

z_diff = z_all_ASD-z_all_control;

% [a,b,p]                 = ROInets.Fisher_r_to_z(z_diff,   ...
%         598,         ...
%         0, ...
%         0.4);

figure
imagesc(z_diff+diag(nan(39,1)))
axis square
hcb=colorbar
lim = max(max(z_all_ASD+diag(nan(39,1))));
zlim([-lim lim])
caxis([-lim lim])
title(hcb,'Group r-Value')
colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
set(gca,'FontSize',20);



%% Group Analyis








%% NBS
% Make design matrix separatly to data
cd('/Users/44737483/Dropbox/MEG_data/rs_data')

total = 30;
ASD = 16; control = 14;

design = zeros(total,2);
design(1:ASD,1) = 1; 
design(total-control+1:end,2) = 1;
dlmwrite('designMatrix.txt',design,'delimiter',' ');
save designMatrix design



cd('/Users/44737483/Documents/scripts_mcq/resting_state/stats')




    
%% Group analysis

subject = {'2660','2708','2721','2735','2842','2850','2851','2852',...
    '2857','2870','2877','2880','2881','2883'};

for i = 1:length(subject)
    % Go to directory and load VE 
    dir_name = ['/Users/44737483/Documents/scripts_mcq/resting_state/' subject{i}];
    cd(dir_name); load('mat2.mat');
    z = mat2.envCorrelation_z+diag(nan(90,1));
    %z(z<2.1) = 0; 
    z_all = z_all+z;
end

z_all = zeros(90);

z_all(z_all<2) = 0; 
z_all = z_all./14;

figure; circularGraph(z_all,'Label',atlas.tissuelabel(1:90));

    
