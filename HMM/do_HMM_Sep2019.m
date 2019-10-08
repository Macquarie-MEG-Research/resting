%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up OSL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath

cd('/Volumes/Robert T5/osl/osl-core');
osl_startup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HCP_lib = '/Volumes/Robert T5/Robert_HCP/';
save_path = '/Volumes/Robert T5/RS_HCP_VE/';

addpath('/Volumes/Robert T5/HCP_MEG_RS_process');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subject List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subject = {'100307','102816','105923','106521','108323','109123',...
    '111514','112920','113922','116524','116726','133019','140117',...
    '146129','149741','153732','154532','158136',...
    '162026','162935','164636','166438','169040','172029','174841',...
    '175237','175540','177746','179245','181232',...
    '187547','189349','191033','191437','191841','192641'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load all three runs of data and save to SPM format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for sub = 1:length(subject)
    fprintf('Reorganising data for %s\n',subject{sub});
    try
        cd([save_path subject{sub}]);
        
        for i = 1:1
            load(['VE_rs_' num2str(i) '.mat']);
            disp('FT --> SPM...');
            spm_eeg_ft2spm(VE,['VE_SPM_' num2str(i)]);
            clear VE
        end
        
    catch
        warning('!!!');
        fprintf('Subject %s could not be processed\n',subject{sub});
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create cell array with path to data-files (only first run)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 1;
data_files = [];
for sub = 1:length(subject)
    for i = 1:1
        data_files{count} = [save_path subject{sub} '/VE_SPM_' ...
            num2str(i) '.mat'];
        count = count + 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up options for the rest of the analysis in OSL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmmdir = '/Volumes/Robert T5/HMM_test';

OSLDIR = getenv('OSLDIR');

todo.envelope = 1;
todo.concat   = 1;
todo.infer    = 1;
todo.output   = 1;
    
options.envelope.windowsize = 0.1;
options.concat.log          = 0;
options.concat.norm_subs    = 1;
options.concat.pcadim       = 40;
options.hmm.nstates         = 10;
options.hmm.nreps           = 1;
options.output.method       = 'pcorr';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Again this is some housekeeping for later, including where to save data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HMMresults = [];
statemaps  = [];

if ~isdir(hmmdir)
  mkdir(hmmdir); 
end

filenames = [];
filenames.hmm = '/Volumes/Robert T5/HMM_test/hmm.mat';
filenames.output = '/Volumes/Robert T5/HMM_test/hmm';
filenames.concat = '/Volumes/Robert T5/HMM_test/env_concat.mat';
mkdir([filenames.output '/envelope/']);
for i = 1:length(subject)
    filenames.envelope{1,i} = [filenames.output ...
        '/envelope/VE_SPM_1_subject_' num2str(subject{i})];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default todo settings
try todo.envelope = todo.envelope; catch, todo.envelope = 1; end
try todo.concat   = todo.concat;   catch, todo.concat   = 1; end
try todo.hmm      = todo.hmm;      catch, todo.hmm      = 1; end

% Default prepare settings
try envelope_do = options.prepare.envelope;   catch, envelope_do = 1;   end

% Default envelope settings
try windowsize = options.envelope.windowsize; catch, windowsize = 0.1;  end
try multiband  = options.envelope.multiband;  catch, multiband  = [];   end

% Default concatenation settings
try logtrans  = options.concat.log;        catch, logtrans   = 0;  end
try norm_subs = options.concat.norm_subs;  catch, norm_subs  = 1;  end
try pcadim    = options.concat.pcadim;     catch, pcadim     = 40; end
try whiten    = options.concat.whiten;     catch, whiten     = 1;  end

% Default HMM settings
try nstates = options.hmm.nstates; catch, nstates = 8; end
try nreps   = options.hmm.nreps;   catch, nreps   = 5; end
try use_old_tbx = options.hmm.use_old_hmm_tbx;  catch, use_old_tbx = 0; end
% Default output settings
try output_method = options.output.method; catch, output_method = 'pcorr'; end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Envelopes using OSL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.envelope
  
  for subnum = 1:length(data_files)
    
    [pathstr,filestr] = fileparts(data_files{subnum});
  
    disp(['Computing envelope data for ' filestr]);
   
    % ---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---
    if ~isempty(multiband)
        for f = 1:numel(multiband)
            D = spm_eeg_load(data_files{subnum});
            current_montage = montage(D,'getindex');
            D = montage(D,'switch',0);
            spm_eeg_filter(struct('D',D,'band','bandpass','freq',multiband{f}));
            D = spm_eeg_load(prefix(data_files{subnum},'f'));
            D = montage(D,'switch',current_montage); D.save;
            
            S = [];
            S.D = fullfile(D.path,D.fname);
            S.winsize = windowsize;
            
            D = osl_hilbenv(S);
            move(D,strrep(filenames.envelope{subnum},'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']));
            D = spm_eeg_load(prefix(data_files{subnum},'f'));
            D.delete;
            disp(['Saving envelope data for ' filestr ' to ' filenames.envelope{subnum} ' band ' num2str(f)])
            clear D
        end
    else % Single frequency band
        D = spm_eeg_load(data_files{subnum});
        S = [];
        S.D = fullfile(D.path,D.fname);
        S.winsize = windowsize;
        
        D = osl_hilbenv(S);
        move(D,filenames.envelope{subnum});  
        disp(['Saving envelope data for ' filestr ' to ' filenames.envelope{subnum}])
        clear D
    end
    % ---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---8<---
    
    
  end
  
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatenate the data for HMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.concat || (todo.infer && ~exist(filenames.concat,'file'))
    
    for f = 1:max(numel(multiband),1)
        
        % Load subjects and concatenate:
        env_concat = [];
        subj_inds = [];
        C = 0;
        
        for subnum = 1:length(data_files)
            
            if ~isempty(multiband)
                D = spm_eeg_load(strrep(filenames.envelope{subnum},'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']));
            else
                D = spm_eeg_load(filenames.envelope{subnum});
            end
           
            % Find timepoints without NaNs (could be improved)
           tttt = isnan(D(:,:));
           [row,col] = find(double(tttt(1,:))== 0);
           samples2use = squeeze(col);
           clear row col tttt
           
%             tbad = ~good_samples(D);
%             samples2use = find(~tbad);
            
            env = D(:,samples2use); %#ok - SPM doesn't like logical indexing
            
            if logtrans
                env = log10(env);
            end
            
            if norm_subs
                env = demean(env,2)./std(env(:));
                %env = normalise(env,2);
            end
            
            env_concat = [env_concat,env];
            subj_inds = [subj_inds,subnum*ones(1,size(env,2))];
            C = C + env*permute(env,[2,1]);
        end
        C = C ./ (length(subj_inds)-1);
        clear env
        
        
        % PCA + whitening - below is equivalent to fastica whitening code but much faster
        [allsvd,MixingMatrix] = eigdec(C,39);
        if whiten
            MixingMatrix = diag(1./sqrt(allsvd)) * MixingMatrix';
        end
        hmmdata =  (MixingMatrix * env_concat)';
        fsample = D.fsample;
        
        if ~isempty(multiband)
            savestr = strrep(filenames.concat,'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']);
            disp(['Saving concatenated envelope data to ' savestr])
            save(savestr,'hmmdata','MixingMatrix','fsample','subj_inds')
        else
            disp(['Saving concatenated envelope data to ' filenames.concat])
            save(filenames.concat,'hmmdata','MixingMatrix','fsample','subj_inds')
        end
        
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN HMM (!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if todo.infer
    
    if ~isempty(multiband)
        hmmdata = [];
        freq_inds = [];
        MixingMatrix = cell(numel(multiband),1);
        for f = 1:numel(multiband)
            tmp = load(strrep(filenames.concat,'.mat',['_',strrep(num2str(multiband{f}),'  ','_') 'Hz.mat']));
            hmmdata = [hmmdata tmp.hmmdata];
            MixingMatrix{f} = tmp.MixingMatrix;
        end
        subj_inds = tmp.subj_inds;
        fsample = tmp.fsample;
    else
        load(filenames.concat)
    end
    
    %hmm = osl_hmm_infer(hmmdata,struct('K',nstates,'order',0,'Ninits',nreps,'Hz',fsample,'zeromean',0));
    
    % Switch between Iead's & Diego's HMM toolboxes
    if use_old_tbx
        rmpath(genpath(fullfile(OSLDIR,'osl2/osl_hmm_toolbox/HMM-MAR')));
        addpath(genpath(fullfile(OSLDIR,'hmmbox_4_1')));

        if envelope_do 
            hmm = ABhmm_infer(hmmdata,nstates,nreps);
        else
            hmm = ABhmm_infer(hmmdata,nstates,nreps,'constrain_mean');
        end;
        addpath(genpath(fullfile(OSLDIR,'osl2/osl_hmm_toolbox/HMM-MAR')));
        rmpath(genpath(fullfile(OSLDIR,'hmmbox_4_1')));

    else

        rmpath(genpath(fullfile(OSLDIR,'osl2/hmmbox_4_1')));
        addpath(genpath(fullfile(OSLDIR,'osl_hmm_toolbox/HMM-MAR')));

        hmm = osl_hmm_infer(hmmdata,struct('K',nstates,'order',0,...
            'Ninits',nreps,'Hz',fsample,'zeromean',~envelope_do));
        %hmm = osl_hmm_infer(hmmdata,struct('K',nstates,'order',0,'Ninits',nreps,'Hz',fsample,'zeromean',false));
        addpath(genpath(fullfile(OSLDIR,'osl2/hmmbox_4_1')));
        rmpath(genpath(fullfile(OSLDIR,'osl_hmm_toolbox/HMM-MAR')));

    end
    
    
    hmm.MixingMatrix = MixingMatrix;
    hmm.fsample = fsample;
    
    % Save results
    disp(['Saving inferred HMM to ' filenames.hmm])
    save(filenames.hmm,'hmm','subj_inds')
    
    HMMresults = filenames.hmm;

end

%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       O U T P U T   R E S U L T S                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hmmdir = '/Volumes/Robert T5/HMM_test';
cd(hmmdir);
load('hmm.mat');

statemaps = [filenames.output,'_',output_method];

% We first need to load voxel data and reorganise into FT --> SPM format

for sub = 1:length(subject)
    fprintf('Processing subject %s\n',subject{sub});
    cd([save_path subject{sub}]);
    disp('loading voxel data...');
    load('VE_wb_1.mat');
    
    disp('FT --> SPM...');
    spm_eeg_ft2spm(VE_wb,'VE_SPM_wholebrain');
    clear VE_wb
end


% Now we need to compute envelope as per ROI-analysis...
% This takes a *while*

for sub = 1:1:length(subject)
    fprintf('Processing subject %s\n',subject{sub});
    
    cd([save_path subject{sub}]);
    
    D = spm_eeg_load('VE_SPM_wholebrain');
    S = [];
    S.D = fullfile(D.path,D.fname);
    S.winsize = windowsize;
    
    disp('Computing envelope...');
    D = osl_hilbenv(S);
    
    move(D,'VE_SPM_wholebrain.mat');
    
end

% Load data and partially correlate with whole-brain envelopes
stat = zeros(5798,hmm.K);

for sub = 1:length(subject)
    fprintf('Processing subject %s\n',subject{sub});
    cd([save_path subject{sub}]);
    
    D = spm_eeg_load('VE_SPM_wholebrain');

    %
    disp('Finding NaNs');
    tttt = isnan(D(1,:));
    [row,col] = find(double(tttt(1,:))== 0);
    samples2use = squeeze(col);
    clear row col tttt
    
    %             tbad = ~good_samples(D);
    %             samples2use = find(~tbad);
    
    env = D(:,samples2use); %#ok - SPM doesn't like logical indexing
    
    if logtrans
        env = log10(env);
    end
    
    if norm_subs
        env = demean(env,2)./std(env(:));
        %env = normalise(env,2);
    end
    
    hmm_sub = hmm; hmm_sub.statepath = hmm.statepath(subj_inds==sub);
    stat = stat + osl_hmm_statemaps(hmm_sub,env,0,'pcorr');
end

% Divide by number of subjects
stat = stat./length(subject);

%% Plot on SPM Template Brain using Fieldtrip
restoredefaultpath
addpath('/Volumes/Robert T5/fieldtrip');
ft_defaults

% Load some summy data
load('/Volumes/Robert T5/processed_Adult_MMN/2589/sourceall.mat');

load(['/Volumes/Robert T5/fieldtrip/template/'...
    'sourcemodel/standard_sourcemodel3d8mm.mat']);
template_grid = sourcemodel;
template_grid = ft_convert_units(template_grid,'mm');
clear sourcemodel;

% Load the SPM T1 brain
mri = ft_read_mri(['/Volumes/Robert T5/fieldtrip/'...
    'template/anatomy/single_subj_T1.nii']);

source = sourceall;
source.pos = template_grid.pos;

cd('/Volumes/Robert T5/HMM_test/spatial_results');

% For each state
for i = 1:10

    % Replace power with partial correlation value for the state
source.avg.pow(source.inside,:) = stat(:,i);

% Interpolate onto SPM brain
cfg              = [];
cfg.voxelcoord   = 'no';
cfg.parameter    = 'pow';
cfg.interpmethod = 'nearest';
sourceI  = ft_sourceinterpolate(cfg, source, mri);

% Write to .nii
disp('Saving results to .nii ...');
cfg = [];
cfg.filename = ['state' num2str(i)];
cfg.parameter = 'pow';
cfg.filetype = 'nifti';
cfg.scaling = 'no';
cfg.datatype = 'float';
ft_volumewrite(cfg, sourceI);

fff = ft_read_mri('state10.nii');

max(fff.anatomy(:))

% sourceI_spare = source;
% sourceI_spare.avg.pow = abs(sourceI_spare.avg.pow);
% 
% find_max_MNI(sourceI_spare, template_grid, 'yes');
% 
% cfg = [];
% cfg.funparameter = 'pow';
% ft_sourceplot(cfg,sourceI);
% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

% title(['State ',num2str(i)]);
% 
% disp('Producing 3D Plot');
% cfg = [];
% cfg.method         = 'surface';
% cfg.funparameter   = 'pow';
% cfg.projmethod     = 'nearest';
% cfg.surfinflated   = 'surface_inflated_both_caret.mat';
% %cfg.surfdownsample = 10
% cfg.projthresh     = 0.5;
% cfg.camlight       = 'no';
% ft_sourceplot(cfg, sourceI);
% ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
% colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
% %view ([70 0 50])
% %delete(findall(gcf,'Type','light'))
% light ('Position',[0 0 50])
% light ('Position',[0 -50 0])
% material dull;
% drawnow;
% view([0 0]);
% title(['State ',num2str(i)]);
% set(gca,'FontSize',14);
% drawnow;
% 

end

% Export to .gii file
addpath(genpath('/Users/rseymoue/Documents/GitHub/MQ_MEG_Scripts'));

d = dir(cd);

for i = 1:10
    nii_file = [cd '/' d(i+2).name];
    nii_to_workbench(nii_file);
end

%% Plot HMM data using OSL

restoredefaultpath
cd('/Volumes/Robert T5/osl/osl-core');
osl_startup

hmm.subj_inds = subj_inds;

[hmmstats, hmmstats_subj] = osl_hmm_stats_by_subject(hmm)

%% Make raincloud plots

addpath('/Users/rseymoue/Documents/scripts/RainCloudPlots-master/tutorial_matlab');
addpath(genpath('/Users/rseymoue/Documents/GitHub/MQ_MEG_Scripts'));
addpath('/Users/rseymoue/Documents/scripts/distinguishable_colors');

% Generate Colors
cl = distinguishable_colors(10,{'w'});

% Figure Position
fig_position = [200 200 1200 600]; % coordinates for figures

% Number of Occurences
ddd = [];
f4 = figure('Position', fig_position);
for i = 1:10
    ddd{1,1} = hmmstats_subj.nOccurrences(i,:);
    subplot(2,5,i);
    rm_raincloud(ddd,cl(i,:));
    xlim([-200 1000]);
    
    if i == 1  
        xlabel('nOccurrences');
        
    elseif i == 6
        xlabel('nOccurrences');
    end
    
    set(gca,'YTickLabel',[0:200:1000]);

    set(gca,'FontSize',14);
    
    set(gca,'YTickLabel',[]);
    
    title(['State ' num2str(i)]);

    
end

print('Noccurances','-dpng','-r300');


% Fractional Occupancy
ddd = [];
f4 = figure('Position', fig_position);
for i = 1:10
    ddd{1,1} = hmmstats_subj.FractionalOccupancy(i,:);
    subplot(2,5,i);
    rm_raincloud(ddd,cl(i,:));
    xlim([-0.1 0.5]);
    
    set(gca,'YTickLabel',[-0.1:0.2:0.5]);
    
    set(gca,'FontSize',14);
    
    set(gca,'YTickLabel',[]);
    
    if i == 1
        xlabel({'Fractional';'Occupancy'});
        
    elseif i == 6
        xlabel({'Fractional';'Occupancy'});
    end
    
    title(['State ' num2str(i)]);
    
    
end
print('fract_occupancy','-dpng','-r300');


% MeanLifeTime (s)
ddd = [];
f4 = figure('Position', fig_position);
for i = 1:10
    ddd{1,1} = hmmstats_subj.MeanLifeTime(i,:);
    subplot(2,5,i);
    rm_raincloud(ddd,cl(i,:));
    xlim([0 0.3]);
    
    set(gca,'YTickLabel',[0:0.05:0.3]);
    
    set(gca,'FontSize',14);
    
    set(gca,'YTickLabel',[]);
    
    if i == 1
        xlabel('Mean LifeTime (s)');
        
    elseif i == 6
        xlabel('Mean LifeTime (s)');
    end
    
    title(['State ' num2str(i)]);
    
    
end
print('mean_life_time','-dpng','-r300');

% INterval Length (s)
ddd = [];
f4 = figure('Position', fig_position);
for i = 1:10
    ddd{1,1} = hmmstats_subj.MeanIntervalLength(i,:);
    subplot(2,5,i);
    rm_raincloud(ddd,cl(i,:));
    xlim([0 9]);
    
    set(gca,'YTickLabel',[0:1:9]);
    
    set(gca,'FontSize',14);
    
    set(gca,'YTickLabel',[]);
    
    if i == 1
        xlabel('Interval Length (s)');
        
    elseif i == 6
        xlabel('Interval Length (s)');
    end
    
    title(['State ' num2str(i)]);
    
    
end
print('interval_length','-dpng','-r300');
