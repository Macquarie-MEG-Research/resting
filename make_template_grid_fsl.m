%% Make template grid from FSL brain

% NOTE: the path to the template file is user-specific
template = ft_read_mri('/Users/44737483/Documents/scripts_mcq/osl/std_masks/MNI152_T1_1mm.nii.gz');
template.coordsys = 'fsl'; % so that FieldTrip knows how to interpret the coordinate system
template = ft_convert_units(template,'cm');

ft_sourceplot(cfg,template)

% segment the template brain and construct a volume conduction model (i.e. head model): 
% this is needed to describe the boundary that define which dipole locations are 'inside' the brain.
cfg          = [];
template_seg = ft_volumesegment(cfg, template);

cfg          = [];
cfg.method   = 'singleshell';
template_headmodel = ft_prepare_headmodel(cfg, template_seg);
template_headmodel = ft_convert_units(template_headmodel, 'cm'); % Convert the vol to cm, because the CTF convenction is to express everything in cm.
 
% construct the dipole grid in the template brain coordinates
% the negative inwardshift means an outward shift of the brain surface for inside/outside detection
cfg = [];
cfg.grid.resolution = 0.8;
cfg.grid.tight      = 'yes';
cfg.inwardshift     = -1.5;
cfg.headmodel       = template_headmodel;
template_grid       = ft_prepare_sourcemodel(cfg);
 
% make a figure with the template head model and dipole grid
figure
hold on
ft_plot_vol(template_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));

% SAVE
cd('/Users/44737483/Documents/scripts_mcq/resting_state');
save template_grid template_grid


