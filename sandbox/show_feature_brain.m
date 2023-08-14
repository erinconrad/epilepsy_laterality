%% Get file locs
locations = fc_toolbox_locs;
results_folder = [locations.main_folder,'results/'];
inter_folder = [results_folder,'analysis/new_outcome/data/'];
atlas_folder = [locations.main_folder,'data/atlas/'];

% add script folder to path
scripts_folder = locations.script_folder;
addpath(genpath(scripts_folder));

%% Load data file
data = load([inter_folder,'main_out.mat']);
data = data.out;

%% Load MNI brain file
mni_brain = load([atlas_folder,'ICBM152.mat']);
vertices = mni_brain.vertices;
faces = mni_brain.faces;

%% Get variables
all_locs = data.all_locs;
all_names = data.all_names;
all_labels = data.all_labels;
all_spike_rates = data.all_spikes;
npts = length(all_names);

%% Show feature on a brain
all_feature = all_spike_rates;



for p = 1:npts

locs = all_locs{p};
labels = all_labels{p};
feature = all_feature{p};

h = trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3));
h.LineStyle = 'none';
h.FaceAlpha = 0.1;
h.FaceColor = [0.7 0.7 0.7];    
hold on

a1 = scatter3(locs(:,1),locs(:,2),locs(:,3),60,'k','LineWidth',2);
a2 = scatter3(locs(:,1),locs(:,2),locs(:,3),50,feature,'filled');
if ~all(isnan(feature))
    clim([min(feature) max(feature)])
end
%text(locs(:,1),locs(:,2),locs(:,3),labels,'fontsize',15);
view(-180,33)
hold off
pause
end


