%% Set up path
addpath('./_classes')
addpath('./_helpers')
% Also ensure that the standard "Experiment" classes are added, or else
% the objects won't import properly

%% Load data and define basic parameters
expts = importdata('./experiment_objects.mat'); % Load in experiment objects here, can be just a single one if you want

n_recordings = numel(expts);
n_interp_bins = 180; % number of bins to interpolate each trial to

correct_thresh = 18;

%% Preprocess data for combination across recordings
[spks_cell, tc_cell, heading_cell, is_head_direction] = deal(cell(n_recordings, 1));
for rec = 1:n_recordings
	fprintf('Processing recording %d/%d\n', rec, n_recordings);
	hdp = HeadingDecoderPreprocessor(expts(rec));
	[spks_cell{rec}, tc_cell{rec}, heading_cell{rec}] = hdp.process('backward'); % 'forward' or 'backward' denotes whether to do light-on to light-off or vice versa
	is_head_direction{rec} = expts(rec).is_head_direction;
end

% Meaning across repeats for the recordings which have repeats
all_spks_aligned_m = cellfun(@(x) nanmean(x, 4), spks_cell, 'UniformOutput', false);
leave_one_tc_m = cellfun(@(x) nanmean(x, 4), tc_cell, 'UniformOutput', false);
heading_aligned_m = cellfun(@(x) nanmean(x, 3), heading_cell, 'UniformOutput', false);

%% Decode from the pseudopopulations
bin_centers = [-180 : 6 : 180];
bin_centers = bin_centers + 6/2;
bin_centers(end) = [];

th = linspace(-180, 180, n_interp_bins)';
tc = cat(1, leave_one_tc_m{:});
ts = cat(1, all_spks_aligned_m{:});
hd = HeadingDecoder(tc, ts, bin_centers);
hd.calculateHeadingDistribution();
hd.chooseHeading();
predicted_heading_pop = cat(2, hd.predicted_heading{:});

%% Calculate decoder performance metrics
de = [];
pc = [];
cm = [];
for ii = 1:numel(hd.predicted_heading)
	pc(ii, :) = calculatePercentCorrect(hd.predicted_heading{ii}, th, correct_thresh);
	de(ii, :) = calculateDecoderError(hd.predicted_heading{ii}, th);
	cm(:, :, ii) = getConfusionMatrix(hd.predicted_heading{ii}, th, 6);
end

%% Visualize results
mid_pt = length(predicted_heading_pop)/2;
f = figure;
set(gcf, 'Units', 'normalized', 'Position', [0.2542, 0.3650, 0.4914, 0.3500])
rep_th = repmat(th, [numel(hd.predicted_heading), 1]);
subplot(3, 1, 1)
plot(predicted_heading_pop, 'bx')
hold on
plot(rep_th, 'Color', [0.7, 0.7, 0.7])
title('predicted heading')
ylabel('heading')
xline(mid_pt, 'g:', 'LineWidth', 3);
xticks([mid_pt - n_interp_bins * 6: n_interp_bins:mid_pt + n_interp_bins * 6])
xticklabels([mid_pt - n_interp_bins * 6: n_interp_bins:mid_pt + n_interp_bins * 6] - mid_pt)
xlim([mid_pt - n_interp_bins * 6, mid_pt + n_interp_bins * 6]);

subplot(3, 1, 2)
scatter(1:size(de, 1) * size(de, 2), reshape(abs(de)', 1, []), 'filled', 'MarkerFaceAlpha', 0.3, 'MarkerFaceColor', [0.7, 0.7, 0.7]);
hold on
plot(n_interp_bins/2:n_interp_bins:length(predicted_heading_pop), nanmedian(abs(de), 2), 'ro:', 'LineWidth', 2)
ylabel('|decoder error|')
xlabel('frames from light off');
xline(mid_pt, 'g:', 'LineWidth', 3);
yline(90, '--')
xlim([mid_pt - n_interp_bins * 6, mid_pt + n_interp_bins * 6]);
xticks([mid_pt - n_interp_bins * 6: n_interp_bins:mid_pt + n_interp_bins * 6])
xticklabels([mid_pt - n_interp_bins * 6: n_interp_bins:mid_pt + n_interp_bins * 6] - mid_pt)
xlim([mid_pt - n_interp_bins * 6, mid_pt + n_interp_bins * 6]);

subplot(3, 1, 3)
plot(n_interp_bins/2:n_interp_bins:length(predicted_heading_pop), nanmean(pc, 2), 'mo:', 'LineWidth', 2)
ylim([0, 1])
xline(mid_pt, 'g:', 'LineWidth', 3);
ylabel('fraction correct')
yline(correct_thresh/180, '--')
xticks([mid_pt - n_interp_bins * 6: n_interp_bins:mid_pt + n_interp_bins * 6])
xticklabels([mid_pt - n_interp_bins * 6: n_interp_bins:mid_pt + n_interp_bins * 6] - mid_pt)
xlim([mid_pt - n_interp_bins * 6, mid_pt + n_interp_bins * 6]);
