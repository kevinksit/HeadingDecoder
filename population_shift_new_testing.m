%% Load data and define basic parameters
expts = importdata('./experiment_objects.mat'); % Load in experiment objects here

n_recordings = numel(expts);
n_interp_bins = 180; % number of bins to interpolate each trial to

is_bad_recording = false(1, n_recordings);
all_spks_aligned = cell(1, n_recordings);
leave_one_tc = cell(1, n_recordings);
heading_aligned = cell(1, n_recordings);
is_head_direction = cell(1, n_recordings);

%% Preprocess data for combination across recordings
for rec = 1:n_recordings
	fprintf('Processing recording %d/%d\n', rec, n_recordings);
    hdp = HeadingDecoderPreprocessor(expts(rec));
    [spks{rec}, tc_cell{rec}, head{rec}] = hdp.process('forward');
    is_bad_recording(rec) = hdp.is_bad_recording;
    is_head_direction{rec} = expts(rec).is_head_direction;
end

all_spks_aligned_m = cellfun(@(x) nanmean(x, 4), spks, 'UniformOutput', false);
leave_one_tc_m = cellfun(@(x) nanmean(x, 4), tc_cell, 'UniformOutput', false);
heading_aligned_m = cellfun(@(x) nanmean(x, 3), head, 'UniformOutput', false);

%% Load cell clusters from data

if ~exist('cluster_id')
    load('clustering_snapshot_05Jan2022.mat');
end

%% get rid of "bad" recordings
n_head_direction_cells = cellfun(@sum, is_head_direction);
is_usable_recording = (n_head_direction_cells > 0) & (~is_bad_recording); % tweak this parameter

%% split cell cluster id's into individual recordings (important for later)
num_cells = cellfun(@length, is_head_direction);
tally = cumsum(num_cells);
tally = cat(2, 0, tally);
for ii = 1:numel(tally) - 1
    cluster_id_cell{ii} = cluster_id(tally(ii) + 1 : tally(ii + 1));
end

%% decode time baybeee
bin_centers = [-180 : 6 : 180];
bin_centers = bin_centers + 6/2;
bin_centers(end) = [];

correct_thresh = 18;
u_clusters = unique(cluster_id);

save_flag = false;
controlled_rand = false;
if controlled_rand
    cell_idx_all = importdata('E:\_HeadingProjectFigures\Figure6\cell_idx_all.mat');
    fns = fieldnames(cell_idx_all);
end

clust_tab = tabulate(cluster_id);
sample_size = min(cat(1, 200, clust_tab(:, 2)));

for c = u_clusters
    
    th = linspace(-180, 180, n_interp_bins)';
    tc = cat(1, leave_one_tc_m{is_usable_recording});
    ts = cat(1, all_spks_aligned_m{is_usable_recording});
    if c == -1
        is_cluster = cat(2, cluster_id_cell{is_usable_recording}) ~= c;
        tc = tc(is_cluster, :, :);
        ts = ts(is_cluster, :, :);
    else
        is_cluster = cat(2, cluster_id_cell{is_usable_recording}) == c;
        tc = tc(is_cluster, :, :);
        ts = ts(is_cluster, :, :);
        if controlled_rand
            cell_idx = cell_idx_all.(fns{c});
        else
            cell_idx = randperm(size(tc, 1), round(1 * size(tc,1)));
        end
        tc = tc(cell_idx(1:min(sample_size, length(cell_idx))), :, :);
        ts = ts(cell_idx(1:min(sample_size, length(cell_idx))), :, :);
    end
    
    hd = HeadingDecoder(tc, ts, bin_centers);
    
    hd.calculateHeadingDistribution();
    hd.chooseHeading();
    predicted_heading_pop = cat(2, hd.predicted_heading{:});
    de = [];
    pc = [];
    cm = [];
    for ii = 1:numel(hd.predicted_heading)
        pc(ii, :) = calculatePercentCorrect(hd.predicted_heading{ii}, th, correct_thresh);
        de(ii, :) = calculateDecoderError(hd.predicted_heading{ii}, th);
        cm(:, :, ii) = getConfusionMatrix(hd.predicted_heading{ii}, th, 6);
    end
    
    %% visualize
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
    prettyPlot
    switch c
        case -1
            output_struct.all.de = de;
            output_struct.all.pc = pc;
            output_struct.all.th = th;
            output_struct.all.predicted_heading = predicted_heading_pop;
        case 1
            output_struct.heading.de = de;
            output_struct.heading.pc = pc;
            output_struct.heading.th = th;
            output_struct.heading.predicted_heading = predicted_heading_pop;
            cell_idx_all.heading = cell_idx;
        case 2
            output_struct.visual.de = de;
            output_struct.visual.pc = pc;
            output_struct.visual.th = th;
            output_struct.visual.predicted_heading = predicted_heading_pop;
            cell_idx_all.visual = cell_idx;
        case 3
            output_struct.multimodal.de = de;
            output_struct.multimodal.pc = pc;
            output_struct.multimodal.th = th;
            output_struct.multimodal.predicted_heading = predicted_heading_pop;
            cell_idx_all.multimodal = cell_idx;
    end
    if save_flag
        save(sprintf('%s_off2on_outputs.mat', date), 'output_struct')
    end
    pause
end
