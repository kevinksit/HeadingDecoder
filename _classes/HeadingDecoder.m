classdef HeadingDecoder < handle
	properties
		tuning_wts
		time_series

		memory = 5 % Samples looking back for fitting a line
		window = 0 % Samples around a single point for heading distribution

		n_segments  % number of identified segments in the data

		heading_distribution % likelihood across time that's used for ultimately determining heading
		predicted_heading % predicted heading at each time bin

		all_peak_centers
		fit_distribution
		target_idx
		scaling

		bin_centers
	end

	methods
		function obj = HeadingDecoder(tuning_curve, time_series, bin_centers)
			% Concatenating the whole tuning curve for normalization, so that we don't have each phase independently normalized
			obj.tuning_wts = tuning_curve; %rescaled_weights; %smoothdata(rescaled_weights, 2, 'sgolay'); %obj.getTuningWeights(rescaled_weights); %rescaled_weights; % obj.getTuningWeights(rescaled_weights);
			obj.time_series = time_series;
			obj.bin_centers = bin_centers;
			obj.n_segments = size(time_series, 3);
		end

		function decode(obj)
			obj.calculateHeadingDistribution();
			obj.chooseHeading();
		end

		function out = getHeadingTrace(obj, rescale_flag)
			if nargin < 2 || isempty(rescale_flag)
				rescale_flag = 0;
			end

			for ii = 1:length(obj.predicted_heading)
				if isnan(obj.predicted_heading(ii))
					out(ii) = NaN;
				else
					out(ii) = obj.bin_centers(obj.predicted_heading(ii));
				end
			end
		end

		function out = getHeadingDistribution(obj)
			if isempty(obj.heading_distribution)
				obj.calculateHeadingDistribution();
			end
			out = obj.heading_distribution;
		end
	end

	methods (Access = public)
		function calculateHeadingDistribution(obj)
			for s = 1:obj.n_segments
				% get the time series and tuning weights for each segment
				time_series = obj.getTs(s);
				tuning_wts = obj.getTw(s);
				tuning_wts = rescale(tuning_wts, 'InputMin', min(tuning_wts, [], 2), 'InputMax', max(tuning_wts, [], 2)); % rescale 0 - 1 by row, if you're getting an error here make sure you cloned GeneralHelperCode
				% Get rid of nans
				nan_row = any(isnan(tuning_wts)') | any(isnan(time_series)');
				tuning_wts(nan_row, :) = [];
				time_series(nan_row, :) = [];
				% Calculate heading distribution
				heading_distribution = (time_series' * tuning_wts)./sum(time_series' * tuning_wts);
				% Temporally filter
				out = obj.temporalFilter(heading_distribution', 6)';
				heading_distribution_all{s} = out;
			end

			% Assign heading distribution
			obj.heading_distribution = heading_distribution_all;
		end

		function [tw] = getTw(obj, s)
			tw = obj.tuning_wts(:, :, s);
		end

		function [ts] = getTs(obj, s)
			ts = obj.time_series(:, :, s);
		end

		function [filt_decoded_probabilities] = temporalFilter(obj, decoded_probabilities, window_size)
			%BAYESIAN_TEMPORAL_FILTER Summary of this function goes here
			%   Detailed explanation goes here
			% keyboard
			% 	ca_time = 0.1:0.1:size(decoded_probabilities, 2);
			% 	Fs = 1/mode(diff(ca_time));
			% window_size = round(window_size*Fs);
			half_win = round(window_size/2);

			%% Pad matrix with zeros
			zero_pad = zeros(size(decoded_probabilities,1),window_size);
			padded_decoded_probabilities = [zero_pad decoded_probabilities zero_pad];

			for step_i = 1:size(decoded_probabilities,2)
				current_window = padded_decoded_probabilities(:,step_i+half_win: step_i-1+window_size+half_win);
				filt_decoded_probabilities(:,step_i) = expm1(sum(log1p(current_window),2,'omitnan'));
				filt_decoded_probabilities(:,step_i) = filt_decoded_probabilities(:,step_i)./sum(filt_decoded_probabilities(:,step_i));
			end

		end

		function chooseHeading(obj, choose_method)
			if nargin < 2 || isempty(choose_method)
				choose_method = 'max';
			end

			switch choose_method
				case 'max'
					for s = 1:obj.n_segments
						[~, predicted_idx] = max(obj.heading_distribution{s}, [], 2);
						predicted_heading{s} = obj.bin_centers(predicted_idx);
					end
					obj.predicted_heading = predicted_heading;
				case 'predict'
					for s = 1:obj.n_segments
						heading_distribution = obj.heading_distribution{s};
						predicted_heading = [];
						for ii = 1:size(heading_distribution, 1)
							candidates = islocalmax(rescale(movmean(heading_distribution(ii, :), 5)), 'MinProminence', 0.5);
							peak_centers = find(candidates);
							if ~isempty(peak_centers)
								if ii <= obj.memory
									predicted_heading(ii) = min(peak_centers);
								else
									x = (ii - obj.memory) : (ii - 1);
									prediction = obj.predictNextPeak(x, predicted_heading);
									if numel(peak_centers) > 1
										target_idx = obj.minCircularDistance(peak_centers, prediction);% min(abs(peak_centers - prediction)); % account for wrapping
									else
										target_idx = 1;
									end
									predicted_heading(ii) = peak_centers(target_idx);
									obj.target_idx(ii) = target_idx;
									temp(ii) = prediction;
								end

							else
								predicted_heading(ii) = NaN;
							end
							obj.all_peak_centers{ii} = peak_centers;
						end
						predicted_heading_all{s} = nan(1, length(predicted_heading));
						predicted_heading_all{s}(~isnan(predicted_heading)) = obj.bin_centers(predicted_heading(~isnan(predicted_heading)));
					end
					obj.predicted_heading = predicted_heading_all;
			end    
		end

	end
	methods % These are some old methods that used to be used with predicting heading, but are now deprecated
		function prediction = predictNextPeak(obj, x, predicted_heading)
			peak_fit = polyfit(x, predicted_heading(x), 1);
			if peak_fit(:, 1) >= -0.1 && peak_fit(:, 1) <= 0.1
				prediction = mean(predicted_heading(x));
			else
				prediction = polyval(peak_fit, x(end) + 1);
			end
			% wrap
			if prediction > length(obj.bin_centers)
				prediction = prediction - length(obj.bin_centers);
			elseif prediction < 0
				prediction = prediction + length(obj.bin_centers);
			end
		end

		function peak_centers = getPeakCenters(obj, candidate_peaks)
			% Find rises and falls (because peaks can be consecutive values above threshold, we don't want to count each point as an individual peak)
			change = diff(candidate_peaks(:));
			peak_locs = find(candidate_peaks);
			split = find([0 (diff(peak_locs)) > 1]);
			n_peaks = numel(split) + 1;
			if n_peaks > 1
				for p = 1:n_peaks
					if p == 1
						peak_group{p} = peak_locs(1:split(p) - 1);
					elseif p == numel(split) + 1
						peak_group{p} = peak_locs(split(p - 1):end);
					else
						peak_group{p} = peak_locs(split(p - 1) : split(p) - 1);
					end
					peak_centers = cell2mat(vertcat(cellfun(@(x) mean(x(:)), peak_group, 'UniformOutput', false)));
				end
			else
				peak_centers = mean(peak_locs);
			end
		end

		function idx = minCircularDistance(obj, x, y)
			dist = y - x + 1;
			dist2 = x - 1 + size(obj.heading_distribution, 2) - y;
			for ii = 1:length(dist) % for each possible peak, choose the smaller
				candidates(ii) = min(abs([dist(ii), dist2(ii)]));
			end
			[~, idx] = min(candidates);
		end
	end
end
