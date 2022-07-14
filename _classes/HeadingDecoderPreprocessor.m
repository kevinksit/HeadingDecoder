classdef HeadingDecoderPreprocessor < handle
    properties
        expt
        is_bad_recording = false;
        
        n_interp_bins
        n_cells
        window
    end
    
    methods
        function obj = HeadingDecoderPreprocessor(expt)
            if ~isa(expt, 'Cue') || ~isa(expt, 'Control')
                error('Input needs to be an Experiment object')
            end
            obj.expt = expt; % this needs to be an "Experiment" object
            obj.n_interp_bins = 180;
            obj.window = 15;
            obj.n_cells = numel(obj.expt.is_head_direction);
        end
        
        function [spks_aligned, leave_one_tc, heading_aligned] = abort(obj)
            obj.is_bad_recording = true;
            spks_aligned = nan(obj.n_cells, obj.n_interp_bins, 2*obj.window);
            leave_one_tc = nan(obj.n_cells, 360/obj.expt.bin_width, 2*obj.window);
            heading_aligned = nan(obj.n_interp_bins, 2*obj.window);
        end
        
        function [spks_aligned, leave_one_tc, heading_aligned] = process(obj, direction)
            % Direct controls if it's going to be light-on to light-off or
            % light-off to light-on
            switch direction
                case 'forward'
                    last = numel(obj.expt.light_data);
                case 'backward'
                    last = numel(obj.expt.light_data) - 1;
                otherwise
                    error('Did not provide a usable direction: ''forward'' or ''backward''');
            end
            
            if last == 0 % single repeat? no dark -> light
               [spks_aligned, leave_one_tc, heading_aligned] = obj.abort();
               return;
            end
            
            leave_one_tc = nan(obj.n_cells, 360/obj.expt.bin_width, 2*obj.window, last);
            heading_aligned = nan(obj.n_interp_bins, 2*obj.window, last);
            spks_aligned = nan(obj.n_cells, obj.n_interp_bins, 2*obj.window, last);
            for ii = 1:last % should go through each rep
                switch direction
                    case 'forward'
                        light_data = obj.expt.light_data(ii);
                        dark_data = obj.expt.dark_data(ii);
                    case 'backward'
                        light_data = obj.expt.light_data(ii+1);
                        dark_data = obj.expt.dark_data(ii);
                end
                
                [light_spks, dark_spks] = deal(nan(obj.n_cells, obj.n_interp_bins, obj.window));
                
                [heading_final_l, heading_final_d] = deal(nan(obj.n_interp_bins, obj.window));
                
                [tuning_curve_l, tuning_curve_d] = deal(nan(obj.n_cells, 360/obj.expt.bin_width, obj.window));
                heading_l = light_data.get('heading');
                heading_d = dark_data.get('heading');
                
                [heading_reference_l, heading_reference_d] = obj.findReferenceHeading(heading_l, heading_d, direction);
                
                % For light
                spks = light_data.get('spikes');
                heading_l = obj.cleanHeading(heading_l);
                segment_light = obj.segmentTrials(heading_l, heading_reference_l);
                if ~isempty(segment_light)
                    [spks_resampled_light, heading_resampled_light] = obj.resampleSpikes(spks, heading_l, segment_light);
                else
                    [spks_aligned, leave_one_tc, heading_aligned] = obj.abort();
                    return;
                end
                
                % For dark
                spks = dark_data.get('spikes');
                
                heading_d = obj.cleanHeading(heading_d);
                segment_dark = obj.segmentTrials(heading_d, heading_reference_d);
                if ~isempty(segment_dark)
                    [spks_resampled_dark, heading_resampled_dark] = obj.resampleSpikes(spks, heading_d, segment_dark);
                else
                    [spks_aligned, leave_one_tc, heading_aligned] = obj.abort();
                    return;
                end

                % lift binData from the experiments, yeah
                light_tuning = nan(obj.n_cells, 60, numel(segment_light));
                for s = 1:numel(segment_light)
                    light_tuning(:, :, s) = obj.expt.binData(spks_resampled_light(:, :, s), heading_resampled_light(:, s));
                end
                
                dark_tuning = nan(obj.n_cells, 60, numel(segment_dark));
                for s = 1:numel(segment_dark)
                    dark_tuning(:, :, s) = obj.expt.binData(spks_resampled_dark(:, :, s), heading_resampled_dark(:, s));
                end
                
                [light_tc, dark_tc] = obj.getTuningCurves(light_tuning, dark_tuning);
                    
                % insert properly
                switch direction
                    case 'forward'
                        if ~isempty(segment_light) % not super pretty, think about a better way here...
                            light_spks(:, :, end - size(spks_resampled_light, 3)+1:end) = spks_resampled_light;
                            heading_final_l(:, end - size(heading_resampled_light, 2) + 1 : end) = heading_resampled_light;
                            tuning_curve_l(:, :, end - size(light_tc, 3) + 1 : end) = light_tc;
                        end
                        
                        if ~isempty(segment_dark)
                            dark_spks(:, :, 1:size(spks_resampled_dark, 3)) = spks_resampled_dark;
                            heading_final_d(:, 1:size(heading_resampled_dark, 2)) = heading_resampled_dark;
                            tuning_curve_d(:, :, 1:size(dark_tc, 3)) = dark_tc;
                        end                        
                        spks_aligned(:, :, :, ii) = cat(3, light_spks, dark_spks);
                        heading_aligned(:, :, ii) = cat(2, heading_final_l, heading_final_d);
                        leave_one_tc(:, :, :, ii) = cat(3, tuning_curve_l, tuning_curve_d);
                        
                    case 'backward'
                        if ~isempty(segment_light)
                            light_spks(:, :, 1:size(spks_resampled_light, 3)) = spks_resampled_light;
                            heading_final_l(:, 1:size(heading_resampled_light, 2)) = heading_resampled_light;
                            tuning_curve_l(:, :, 1:size(light_tc, 3)) = light_tc;
                        end
                        if ~isempty(segment_dark)
                            dark_spks(:, :, end - size(spks_resampled_dark, 3) + 1 : end) = spks_resampled_dark;
                            heading_final_d(:, end - size(heading_resampled_dark, 2) + 1 : end) = heading_resampled_dark;
                            tuning_curve_d(:, :, end - size(dark_tc, 3) + 1 : end) = dark_tc;
                        end
                        
                        fprintf("Light trials; %d | Dark trials: %d\n", size(light_tc, 3), size(dark_tc, 3));
                        spks_aligned(:, :, :, ii) = cat(3, dark_spks, light_spks);
                        heading_aligned(:, :, ii) = cat(2, heading_final_d, heading_final_l);
                        leave_one_tc(:, :, :, ii) = cat(3, tuning_curve_d, tuning_curve_l);
                end
            end
        end
        
        
        function [heading_reference_l, heading_reference_d] = findReferenceHeading(obj, heading_light, heading_dark, direction)
            switch direction
                case 'forward'
                    front_heading = heading_light;
                    back_heading = heading_dark;
                case 'backward'
                    front_heading = heading_dark;
                    back_heading = heading_light;
            end
            slow_idx = find(diff(front_heading(end-100:end)) < 1); % find times when the speed has decreased, this is because the last trial isn't always exactly right
            slow_idx = slow_idx(diff(slow_idx) == 1); % consecutive only (tsries to get rid of issues)
            front_reference = front_heading(end-(100-slow_idx(end))); % last slow trial is ideal, because it will get rid of junk
            
            back_reference = back_heading(1);
            
            switch direction
                case 'forward'
                    heading_reference_l = front_reference;
                    heading_reference_d = back_reference;
                case 'backward'
                    heading_reference_l = back_reference;
                    heading_reference_d = front_reference;
            end
        end
        
        function heading = cleanHeading(obj, heading)
            is_jump = abs(diff(heading)) > 20 & abs(diff(heading)) < 350;
            ct = 1;
            while any(is_jump)
                if ct == 10
                    break
                end
                heading(is_jump) = heading(max(find(is_jump) - 1, 1));
                is_jump = abs(diff(heading)) > 20 & abs(diff(heading)) < 350;
                ct = ct + 1;
            end
        end
        
        function segment = segmentTrials(obj, heading, heading_reference)
            threshold = 1;
            breaks = find(abs(heading - heading_reference)<threshold);
            segment_length = diff(breaks);
            segment_length(segment_length == 1) = [];
            std_sl = std(segment_length);
            while any(diff(breaks) > 200) || std_sl > 20
                if threshold > 5
                    breaks = 0;
                    disp('o no')
                    break
                end
                threshold = threshold + 0.5;
                breaks = find(abs(heading - heading_reference)<threshold);
                segment_length = diff(breaks);
                segment_length(segment_length == 1) = [];
                std_sl = std(segment_length);
            end
            
            segment = [];
            for ii = 1:length(breaks) - 1
                segment{ii} = breaks(ii) + 1:breaks(ii+1);
            end
            if isempty(segment) % no segments
                return; % will error on the next line
            end
            segment(cellfun(@length, segment) < 10) = []; % get rid of small chunks
        end
        
        function [spks_resampled, heading_resampled] = resampleSpikes(obj, spks, heading, segment)
            if isempty(segment)
                disp('Something not good right now...')
            else
                % circshift the data
                spks_resampled = nan(size(spks, 1), obj.n_interp_bins, numel(segment));
                heading_resampled = nan(obj.n_interp_bins, numel(segment));
                for ii = 1:numel(segment)
                    % differences
                    [~, min_idx] = min(abs(heading(segment{ii}) - (-180))); % shift the data so that it "starts" at 0 
                    shift = min_idx;
                    
                    spks_resampled(:, :, ii) = resample(circshift(spks(:, segment{ii}), 1 - shift, 2), obj.n_interp_bins, length(segment{ii}), 0, 'Dimension', 2);
                    heading_resampled(:, ii) = resample(circshift(heading(segment{ii}), 1 - shift), obj.n_interp_bins, length(segment{ii}), 0);
                end
            end
        end
        
        function [light_tc, dark_tc] = getTuningCurves(obj, light_tuning, dark_tuning, style)
            if nargin < 4 || isempty(style)
                style = 'local';
            end
            
            light_tc = nan(size(light_tuning));
            dark_tc = nan(size(dark_tuning));
            switch style
                case 'local'
                    for ii = 1:size(light_tuning, 3)
                        light_tc(:, :, ii) = nanmean(light_tuning(:, :, 1:end~=ii), 3); % temp use all
                    end
                    for ii = 1:size(dark_tuning, 3)
                        dark_tc(:, :, ii) = nanmean(dark_tuning(:, :, 1:end~=ii), 3);
                    end
                case 'global'
                    for ii = 1:size(light_tuning, 3)
                        light_tc(:, :, ii) = nanmean(cat(3, light_tuning(:, :, 1:end~=ii), dark_tuning), 3);
                    end
                    
                    for ii = 1:size(dark_tuning, 3)
                        dark_tc(:, :, ii) = nanmean(cat(3, light_tuning, dark_tuning(:, :, 1:end~=ii)), 3);
                    end
            end
        end
    end
end
