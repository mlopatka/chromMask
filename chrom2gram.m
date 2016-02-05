classdef chrom2gram
    % a two-dimensional comprehensive chromatogram object class for Matlab
    properties
        d1 % first chromatographic dimension in data points
        d2 % second chromatographic dimension in data points
        mass_range_max % maximum m/z value observed
        mass_range_min % minimum m/z observed
        masses % the vector of niominal m/z values in the MS channel
        scan_acquisition_time % retention time 
        scan_duration % the duration of one scanning event in seconds
        scan_index % scan index of each point
        time_idx % time index of each data point
        xyz % the raw chromatographic data
        tile_idx % for the tile representation of the chrom2gram coordinates for 2-d representation
        p % probabilistic peak detection vector or matrix
        baseline % baseline from row-wise baseline correction
        peakAreas % list of peaks
        rt_idx % corresponding retention times for peak areas 
        anchor_pts % any kind of anchor points available (for example deuterated compunds
    end
    
    methods
        %% constructor
        function obj = chrom2gram(c) % object constructor from struct
            c_idx = fieldnames(c);  
            for i = 1:numel(c_idx)
                if strcmpi(c_idx{i}, 'rt_d2') % special conditional for legacy format
                    obj.rt_idx = [c.rt_d1,c.rt_d2];
                elseif strcmpi(c_idx{i}, 'chrom2gram') 
                    obj.xyz = c.chrom2gram;
                else
                    try
                        eval(['obj.', c_idx{i}, ' = c.', c_idx{i}, ';']);
                    catch
                        disp(['discarding ', c_idx{i}, ' field to adhere to chrom2gram schema']);
                    end
                end
            end
            obj.tile_idx = generateTiles(obj);
        end 
        %% return the total intensity (summed along masses)
        function total_intensity = getTotalIntensity(obj)
            total_intensity = (sum(obj.xyz,2));
        end
        %% reshape the data to a 2 dimensional tensor
        function raveled = getRavel(obj)
            if and(size(obj.xyz,1)== obj.d1*obj.d2, size(obj.xyz,2)>1)               
                raveled = reshape(obj.xyz, [obj.d2, obj.d1, numel(obj.masses)]);
            else
                raveled = obj.xyz;
            end
        end
        %% overloading the imagesc function to get nice plots
        function imagesc(obj)% we can overload the plot operator
            figure; imagesc(sum(obj.getRavel,3)); axis xy
        end
        %% issolate specific ion channel(s)
        function [ion_layers,ion_mask] = getIonChannel2D(obj,ion_mask)
            if ~isequal(intersect(obj.masses, ion_mask), ion_mask)
                error('requested masses are not all present in the chromatogram')
            end
            t = obj.getRavel;
            ion_layers = t([1:size(t,1)], [1:size(t,2)], [ion_mask]);
        end
        %% convert raveled representation tile indeces to xyz tile indeces
        function vecTilesIDX = getxyzTiles(obj)
            vecTilesIDX = cell(1,size(obj.tile_idx,1));
            for vert = 1:size(obj.tile_idx,1) % maybe get rid of this loop somehow
                [X,Y] = meshgrid(obj.tile_idx(vert,1):obj.tile_idx(vert,2),...
                    obj.tile_idx(vert,3):obj.tile_idx(vert,4));
                temp = ((X.*obj.d2)+(Y))-obj.d2;
                vecTilesIDX{vert} = temp(:);
            end
        end
        %% generate 2-D tile indeces based ont he hardcoded parameters inside this function.
        function tile_idx = generateTiles(obj)
            m = obj.d1; n = obj.d2; w = 40; h = 40; ov = 20; oh = 20;
            xind = round(linspace(1, m-w, ceil((m-w)/(w-oh))));
            yind = round(linspace(1, n-h, ceil((n-h)/(h-ov))));
            [q2,q] = meshgrid(xind, yind);
            clearvars xind yind
            pairs = [q2(:) q(:)];
            clearvars q2 q
            tile_idx(:,1) = pairs(:,1);
            tile_idx(:,2) = pairs(:,1) + w;
            tile_idx(:,3) = pairs(:,2);
            tile_idx(:,4) = pairs(:,2) + h; 
        end
        %% go to the local summed tile representation
        function kalForm = getLiscaForm(obj)
            kalForm = zeros(size(obj.tile_idx,1), numel(obj.masses));
            if ~isempty(obj.masses)
                temp_idx = obj.getxyzTiles;
                for iterable = 1:size(kalForm,1)
                    kalForm(iterable,:) = sum(obj.xyz(temp_idx{iterable},:));    
                end
            end
        end
        %% suggest an ion mask based on low standard deviation in a particular mass channel
        function ion_mask = getIonMask(obj)
            if isempty(obj.masses)
                error('you can not perform on filtering on a first order signal!');
            end
            std_vals = (std(double(obj.xyz)));
            ion_mask = (std_vals-mean(std_vals(200:end)))<10e-10;
        end
        %% plot ion-masked tic from GCxGC-MS data
        function h = showIonMaskedTIC(obj)
            t = obj.xyz(:,~obj.getIonMask);
            t = reshape(t, obj.d2,obj.d1,size(t,2));
%             h = imagesc(sum(t,3),[0, 10^(5)]); axis xy;
            h = pcolor(sum(t,3)); shading interp; axis xy;
        end
    end % end of methods for this object type
end % end of class definition