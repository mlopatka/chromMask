function deut_coords = isolate_deuterated_signatures(c)
%% function to isolate the location of 14 deuterated alcanes and aromatics 
% added to fire debris samples for GCxGC-MS analysis.
% example useage: deut_coords = isolate_deuterated_signatures(c)
% - where c is a chrom2gram object type
% - where deut_coords is a 14x2 matrix with d1 and d2 indeces locating the
% deuterated markers.
%
% hardcoded parameters include: 
    % search window size
    % distance distribution,
    % mass correlation distribution.
    % lambda param weigting emphasis of distance versus mass similarity
% 
% Inspired by the work published by Woldegebriel & Truyols.
% adapted from the original manuscript to use manual distribution 
% definitions (please cite if using this for academic research).
% M.T. Woldegebriel & G. Vivo Truyols (2015).
% Probabilistic Model for Untargeted Peak Detection in LC-MS Using Bayesian Statistics. 
% Analytical Chemistry, 87 (14), 7345-7355. doi: 10.1021/acs.analchem.5b01521
% 
% Implementation by Martin Lopatka [martin.lopatka@gmail.com]
% Copyright (c) 2016, Martin Lopatka
% All rights reserved.

%% original code to get the ion signatures using tight lisca window on standard ladder analysis 
%     d1_smudge = 2; d2_smudge = 4;
%     % scale up to multidimensional without a loop
%     win_p(:,1) = l(:,1)-d1_smudge;
%     win_p(:,2) = l(:,1)+d1_smudge;
%     win_p(:,3) = l(:,2)-d2_smudge;
%     win_p(:,4) = l(:,2)+d2_smudge;
    %scatter(win_p(:,1), win_p(:,3), 'rp'); hold on; scatter(win_p(:,2), win_p(:,4), 'rp');
    %hold on; scatter(l(:,1), l(:,2), 'bp') 
%     for i = 1:size(win_p,1)
%         ion_anchors(i,:) = (squeeze(sum(sum(squeeze(t([win_p(i,3):win_p(i,4)],[win_p(i,1):win_p(i,2)], [1:size(t,3)])),2),1)));
%     end
%    % sum to unity
%    ion_anchors = ion_anchors./repmat(sum(ion_anchors,2),[1, size(ion_anchors,2)]);
%    save('anchors.mat', {'l', 'ion_anchors', 'D'});
%% now we just load the saved anchor profiles
    load anchors.mat 
    % These should be 14 by 4 cell array with the following fields:
    % {'compound name',[locations],[ions],[reference spectra]}
    %TODO replace with cleaner spectra from GC-MS data    
    
    l = reshape([anchors{:,2}],[],2);
    d1_search = 30; d2_search = 80;
    win_p(:,1) = l(:,1)-d1_search;
    win_p(:,2) = l(:,1)+d1_search;
        
    if isa(c,'chrom2gram')
        t = c.getRavel;
        win_p(:,3) = l(:,2)-d2_search;
        win_p(:,4) = l(:,2)+d2_search;
        deut_coords = zeros([size(win_p,1),2]);
    else % must be GC-MS data
        deut_coords = zeros([size(win_p,1),1]);
    end
    
    clearvars l d1_search d2_search

    
    %% begin search loop if 1D GC-MS data
    if ~isa(c,'chrom2gram')
        for i = 1:size(win_p,1)
            search_window = c(win_p(i,1):win_p(i,2), anchors{i,3}); % 1D search window
            x_idx_rt = size(search_window,1);
            y_p = zeros([x_idx_rt, 2]); % make sure no contaminations from last iteration
            [x,y] = find(sum(search_window,2)>0); % this is HELLA sloppy , could be faster, but im tired now... optimize later with arithmetic.
            dist_score = pdist2([x,y],[round(x_idx_rt/2),1], 'euclidean'); % ewuclidean distance by default.
            y_p(:,1) = (wblpdf(dist_score,25,1.05));
            y_p(dist_score==0,1) = max(y_p(:,1)); % override little cliff in the wbl pdf at 0.
            y_p(or(isinf(dist_score), isnan(dist_score)),1) = 0;
            y_p(:,1) = y_p(:,1)./max(y_p(:,1)); 
            % likelihood of that distance from centroid

            search_window = search_window./repmat(sum(search_window,2),[1, size(search_window,2)]);
            corr_score = pdist2(anchors{i,4}(anchors{i,3}),search_window, 'cosine'); % must match in dimensionality with the masked t variable
            %corr_score = pdist2(ion_anchors(i,:),search_window, 'seuclidean', sum(abs(D))); % this is interesting but requires a different distribution
            y_p(:,2) = (pdf('normal',corr_score,0,0.1));
            y_p(:,2) = y_p(:,2)./max(y_p(:,2));
            %likelihood of that correlation in the mass channel 
            lambda = 0.5; % weight the mass correlation more than the distance
            [~, idx] = max((y_p(:,1)*lambda).*y_p(:,2)); % isolated most likely location of deuterated compound apex
            deut_coords(i) = win_p(i,1)+idx;
        end
    %% begin search loop for 2D comprehensive chromatogram    
    else
        for i = 1:size(win_p,1)
            search_window = t([win_p(i,3):win_p(i,4)],[win_p(i,1):win_p(i,2)], anchors{i,3});
            [x_idx_w, y_idx_w, z_idx_w] = size(search_window);
            y_p = zeros([x_idx_w * y_idx_w, 2]); % make sure no contaminations from last iteration
            [x,y] = find(sum(search_window,3)>0); % this is HELLA sloppy , could be faster, but im tired now... optimize later with arithmetic.
            dist_score = pdist2([x,y],[round(x_idx_w/2),round(y_idx_w/2)], 'euclidean'); % ewuclidean distance by default.
            y_p(:,1) = (wblpdf(dist_score,25,1.05));
            y_p(dist_score==0,1) = max(y_p(:,1)); % override little cliff in the wbl pdf at 0.
            y_p(or(isinf(dist_score), isnan(dist_score)),1) = 0;
            y_p(:,1) = y_p(:,1)./max(y_p(:,1)); 
            % likelihood of that distance from centroid

            search_window= double(reshape(search_window, [x_idx_w*y_idx_w,z_idx_w]));
            search_window = search_window./repmat(sum(search_window,2),[1, size(search_window,2)]);
            corr_score = pdist2(anchors{i,4}(anchors{i,3}),search_window, 'cosine'); % must match in dimensionality with the masked t variable
            %corr_score = pdist2(ion_anchors(i,:),search_window, 'seuclidean', sum(abs(D))); % this is interesting but requires a different distribution
            y_p(:,2) = (pdf('normal',corr_score,0,0.1));
            y_p(:,2) = y_p(:,2)./max(y_p(:,2));
            %likelihood of that correlation in the mass channel 
            lambda = 0.5; % weight the mass correlation more than the distance
            [~, idx] = max((y_p(:,1)*lambda).*y_p(:,2)); % isolated most likely location of deuterated compound apex
            [win_idx_2, win_idx_1] = ind2sub([x_idx_w,y_idx_w],idx);
            deut_coords(i,1) = win_p(i,1)+win_idx_1;
            deut_coords(i,2) = win_p(i,3)+win_idx_2;
        end
    end
       %% define a distribution for distance from ideal anchor location 2-D
%     plot([0:1:100],wblpdf([0:1:100],25,1.05)) % lots of tail baaby!!!!

    %% define a distribution for the ideal anchor mz correlation score 2-D
%     plot([-1:0.01:2],(pdf('normal',[-1:0.01:2],1,0.25)))
