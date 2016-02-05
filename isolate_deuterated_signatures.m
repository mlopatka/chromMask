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
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

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
    t = c.getRavel;
    
    %define skewed euclidean distance equation
    %euclidDistance = @(x,y)  sqrt(sum((x'-y').^2));
    l = reshape([anchors{:,2}],[],2);
    
    d1_search = 30; d2_search = 80;
    win_p(:,1) = l(:,1)-d1_search;
    win_p(:,2) = l(:,1)+d1_search;
    win_p(:,3) = l(:,2)-d2_search;
    win_p(:,4) = l(:,2)+d2_search;
    
    clearvars l d1_search d2_search
%     scatter(win_p(:,1), win_p(:,3), 'rp'); hold on; scatter(win_p(:,2), win_p(:,4), 'rp');
%     hold on; scatter(l(:,1), l(:,2), 'bp')
    
    %% define a distribution for distance from ideal anchor location
%     plot([0:1:100],wblpdf([0:1:100],25,1.05)) % lots of tail baaby!!!!

    %% define a distribution for the ideal anchor mz correlation score
%     plot([-1:0.01:2],(pdf('normal',[-1:0.01:2],1,0.25)))
    deut_coords = zeros([size(win_p,1),2]);
    
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