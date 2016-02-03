function pmat = maskedPeakDetect(chrom)
% hardcoding some essential things that should stay relatively constant
% over a batch processing paradigm, so I dont want to have to pass them in
% and out fo the function.
addpath '../../lisca/src/'
pparams.sigma_peak = 0.8;
pparams.sigma_noise = 8.7076e+03;
pparams.numpy= 1;
pparams.alpha = 0.75;
pmat = [];

mask_holder = genMask(chrom, 'ion', 100);
    %% all error checking block before we attempt to do anything.
    if isa(chrom, 'chrom2gram') % chrom2gram object passed.
        if ~isempty(chrom.p)
            disp('peak detection already performed for this sample, returning original p vector');
            pmat = chrom.p;
        else
            chrom = chrom.getRavel;
            szc = size(chrom);
        end
    else % tensor structure passed.
        szc = size(chrom);
    end
    % either way if we get here we have a tensor and we know its size
    if isempty(pmat)
       pmat = zeros(szc);
       for i = 1:size(mask_holder,1) % iterate over the masks.
       % check dimensions of the output and the perform appropriate peak
       % detection or transform first
            mini_tic = squeeze(sum(sum(reshape(chrom(mask_holder(i,:)), [szc(1)/size(mask_holder,1),szc(2), szc(3)]),1),3));
            p = getPeaksConv(1:numel(mini_tic), mini_tic, pparams.sigma_peak, pparams.sigma_noise, pparams.alpha, pparams.numpy, 0);
            pmat(mask_holder(i,:)) =  repmat(p, [szc(3),1,szc(1)/size(mask_holder,1)]);
       end
       disp('peak detection performed.');
    end
end % end function maskedPeakDetection

function m = genMask(c, mtype, param)
    tic;
    sz = size(c.getRavel);
    m.original_size = sz;
    m.type = mtype;
    switch mtype
%%%%%%%%%%%%%%%%%%%%%%%%
        case 'rows'
        %% 1-D rows pattern - for GCxGC-MS
        disp('row-wise mini tics for GCxGC-MS data selected');
        MASK_LEVELS = param;
        mask_idx = round(linspace(0,sz(1),MASK_LEVELS+1));
        mask_holder = int32(zeros([MASK_LEVELS, mean(diff(mask_idx))*sz(2)*sz(3)])); % preallocate

        for i  = 2:numel(mask_idx)
            [X,Y,Z] = meshgrid([mask_idx(i-1)+1:mask_idx(i)],[1:sz(2)], [1:sz(3)]);
            mask_holder(i-1,:) = int32(sub2ind(sz,X(:),Y(:),Z(:)));
        end
%%%%%%%%%%%%%%%%%%%%%%%%
        case 'venetian'
            %% venetian blinds pattern - for GCxGC-MS
            blind_spacing  = param;
            disp('staggered strips - venetian blind mask specified');
%%%%%%%%%%%%%%%%%%%%%%%%        
        case 'tile'
            %% Tile based, using 2-D tile structures.
            tile_idx = c.gentrerateTiles;
%%%%%%%%%%%%%%%%%%%%%%%%
        case 'ion'
            %% unravelled chrom2gram once per nominal mass channel
            mask_holder = int32(reshape([1:prod(sz)], [sz(2)*sz(1), sz(3)])');
    end% end of the switch block
    b = toc;
    disp([mtype, ' mask generated... ',num2str(b), ' seconds elapsed.']);
end