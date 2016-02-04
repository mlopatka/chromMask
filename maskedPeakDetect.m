function pmat = maskedPeakDetect(chrom)
%
addpath '../../lisca/src/'; % path to the probabilistic peak detection code.
% getPeakConv is found at https://github.com/mlopatka/getPeakConv
%
% hardcoding some essential things that should stay relatively constant over a batch processing paradigm.
%% batch processing parameter choices
pparams.sigma_peak = 0.95;
pparams.sigma_noise = 154.0573;
pparams.numpy= 2;
pparams.alpha = 0.75;
pmat = [];

%% begin function
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
       pmat = int16(zeros(szc)); % trade memory for processing and lose all precision beyond 3rd decimal place
       
       %% different behaviour for different chromatography and mask combinations
       if numel(szc) > 3 % nope nope nope nope
           error('4-D tensor passed to peak detection... noping the F out!');
        
       elseif numel(szc) == 3 % can only be GCxGC-MS
           if strcmpi(mask_holder.type, 'row') % row masking, need to sum accros multiple dimensions
               for i = 1:size(mask_holder.mask,1) % iterate over the row masks.
                   mini_tic = (sum(sum(reshape(chrom(mask_holder.mask(i,:)), [szc(1)/size(mask_holder.mask,1),szc(2),szc(3)]),1),3));
                   p = getPeaksConv(1:numel(mini_tic), mini_tic, pparams.sigma_peak, pparams.sigma_noise, pparams.alpha, pparams.numpy, 0);
                   pmat(mask_holder.mask(i,:)) = repmat(repmat(int16(1000*p),[(szc(1)/size(mask_holder.mask,1)),1]), [1,1,szc(3)]);
               end % end iterate over masks
           elseif strcmpi(mask_holder.type, 'ion')
               for i = 1:size(mask_holder.mask,1) % iterate over the row masks.
                   tic;
                   if i == 19
                       keyboard
                   end
                   mini_tic = chrom(mask_holder.mask(i,:));
                   % need to downsample for feasibility
                   mini_tic = downsample(mini_tic,10);
                   p = getPeaksConv(1:numel(mini_tic), mini_tic, pparams.sigma_peak, pparams.sigma_noise, pparams.alpha, pparams.numpy, 0);
                   p = interp1([1:numel(p)],p,linspace(1,numel(p),10*numel(p)));% upsample with rough interp.
                   pmat(mask_holder.mask(i,:)) = reshape(int16(1000*p), [szc(1), szc(2)]);
                   toc
               end % end iterate over masks
           elseif strcmpi(mask_holder.type, 'tile')
               
           elseif strcmpi(mask_holder.type, 'venetian')
               
           else
               error('unrecognzed mask type designation, please update genMask function');
           end
       elseif numel(szc) == 2% can be GC-MS or GCxGC-FID
           
           
       else % can only be GC-FID
           % cant really do much masking, process as a vector and return p
           pmat = int16(1000*getPeaksConv(1:numel(chrom), chrom, pparams.sigma_peak, pparams.sigma_noise, pparams.alpha, pparams.numpy, 0));
       end % end condition deciding how to handle data based on dimensions/        
    end % end generation of pmat with masked data
    disp('peak detection performed.');
end % end function maskedPeakDetection

function mask_holder = genMask(c, mtype, param)
    tic;
    sz = size(c.getRavel);
    mask_holder.type = mtype;
    switch mtype
%%%%%%%%%%%%%%%%%%%%%%%%
        case 'row'
        %% 1-D rows pattern - for GCxGC-MS
        disp('row-wise mini tics for GCxGC-MS data selected');
        MASK_LEVELS = param;
        mask_idx = round(linspace(0,sz(1),MASK_LEVELS+1));
        mask_holder.mask = int32(zeros([MASK_LEVELS, mean(diff(mask_idx))*sz(2)*sz(3)])); % preallocate
        for i = 2:numel(mask_idx)
            [X,Y,Z] = ndgrid([mask_idx(i-1)+1:mask_idx(i)],[1:sz(2)], [1:sz(3)]);
            mask_holder.mask(i-1,:) = int32(sub2ind(sz,X(:),Y(:),Z(:)));
        end
%%%%%%%%%%%%%%%%%%%%%%%%
        case 'venetian'
            %% venetian blinds pattern - for GCxGC-MS
            blind_spacing  = param; % now the parameter is the space between blinds
            disp('staggered strips - venetian blind mask specified');
%%%%%%%%%%%%%%%%%%%%%%%%        
        case 'tile'
            %% Tile based, using 2-D tile structures.
            tile_idx = c.gentrerateTiles;
            disp('2-D tiles mask specified');
%%%%%%%%%%%%%%%%%%%%%%%%
        case 'ion'
            %% unravelled chrom2gram once per nominal mass channel
            mask_holder.mask = int32(reshape([1:prod(sz)], [sz(2)*sz(1), sz(3)])');
            disp('ion by ion mask defined in nominal mass channel dimension');
    end% end of the switch block
    b = toc;
    disp([mtype, ' mask generated... ',num2str(b), ' seconds elapsed.']);
end