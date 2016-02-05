%S I (Steffan) will comment below. To distinguish our comments, mine will 
%S be denoted '%S'

function c = parse_UofA_CDF(path2File)
% system specification
% 200 Hz
%S Analytes are passed from the second column to the MS at 200 Hz
% 2 second d2 retention time (400dp) !!! HARDCODED for our system !!! 
c.d2 = 400; %S 2 s retention time * 200 Hz = 400
% 10-500 mz range
% 52 minute run time

ncid = netcdf.open(path2File, 'NOWRITE'); %netcdf file pointer into open file
%S netcdf can be used to call functions in the NetCDF library. In this 
%S case it is used to open the cdf file in a read-only mode and to denote 
%S it ncid

if strcmpi(netcdf.inqFormat(ncid), 'FORMAT_CLASSIC') 
    disp('netCDF classic format detected, all is well...');
    %S inqFormat returns the format of the CDF file. The return is 
    %S compared (case insensitive) to 'FORMAT_CLASSIC'; if they are equal,  
    %S a logical 1 is returned by strcmpi, and a message is displayed
else
    disp('incompatible net CDF format detected, unexpected performance may occur');
    %S If the format is not classic, another message is displayed
end

vars_in_cdf = netcdf.inqVarIDs(ncid);
%S Returns the identifiers (numbers) of the variables in the imported CDF 
%S file
for i = 1:numel(vars_in_cdf) 
    %S Loop runs over all variables in the CDF
    varname = netcdf.inqVar(ncid,vars_in_cdf(i));
    %S Returns the name for a variable
    data = netcdf.getVar(ncid,vars_in_cdf(i));
    %S Returns the values of a variable
    if and(isa(data(1), 'integer'), ~isequal(data(1), data(end)))
        eval(['c.',varname, ' = data;']); 
    end
    %S If the first value in a variable is of the type 'integer', and the 
    %S first and last value in a variable are not equal, create a string 
    %S 'c.VARNAME = data;', e.g. 'c.mass_values = data;' and evaluate it as
    %S code
    
    if isa(data(1), 'double') 
    %S If the first value in a variable is of the type 'double' 
        %if std(data)<1.0e-10
        if range(data)==0    % You could also use this if you want to
        %S know if it is always the same value. On my machine it is faster 
            eval(['c.',varname, ' = data(1);']); %only need one if it is always the same value 
        %S And if the standard deviation of the values in a variable is
        %S very small, i.e. the data is always the same value, create a
        %S string 'c.VARNAME = data(1);', e.g. 'c.mass_values = data(1);' 
        %S and evaluate it as code 
        else
            eval(['c.',varname, ' = data;']); 
        %S If the first value in a variable is of the type 'double' and if 
        %S the standard deviation of the values in a variable is not very 
        %S small, i.e. the data is not constant, create a string 
        %S 'c.VARNAME = data;', e.g. 'c.mass_values = data;', and evaluate
        %S it as code
        
        end
    end
end

netcdf.close(ncid); %close pointer into the file
%S Why would you close the file? If you implement no further code hereafter
%S that changes the netcdf file, is it still better to close the cdf file 
%S for memory reasons?

clearvars -except c
%S remove all variables from the workspace, except for c. You can
%S do this because all variables that are needed have been introduced as 
%S properties of the object 'c'

c.masses = int16(c.mass_range_min:c.mass_range_max);
%S Set the masses property of 'c' to be an array of 16-bit signed integers,
%S from the smallest m/z to the largest m/z obtained from the cdf file. You
%S use 16-bit integers because they are Matlab's data type using least 
%S memory that can accomodate typical m/z values (up until m/z = 32768). 

c.xyz = int32(zeros([numel(c.scan_index), numel(c.masses)]));
%S Set the xyz property of 'c' to be a 32-bit signed integer array sized
%S {number of MS scans} by {number of masses} filled with zeros. By filling
%S the array with zeros, an appropriate amount of memory is reserved, so 
%S any future changes to xyz will not be slowed down because of memory
%S checking 

c.d1 = numel(c.scan_index)/c.d2;
%S Set the d1 property of 'c' (separation dimension 1) to the total number
%S of MS scans divided by the number of  scans made for each 2nd dimension
%S modulation in order to obtain the number of modulations, which is the
%S length of the first dimension

if ~round(c.d1)==c.d1
    error('data does not neatly divide into d2 and d1 dimensions!');
end
%S Check if d1 is an integer or not. It should be since it is a
%S representation of the number of GCxGC modulations. If it is not, an
%S error will be thrown

disp('all variables loaded into memory, reformatting structure...');
%S Progress message

%S scan_index: increases with 488 for each index, 624000
%S total_intensity: TIC, 624000
%S intensity_values: mass intensity, 304570529

for i = 1:numel(c.total_intensity)-1 %S runs over the number of mass scans
    vals = c.intensity_values(c.scan_index(i)+1:c.scan_index(i+1)); 
    %S set vals to be the mass intensities for all the masses that were
    %S scanned in a specific scan (the masses that were scanned are
    %S variable from scan to scan)
    mass_idx = c.mass_values(c.scan_index(i)+1:c.scan_index(i+1));
    %S set mass_idx to be the m/z values that were scanned in a specific 
    %S scan (those are variable from scan to scan)
    
    c.xyz(i, mass_idx) = vals; %
    %S c.xyz is an array of all scans vs all m/z values that could be
    %S detected on the spectrometer. For those m/z that were scanned, the
    %S values are plugged in c.xyz
    
    if ~isequal(sum(c.xyz(i,c.mass_range_min:end)), c.total_intensity(i)) %make sure we are doing things right
         disp('sum of masses does not equal corresponding TIC value');
    %S Check if the TIC for a scan provided by the machine equals the sum 
    %S of the intensities measured for that scan
    end
end

c.mass_range_min = 10; %S Why are they set here again? They have already been used before 
c.mass_range_max = 500; %S And couldnt you just use min(c.mass_values) instead of the code below

% % verify mass range the hard way!
% [bin_counts, bin_centers] = hist(double(c.mass_values), [1:550]);
% mass_idx_histogram = find(bin_counts~=0, 1, 'first');
% c.mass_range_min = bin_centers(mass_idx_histogram);
% % lower end
% 
% mass_idx_histogram = find(bin_counts~=0, 1, 'last');
% c.mass_range_max = bin_centers(mass_idx_histogram);
% % upper end

c.xyz(:,1:c.mass_range_min-1)=[];
%S Delete those columns of c.xyz that fall below the minimum m/z value that
%S is scanned. I would suggest to also delete any columns that relate to 
%S m/z values > mass_range_max, even if it is not the case with our current
%S dataset.
%S c.xyz(:,c.mass_range_max+1:end)=[];

c.masses = [c.mass_range_min:c.mass_range_max];
%S Set c.masses to the range of masses that is scanned

clearvars -except c
c = rmfield(c, {'mass_values', 'intensity_values'});
c = cleanFields(c);
%S Delete those variables of c that are not needed anymore: mass_values 
%S and intensity_values, since these have been moved to c.xyz. Also, call
%S the cleanFields function to be operated on c

c = chrom2gram(c); % call matlab class constructor!
%[c.tile_idx, ~] = gc2d_recommendTiles(c.d1,c.d2,c.d1,c.d2,45,45);
%c.p = gc2d_peakDetect(c);
%c.baseline = gc2d_baseline(c);

end

%S Am I correct that you can introduce multiple functions in one function
%S file?

function c = cleanFields(c)
    c_idx = fieldnames(c); %S set c_idx to be the names of the fields in structure c
    for i = 1:numel(c_idx) %S for each field
        if and(size(c.(c_idx{i}),1)==1, size(c.(c_idx{i}),2)==1) %scalar fields only
            %S You could also do
            %S if size(c.(c_idx{i})) == [1 1]
            if abs(c.(c_idx{i})) == 9999 %only bullshit filler value anyways.
                c = rmfield(c, c_idx(i)); %S remove fields with scalar filler values
            end
        elseif or(and(size(c.(c_idx{i}),1)>1, size(c.(c_idx{i}),2)==1),and(size(c.(c_idx{i}),1)==1, size(c.(c_idx{i}),2)>1))% only looking at vector fields
            %S You could also do, shorter: 
            %S elseif min(size(c.(c_idx{i}))) == 1
            if numel(unique(c.(c_idx{i}))) == 1
                c.(c_idx{i}) = unique(c.(c_idx{i}));
            %S You could also do
            %S if range(c.(c_idx{i})) == 0
            %S    y = c.(c_idx{i});
            %S    c.(c_idx{i}) = y(1);
            end
        else
            disp('skipping multi-dimensional or tensor field')
        end
    end
end

%S What is the difference between using c_idx{i} and c_idx(i)?

