addpath ../../lisca/src/
a = dir('/Users/Martin/SURFdrive/GCxGC-MS_fireDebris/*.cdf');
b = dir('/Users/Martin/SURFdrive/GCxGC-MS_fireDebris/*.mat');
t = regexp([b(:).name], '.mat','split');
t(cellfun(@isempty, t)) = [];

for i = 1:numel(a)
    if ~ismember(a(i).name(1:end-4),t)
        eval(['c = parse_UofA_CDF(','''/Users/Martin/SURFdrive/GCxGC-MS_fireDebris/',a(i).name,''');']);
        c.p = maskedPeakDetect(c);
        % solvent peak mask
        % (and(c.p(:,:,c.masses==84)>300, c.p(:,:,c.masses==49)>300));
        clearvars -except a i c t
        eval(['save(','''/Users/Martin/SURFdrive/GCxGC-MS_fireDebris/',a(i).name(1:end-3),'mat'',','''-v7.3'',','''c'');']);
    else
        disp('sample already processed... moving on.');
    end
end