a = dir('/Users/Martin/SURFdrive/GCxGC-MS_fireDebris/*.cdf');
for i = 1:numel(a)
    eval(['c = parse_UofA_CDF(','''/Users/Martin/SURFdrive/GCxGC-MS_fireDebris/',a(i).name,''');']);
    c.p = maskedPeakDetect(c);
    clearvars -except a i c
    eval(['save(','''/Users/Martin/SURFdrive/GCxGC-MS_fireDebris/',a(i).name(1:end-3),'mat'',','''-v7.3'',','''c'');']);
end