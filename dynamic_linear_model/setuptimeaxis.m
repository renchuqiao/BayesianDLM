% setup time axis for labelling and neat plotting

ti=char(time(2:end,1)); %ti selects out all the dates for 1762 days e.g. "20070717". size(ti)= 1762x8
Vnames=Ynames(1,3:3:end);
BAnames=Ynames(1,4:3:end); 
Ynames=Ynames(1,2:3:end); 
names=char(Ynames); 
 
[iyr,tticks]=unique(str2num(ti(:,3:4)));  % iyr is unique years (in 2000s) %ti is all the dates listed as 8 digits 
%we only interested in the 3rd and 4th digits of ti coz those denote the
%year, say 08, 09, 10, 11,..., 

[imo12,tticks12]=unique(str2num(ti(1:(tticks(2)-1),5:6)));
[imo13,tticks13]=unique(str2num(ti((tticks(2)):(tticks(3) - 1),5:6)));
tticks13 = tticks13 + tticks(2) - 1;
[imo14,tticks14]=unique(str2num(ti((tticks(3)):end,5:6)));
tticks14 = tticks14 + tticks(3) - 1;
tticks_all = [tticks12', tticks13', tticks14']';
tticks_3 = tticks_all(2:2:end);
                                          % tticks is index of first occurrence of each 
                                          
                 % drop first tick mark               
%   i = 1:2:length(tticks); iyr=iyr(i); tticks=tticks(i);  % to use fewer,  equally spaced tick marks 
yr=ti(tticks_3,3:4); 
%tticks is the NUMBER of dates of the first day of each year
%ti(tticks,:) gives the 8-digit dates of all the first days of each year
%yr=ti(tticks,3:4);  gives the unique years: 08, 09, 10, 11, 12, 13, 14. 

mo=ti(tticks_3,5:6); %gives the months of the first days of each year. here mo is a 7x2 and each entry is 01
da=ti(tticks_3,7:8);  % char year, month, day at ticks. %gives the dates of the first days of each year. here da is a 7x2 and [02, 02, 04, 03, 03, 02, 02]
tdates=[]; for i=1:length(tticks_3), tdates=[tdates mo(i,:)  '/' yr(i,:) '|']; end; %'/' da(i,:) '|' ]; end
%size tdates = 1x42    %%%%tdates =01/08|01/09|01/10|01/11|01/12|01/13|01/14|
xa=['set(gca,''Xtick'',tticks_3);set(gca,''XtickLabel'',tdates);xlabel(''month/year'',''FontSize'',22); box off;'];

