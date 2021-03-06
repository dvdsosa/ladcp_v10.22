function prepnav(stn,values)
% function prepnav(stn,values)
%
% prepare navigational data for LADCP
% we an array 'data' containing the 2 columns
% latitude in decimal degrees    longitude in decimal degrees
% and the vector 'timnav' containing the time of the navigational
% data in Julian days
%
% THIS FILE IS CRUISE SPECIFIC
%
% to create a file for your own cruise, modify this file
%
% The navigational data should be at a resolution of 1 per second.
% Lower resolution will lead to worse processing results.

% G.Krahmann, IFM-GEOMAR, Aug 2005

% if you do no have navigational data to be used in the
% LADCP processing, uncomment the next line

%disp('YOU FIRST NEED TO EDIT THE FILE cruise_id/m/prepnav.m !')
%return


% first copy navigational data to the raw NAV data directory
% data/raw_nav
% In our example here we used a Seabird CTD system which also
% recorded navigational data. Thus we use the same data as
% for the CTD-TIME data
%eval(['!cp /nbp/science/NBP0408/ladcp/ctd/p408',int2str0(stn,3),...
%  	'ladcp.cnv data/raw_ctd/p408',int2str0(stn,3),'.tim'])

% load this data and convert to standard format
%
% in this example
% we skip the header of the file and extract the lat and lon columns
% into 'data' and the time vector into 'timctd'
%
% you will have to make sure that the time is stored in Julian days
% 
% in this example the time came from the Seabird CTD software
% which records only day of year within the record, so
% we had to add the year
[hdr,data] = hdrload(['data/raw_ctdtime/msm_04_1_',int2str0(stn,3),'s1s2s3s4s5_1sec.cnv']);
timnav = data(:,2) + julian([2006,1,0,0,0,0]);
data = data(:,[3,4]);

% To reduce the amount of data we crop the navigational data to
% the same time as the CTD-TIME data. In our example case that
% was an unnecessary exercise since they are the same data, but if
% you have independent navigational data (e.g. daily navigational files)
% this will reduce file size.
good = find(timnav>=values.start_cut & timnav<=values.end_cut);
timnav = timnav(good);
data = data(good,:);

% store data in the standard location
if str2num(version('-release'))>=14
  eval(['save data/nav/nav',int2str0(stn,3),' timnav data -v6'])
else
  eval(['save data/nav/nav',int2str0(stn,3),' timnav data'])
end
