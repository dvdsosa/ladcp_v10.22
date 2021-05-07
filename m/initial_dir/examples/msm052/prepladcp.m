function [file] = prepladcp(stn)
% function [file] = prepladcp(stn)
%
% prepare LADCP data for LADCP processing
%
% we need the raw LADCP data to be in the correct place and
% have the correct names.
%
% THIS FILE IS CRUISE SPECIFIC
%
% to create a file for your own cruise, modify this file
%
% you will just need to copy and possibly rename the files
% In case of old BB and NB systems you might need to append
% the raw data files.
%
% the convention for filenames is
%
% xxxDN000.000  and  xxxUP000.000  	with xxx the 3-digit station number
%
% they need to be copied into one directory per station
% data/raw_ladcp/xxx		with xxx the 3-digit station number

% G.Krahmann, IFM-GEOMAR, Aug 2005

%disp('YOU FIRST NEED TO EDIT THE FILE cruise_id/m/prepladcp.m !')
%return


if ~exist(['data/raw_ladcp/',int2str0(stn,3)])
  eval(['!mkdir data/raw_ladcp/',int2str0(stn,3)])
end
eval(['!cp /home/ladcp/msm052/data/raw/',int2str0(stn,3),'/',...
	int2str0(stn,3),'DN*.* data/raw_ladcp/',int2str0(stn,3)])
eval(['!cp /home/ladcp/msm052/data/raw/',int2str0(stn,3),'/',...
	int2str0(stn,3),'UP*.* data/raw_ladcp/',int2str0(stn,3)])

% set file name
file = ['data/raw_ladcp/',int2str0(stn,3)];
