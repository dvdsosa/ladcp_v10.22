function [data,params,values,messages]=calc_soundsp(data,params,values,messages)
% function [data,params,values,messages] =...
%	calc_soundsp(data,params,values,messages)
%
% Calculate the sound speed profile from the best available source
% and correct the ADCP velocity data for it
%
% version 0.11	last change 27.11.2015

% Matin Visbeck 
%  December 2002, LDEO
% G.Krahmann, IFM-GEOMAR, Jul 2005

% wrong sign in handling CTD-Prof                 GK Jun 2007     0.1->0.2
% sc had wrong points set compared with data.rX   GK, 07.03.2012  0.2-->0.3
% using now preferably the CTDprof data to calculate the sound speed
% for every time step. For unkown reasons salinity was fixed to
% 34.5 for the sound speed calculation in the ctdtime data
%                                                 GK, 05.07.2013  0.3-->0.4
% replaced interp1q with interp1                  GK, 11.07.2013  0.4-->0.5
% catch case of all NaN sound speeds              GK, 16.04.2014  0.5-->0.6
% fix small bug when only ctdtime was given       GK, 13.05.2014  0.6-->0.7
% set zbottom to 0 to force recalculation         GK, 07.07.2015  0.7-->0.8
% calculate depth from CTD pressure               GK, 15.07.2015  0.8-->0.9
% catch unknown lat                               GK, 15.09.2015  0.9-->0.10
% fix bug with preset bottom depth                GK, 27.11.2015  0.10-->0.11


%
% general function info
%
disp(' ')
disp('CALC_SOUNDSP:  calculate sound speed profile and correct the raw data')


values.GEN_Sound_sp_calc = '[NA]';
if ~isempty(data.ctdprof)

  disp('    Calculating soundspeed from CTD profile ')
  if isnan(values.lat)
    zctd = sw_dpth(data.ctdprof(:,1),0);
  else
    zctd = sw_dpth(data.ctdprof(:,1),values.lat);
  end
  zctd(1) = -1e4;
  zctd(end) = 1e4;
  data.ctdprof_ss = sw_svel(data.ctdprof(:,3),data.ctdprof(:,2),...
		data.ctdprof(:,1));
  data.ctdprof_z = sw_dpth(data.ctdprof(:,1),values.lat);
  data.ss = interp1(-zctd,data.ctdprof_ss,data.z')';
  values.GEN_Sound_sp_calc = '[CTD]';
  
elseif ~isempty(data.ctdtime_data)

  pp = sw_pres(abs(data.z),values.lat);
  disp('    Calculating soundspeed from CTD pressure and temp timeseries')
  %data.ss = sw_svel(34.5*ones(1,length(pp)),data.ctdtime_data(:,2)',pp);
  data.ss = sw_svel(data.ctdtime_data(:,3)',data.ctdtime_data(:,2)',pp);
  values.GEN_Sound_sp_calc = '[CTD]';

else

  pp = sw_pres(abs(data.z),values.lat);
  disp('    Calculating soundspeed from integrated w and ADCP temp')
  data.ss = sw_svel(34.5*ones(1,length(pp)),data.temp(1,:),pp);
  values.GEN_Sound_sp_calc = '[T-P]';

end


%
% apply sound speed correction to ADCP velocities
%
% later there will be a correction of the bin distance from the ADCP
%
disp('    Correcting all ADCP velocities for sound speed ')
sc = meshgrid(data.ss./data.sv(1,:),data.izd);
if values.up==1
  sc = [sc;meshgrid(data.ss./data.sv(2,:),data.izu)];
  sc = flipud(sc);
end
if any(~isnan(sc(:)))
  data.ru = data.ru.*sc;
  data.rv = data.rv.*sc;
  data.rw = data.rw.*sc;
  if isfield(data,'hbot')
    data.hbot = data.hbot.*sc(end,:);
    data.bvel(:,1:3) = data.bvel(:,1:3).*sc(end-[0:2],:)';
    if isfield(data,'bvel_rdi')
      data.bvel_rdi(:,1:3) = data.bvel_rdi(:,1:3).*sc(end-[0:2],:)';
    end
    if isfield(data,'bvel_own')
      data.bvel_own(:,1:3) = data.bvel_own(:,1:3).*sc(end-[0:2],:)';
    end
  end
  if isfield(data,'hsurf')
    data.hsurf = data.hsurf.*sc(1,:);
  end
  if params.zbottom>0 & params.zbottom_orig==0
    params.zbottom = 0;
  end
end

