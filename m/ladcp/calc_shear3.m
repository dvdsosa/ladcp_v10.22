function [ds,dr,messages] = calc_shear3(d,p,ps,dr,messages)
% function [ds,dr,messages] = calc_shear3(d,p,ps,dr,messages)
%
% - compute shear profiles 
% - use only central difference
% - use 2*std editing
% 
% version 3  last change 16.5.2015

%  Martin Visbeck, LDEO, 3/7/97
% some modifications and code cleanup                          GK, 16.5.2015  2-->3

%
% general function info
%
disp(' ')
disp(['CALC_SHEAR3:  calculate a baroclinic velocity profile based on shears only'])


%
% resolution of final shear profile in meter
%
% I am not sure whether this makes sense, as the common case is
% dz=bin_length  but we use here central differences for the shear which are
% then over 2*dz distances.  GK, May 2015
%
dz = ps.dz;
disp(['    Averaging shear profile over (ps.dz) ',num2str(dz),' [m]'])


%
% Discard a certain amount of data as suspected outliers that
% for relatively small numbers of values will skew the mean.
% In the old calc_shear2.m the allowed range was the std of the shears 
% times the set factor stdf. Default was stdf=2 .
% In calc_shear3.m the outlier determination is iterative and thus a bit
% safer in calculation. The allowed range is set as a fraction of the whole
% population. For gaussian distributions stdf converts to this fraction.
% As usually stdf is not varied much, we have implemented a lookup list
% here.
%
use_new_outlier = 1;
stdf = ps.shear_stdf;
if use_new_outlier==1
  fracs = 1 - [31.8,13.4,4.5,1.3,0.3]/100;
  stdfs = [1,1.5,2,2.5,3];
  [dummy,ind] = nmin(abs(stdf-stdfs));
  frac = fracs(ind);
end
disp(['    Maximum allowed std within calculation intervals : ',num2str(stdf)])
disp(['    Data deviating more from the median will be discarded.'])


% don't know what that is for, so comment it   GK May 2015
%if nargin<4, 
%  dr.dummy = 0; 
%end


%
% check if only one istrument is to be used
%
if ps.up_dn_looker==2
  % down looker only
  d.weight(d.izu,:)=nan;
elseif ps.up_dn_looker==3
  % up looker only
  d.weight(d.izd,:)=nan;
end


%
% weightmin
%
% A problem in the old calc_shear2.m was the usage of the wrong weights to
% reject velocity pairings for the shear. For the central shears it was
% using the weights from the central velocity which is not used in the
% central difference calculation at all. In calc_shear3.m this was fixed by
% using the combined weight of the two velocities that actually enter the
% calculation of the central differences shear. Unfortunately this lead to
% significantly fewer shear values making the whole calculation more
% uncertain. 
% Thus the whole idea of using the weights from the regular inversion was
% not used anymore. Instead ALL water values (as stored in d.izmflag) except 
% for the first bins in up and downlooker are now used. The typical weight 
% minimum from calc_shear2.m had no effect anyway as the default threshold
% was set very low.
%
if 0
  disp(['    Minimum weight  ',num2str(ps.shear_weightmin),' for data to be used in shear calc.'])

  is1 = sum(isfinite(d.weight(:)));
  local_weight = double(d.weight>ps.shear_weightmin);
  local_weight = replace( local_weight, local_weight==0, nan );
  is2 =  sum(isfinite( local_weight(:) ));
  disp(['    Will use ',int2str(is2/is1*100),' % of data.'])
else
  local_weight = d.izmflag+1;
  if isfield(d,'up')
    localweight([d.izu(1),d.izd(1)],:) = nan;
  else
    localweight([d.izd(1)],:) = nan;
  end
end

%
% compute shear
%
% two ways are offered here
% first:   central differences for the shears
% second:  single differences 
% the first is similar to the ways of the old calc_shear2.m
%
% central differences
if 1
  local_weight = [repmat(nan,1,size(local_weight,2));diff2(local_weight);repmat(nan,1,size(local_weight,2))]+1;
  ushear = [NaN*d.ru(1,:);diff2(d.ru(:,:))./diff2(d.izm);NaN*d.ru(1,:)].*local_weight;
  vshear = [NaN*d.rv(1,:);diff2(d.rv(:,:))./diff2(d.izm);NaN*d.rv(1,:)].*local_weight;
  wshear = [NaN*d.rw(1,:);diff2(d.rw(:,:))./diff2(d.izm);NaN*d.rw(1,:)].*local_weight;
  zshear = -d.izm;
  diff_step = 2;
% single differences
else
  ushear = diff( d.ru.*local_weight )./diff(d.izm);
  vshear = diff( d.rv.*local_weight )./diff(d.izm);
  wshear = diff( d.rw.*local_weight )./diff(d.izm);
  zshear = -(d.izm(1:end-1,:)+d.izm(2:end,:))/2;
  diff_step = 1;
end


%
% set depth levels
%
z = dr.z;


%
% prepare shear solution result vectors
%
ds.usm = repmat(nan,length(z),1);
ds.vsm = ds.usm;
ds.wsm = ds.usm;
ds.usmd = ds.usm;
ds.vsmd = ds.usm;
ds.use = ds.usm;
ds.vse = ds.usm;
ds.wse = ds.usm;
ds.nn = ds.usm;
ds.z = z;


%
% loop over depth levels and calculate the average shear at that level
%
for n=[1:length(z)]

  i2 = find( abs( zshear - z(n) ) <= dz/2*diff_step);
  i1 = i2( isfinite( ushear(i2) + vshear(i2) ) );
  ds.nn(n) = length(i1);
  if ds.nn(n) > 2

    % two ways to select outliers
    % first:   select all that are beyond a fixed range around the median
    % second:  iteratively reject the worst (largest distance from mean)
    %          until a fixed fraction is rejected
    % the second is usually the safer calculation but is a bit slower
    if 0
      usmm = median( ushear(i1) );
      ussd1 = std( ushear(i1) );
      vsmm = median( vshear(i1) );
      vssd1 = std( vshear(i1) );
      wsmm = median( wshear(i1) );
      wssd1 = std( wshear(i1) );
      ii1 = i1( find(abs(ushear(i1)-usmm)<stdf*ussd1) );
      ii2 = i1( find(abs(vshear(i1)-vsmm)<stdf*vssd1) );
      ii3 = i1( find(abs(wshear(i1)-wsmm)<stdf*wssd1) );
    else
      [dummy,ii1] = meanoutlier(ushear(i1),frac);
      [dummy,ii2] = meanoutlier(vshear(i1),frac);
      [dummy,ii3] = meanoutlier(wshear(i1),frac);
      ii1 = i1(ii1);
      ii2 = i1(ii2);
      ii3 = i1(ii3);
    end

    % two ways of calculating the mean and std of the selected shears
    % first:  if there is a rejected one in any of u,v,w shears then use it
    %         for non of the calculations
    % second: if there is a rejected one in any of u,v,w shears then use it
    %         only in u,v or w calculations
    % the second one is the one used by the old calc_shear2.m
    % but to me this does not make sense, GK May 2015
    if 1
      dummy = zeros(size(ushear));
      dummy(ii1) = 1;
      dummy(ii2) = dummy(ii2)+1;
      dummy(ii3) = dummy(ii3)+1;
      ii = find(dummy==3);
      if length(ii)>1
        ds.usm(n) = mean(ushear(ii));
        ds.usmd(n) = median(ushear(ii));
        ds.use(n) = std(ushear(ii));
        ds.vsm(n) = mean(vshear(ii));
        ds.vsmd(n) = median(vshear(ii));
        ds.vse(n) = std(vshear(ii));
        ds.wsm(n) = mean(wshear(ii));
        ds.wse(n) = std(wshear(ii));
      end
    else
      if length(ii1)>1
        ds.usm(n) = mean(ushear(ii1));
        ds.usmd(n) = median(ushear(ii1));
        ds.use(n) = std(ushear(ii1));
      end
      if length(ii2)>1
        ds.vsm(n) = mean(vshear(ii2));
        ds.vsmd(n) = median(vshear(ii2));
        ds.vse(n) = std(vshear(ii2));
      end
      if length(ii3)>1
        ds.wsm(n) = mean(wshear(ii3));
        ds.wse(n) = std(wshear(ii3));
      end
    end
  end

end


%
% a debugging figure
%
if 0
sfigure(3);
clf
orient tall
subplot(1,2,1)
plot(ushear,zshear,'b.','markersize',3)
hold on
plot(ds.usm,ds.z,'r')
plot(ds.usmd,ds.z,'k')
inv_shear_u = -diff(dr.u-mean(dr.u))/dz;
plot(inv_shear_u,(z(1:end-1)+z(2:end))/2,'g')
set(gca,'ydir','reverse')

subplot(1,2,2)
plot(vshear,zshear,'b.','markersize',3)
hold on
plot(ds.vsm,ds.z,'r')
plot(ds.vsmd,ds.z,'k')
inv_shear_v = -diff(dr.v-mean(dr.v))/dz;
plot(inv_shear_v,(z(1:end-1)+z(2:end))/2,'g')
set(gca,'ydir','reverse')

sfigure(2)
end



%
% integrate shear profile (from bottom up)
%
ds.usm = replace( ds.usm, isnan(ds.usm), 0 );
ds.vsm = replace( ds.vsm, isnan(ds.vsm), 0 );
ds.wsm = replace( ds.wsm, isnan(ds.wsm), 0 );

ds.ur = flipud(cumsum(flipud(ds.usm)))*dz;
ds.vr = flipud(cumsum(flipud(ds.vsm)))*dz;
ds.wr = flipud(cumsum(flipud(ds.wsm)))*dz;
ds.ur = ds.ur-mean(ds.ur);
ds.vr = ds.vr-mean(ds.vr);
ds.wr = ds.wr-mean(ds.wr);

%
% This is a peculiar place for the single ping error estimate. But 
% as it is based on the variability in the data itself, it makes sense.
% The assumption is that there should be basically zero shear in the
% vertical velocities. At least it is so small as to be not detectable
% here. Thus any variability in the vertical shear is caused by the
% errors/noise of the measurement. Together with an angular conversion
% factor this gives an error/noise value for the horizontal velocities.
%
dz2 = diff_step*abs(mean(diff(d.zd)));
if isfield(d,'down')
  fac = 1/tan(d.down.Beam_angle*pi/180)*sqrt(2)*dz2;
else
  fac = 1/tan(d.up.Beam_angle*pi/180)*sqrt(2)*dz2;
end
ds.ensemble_vel_err = ds.wse*fac;
dr.ensemble_vel_err = ds.wse*fac;


%
% store result and give text output
%
dr.u_shear_method = ds.ur;
dr.v_shear_method = ds.vr;
dr.w_shear_method = ds.wr;
uds = nstd(dr.u-mean(dr.u)-ds.ur);
vds = nstd(dr.v-mean(dr.v)-ds.vr);
uvds = sqrt(uds^2+vds^2);
if uvds>nmean(dr.uerr)*1.5
  error_increase_factor = 1/nmean(dr.uerr)*uvds/1.5;
  warn = ('>   Increasing error estimate because of elevated shear - inverse difference');
  disp(warn)
  disp(['>   by a factor of ',num2str(error_increase_factor)])
  disp(['>   std of difference between regular and shear profile : ',num2str(uvds),' [m/s]'])
  messages.warn = strvcat(messages.warn,warn);
  dr.uerr = dr.uerr * error_increase_factor;
end

%--------------------------------------------------

function x = diff2(x,k,dn)
%DIFF2   Difference function.  If X is a vector [x(1) x(2) ... x(n)],
%       then DIFF(X) returns a vector of central differences between
%       elements [x(3)-x(1)  x(4)-x(2) ... x(n)-x(n-2)].  If X is a
%	matrix, the differences are calculated down each column:
%       DIFF(X) = X(3:n,:) - X(1:n-2,:).
%	DIFF(X,n) is the n'th difference function.

%	J.N. Little 8-30-85
%	Copyright (c) 1985, 1986 by the MathWorks, Inc.

if nargin < 2,	k = 1; end
if nargin < 3,	dn = 2; end
for i=1:k
	[m,n] = size(x);
	if m == 1
                x = x(1+dn:n) - x(1:n-dn);
	else
                x = x(1+dn:m,:) - x(1:m-dn,:);
	end
end

