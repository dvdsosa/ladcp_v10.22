function []=saveres(data,dr,p,ps,f,values)
% function []=saveres(data,dr,p,ps,f,values)
%
% store LADCP result in RODB format
% and store some other info as mat file
%
% version 1	 last change 10.02.2017

% changed values stored as lat/lon            G.K. May 2007   0.2-->0.3
% modified version id                         GK  Jun 2007    0.3-->0.4
% change save command                         GK  Aug 2007    0.4-->0.5
% added some header info                      GK  Jun 2008    0.5-->0.6
% write target strength info into mat file    GK, 16.11.2012  0.6-->0.7
% allow target strength in ASCII file         GK, 22.06.2015  0.7-->0.8
% fix NaN in bottom file position; replace julian/gregorian with
% Matlab datenum/datevec                      GK, 10.02.2017  0.8-->1


%
% make sure there is an error field
%
if ~isfield(dr,'uerr')
 dr.uerr=dr.u*NaN;
end


%
% extract some target strength data and transmit voltage
%
ts.dn.ts = data.ts_edited(data.izd(3),:);
ts.dn.z = data.z-data.zd(3);
ts.dn.xmv = data.xmv(1,:);
if isfield(data,'ts_all_d')
  ts.all.ts = data.ts_all_d;
  ts.all.z = data.izm;
  ts.all.weight = data.weight;
  ts.all.ts_edited = data.ts_edited;
  ts.all.time = data.tim;
end
if values.up==1
  ts.up.ts = data.ts_edited(data.izu(3),:);
  ts.up.z = data.z+data.zu(3);
  ts.up.xmv = data.xmv(2,:);
end
[zmin,indmin] = nmin(ts.dn.z);
zstep = diff(dr.z(1:2));
ea_down = dr.z*nan;
ea_up = dr.z*nan;
tim_down = dr.z*nan;
tim_up = dr.z*nan;
[a,b,c,d,e,f1] = datevec(data.tim);
data.datenum = datenum( [a',b',c',d',e',f1'] );
for n= 1:length(dr.z)
  ind = find( abs(dr.z(n)+ts.dn.z(1:indmin))<=zstep/2 );
  ea_down(n) = nmean(ts.dn.ts(ind));
  tim_down(n) = nmean(data.datenum(ind));
  ind = find( abs(dr.z(n)+ts.dn.z(indmin:end))<=zstep/2 );
  ea_up(n) = nmean(ts.dn.ts(ind+indmin-1));
  tim_up(n) = nmean(data.datenum(ind+indmin-1));
end


%
% store some results as a MAT file
%
save6([f.res,'.mat'],'ts','dr','p','ps','f')


%
% open file
%
fid = fopen([f.res,'.lad'],'wt');


%
% write file header and data
%
fprintf(fid,['Filename    = %s\n'],f.res);
fprintf(fid,['Date        = %s\n'],datestr(p.time_start,26));
fprintf(fid,['Start_Time  = %s\n'],datestr(p.time_start,13));
fprintf(fid,['Latitude    = %s\n'],num2str(values.lat));
fprintf(fid,['Longitude   = %s\n'],num2str(values.lon));
fprintf(fid,['Deviation   = %f\n'],values.magdev);
fprintf(fid,['Version     = %s\n'],p.software);
fprintf(fid,['Processed   = %s\n'],datestr(clock));
if p.save_ascii_target_strength==1
  fprintf(fid,['Units       = m:m/s:m/s:m/s:dB:days:dB:days\n'],[]);
  fprintf(fid,['Columns     = z:u:v:ev:down_ts:down_time:up_ts:up_time\n'],[]);
  fprintf(fid,['%6.1f %6.3f %6.3f %6.3f %5.1f %12.5f %5.1f %12.5f \n'],[dr.z,dr.u,dr.v,dr.uerr,ea_down,tim_down,ea_up,tim_up]');
else
  fprintf(fid,['Units       = m:m/s:m/s:m/s\n'],[]);
  fprintf(fid,['Columns     = z:u:v:ev\n'],[]);
  fprintf(fid,['%6.1f %6.3f %6.3f %6.3f \n'],[dr.z,dr.u,dr.v,dr.uerr]');
end


%
% close file
%
fclose(fid);


%
% in case we have bottom track data, store that in another file
%
if isfield(dr,'ubot')

  % save bottom track data
  % open file
  fid = fopen([f.res,'.bot'],'wt');

  fprintf(fid,['Filename    = %s\n'],f.res);
  fprintf(fid,['Date        = %s\n'],datestr(p.time_start,26));
  fprintf(fid,['Start_Time  = %s\n'],datestr(p.time_start,13));
  fprintf(fid,['Latitude    = %s\n'],num2str(values.lat));
  fprintf(fid,['Longitude   = %s\n'],num2str(values.lon));
  fprintf(fid,['Deviation   = %f\n'],values.magdev);
  fprintf(fid,['Bottom depth= %d\n'],fix(p.zbottom));
  fprintf(fid,['Version     = %s\n'],p.software);
  fprintf(fid,['Processed   = %s\n'],datestr(clock));
  fprintf(fid,['Units       = m:m/s:m/s:m/s\n'],[]);
  fprintf(fid,['Columns     = z:u:v:err\n'],[]);
  fprintf(fid,['%6.1f %6.3f %6.3f %6.3f\n'],...
            [dr.zbot,dr.ubot,dr.vbot,dr.uerrbot]');

  fclose(fid);

end
