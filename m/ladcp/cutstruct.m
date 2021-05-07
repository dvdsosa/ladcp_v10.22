function a=cutstruct(a,ii)
% function a=cutstruct(a,ii)
%
% reduce arrays within a structure
%
% This works only for fields with up to 3 dimensions !
%
% input  :  a              - structure 
%           ii             - vector containing 1 and 0
%                            All dimensions of fields of 'a' that
%                            have the same length as 'ii' will be
%                            reduced to the parts where 'ii' is 1.
%         
% output :  a              - reduced structure
%
% version 2  last change 19.06.2015

% G.Krahmann, GEOMAR
%
% bug fix for fields that have 3 dimensions

lz = length(ii);
iok = find(ii==1);
if isfield(a,'cutindx')
  a.cutindx = a.cutindx(1)-1+[iok(1) iok(end)];
else
  a.cutindx = [iok(1) iok(end)];
end
if isstruct(a)
  fnames = fieldnames(a);
  for n=1:size(fnames,1)
    dummy = getfield(a,fnames{n});
    [ly,lx,lzz]=size(dummy);
    if ly==lz
      if lzz>1
        a=setfield(a,fnames{n},dummy(iok,:,:));
      else
        a=setfield(a,fnames{n},dummy(iok,:));
      end
    elseif lx==lz
      if lzz>1
        a=setfield(a,fnames{n},dummy(:,iok,:));
      else
        a=setfield(a,fnames{n},dummy(:,iok));
      end
    end
  end
else
  error('first argument must be a structure')
end
