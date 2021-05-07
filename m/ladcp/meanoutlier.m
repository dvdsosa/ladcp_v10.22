function [mout,leftover]=meanoutlier(data,frac);
% function [mout,leftover]=meanoutlier(data,frac);
%
% GEOMAR SVN $Id: meanreduce.m 43 2014-06-27 15:16:28Z gkrahmann@geomar.de $
%
% average only the values that deviate least from the mean
% outliers are iteratively removed one by one and the mean is recalculated
% each time
%
% input  :  data          - data vector
%           frac          - fraction of values to be used
%                           i.e. with frac=0.6  the most deviating values
%                           will be removed until only 60% of the original
%                           number of values is left
%                           non finite numbers are not counted in the
%                           fraction determination !
%
% output :  mout          - mean of the remaining values
%           leftover      - index vector of the remaining values
%
% uses :
%
% version 1		last change 17.05.2015

% G.Krahmann, GEOMAR Kiel, May 2015

%
% look for non-finite values
%
leftover = find(isfinite(data));
startlength = length(leftover);

%
% iterate until we have only frac values left over
%
mout = mean(data(leftover));
nleftover = startlength;
while nleftover>frac*startlength
  dis = abs(data(leftover)-mout);
  [dummy,worst] = max(dis);
  if worst>1 & worst<nleftover
    leftover = leftover([1:worst-1,worst+1:nleftover]);
  elseif worst==1
    leftover = leftover(2:end);
  else
    leftover = leftover(1:end-1);
  end
  nleftover = nleftover-1;
  mout = mean(data(leftover));
end

