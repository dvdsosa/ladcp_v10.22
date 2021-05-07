function [] = plot_controls(fig)
% function [] = plot_controls(fig)
%
% reloads the stored figure into the display window
%
% input  :  fig             - figure number
%
% version 2  last change 13.07.2013

% G.Krahmann, GEOMAR

% added clf                                 GK, 13.07.2013  1-->2

global mh

sfigure(2);
clf
figload(['tmp/',int2str(fig),'.fig'])

sfigure(1);
for n=1:length(mh)
  if mh(n)~=0
    set(mh(n),'foregroundcolor',[0,0,0]);
  end
end
set(mh(fig),'foregroundcolor',[1,0,0]);
