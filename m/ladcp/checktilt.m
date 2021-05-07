function r=checktilt(rot,angle)
%
%  angle = [ rol1 ; pit1 ; rol2 ; pit2]
% rotate pit2,rol2 by rot and then compare results

good = find(~isnan(angle(3,:)+angle(4,:)+angle(1,:)+angle(2,:)));
[rol21,pit21] = uvrot(angle(3,good),angle(4,good),rot);

dr = rol21-angle(1,good);
dp = pit21-angle(2,good);

r = (dr-mean(dr)).^2 + (dp-mean(dp)).^2;
r = sum(r);
