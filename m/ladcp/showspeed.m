
tic
a = [];
for n=1:100000
  a = [a,1];
end
toc

tic
a = ones(1,100000);
for n=1:length(a(:))
  disp([int2str(n),' of ',int2str(length(a(:)))])
  a(n) = 2;
end
toc
