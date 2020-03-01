function f = grad_dxy(M, dy,dx)
fx = (M([2:end 1],:,:)-M)/(dy*dy);
fy = (M(:,[2:end 1],:)-M)/(dx*dx);
f = cat(3,fx,fy);
