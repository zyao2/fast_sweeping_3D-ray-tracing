function TT=eikonal_time2d(W,TT,dy,dx)
[ny,nx]=size(W);
TT=fsweep_2d(W,TT,1,1,ny,1,1,nx,dy,dx);
TT=fsweep_2d(W,TT,ny,-1,1,1,1,nx,dy,dx);
TT=fsweep_2d(W,TT,1,1,ny,nx,-1,1,dy,dx);
TT=fsweep_2d(W,TT,ny,-1,1,nx,-1,1,dy,dx);

return;

NDims=2;
TT_old=TT;
[incs,inits,ends]=initial_fs(NDims);
for kk=1:12
    kk
    [incs,inits,ends]=setSweep(incs,inits,ends,[ny,nx]);
    TT=fsweep_old(W,TT_old,inits(1),incs(1),ends(1),inits(2),incs(2),ends(2),dy,dx);
    error=max(max(abs(TT-TT_old)))
    if(error<0.01)
        break;
    end
    TT_old=TT;
end
return;






