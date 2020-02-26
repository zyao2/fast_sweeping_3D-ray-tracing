function TT=eikonal_time2d_h(W,TT,h)
[ny,nx]=size(W);
TT=fsweep_2d_h(W,TT,1,1,ny,1,1,nx, h);
TT=fsweep_2d_h(W,TT,ny,-1,1,1,1,nx,h);
TT=fsweep_2d_h(W,TT,1,1,ny,nx,-1,1,h);
TT=fsweep_2d_h(W,TT,ny,-1,1,nx,-1,1,h);
return;

NDims=2;
TT_old=TT;
[incs,inits,ends]=initial_fs(NDims);
for kk=1:18
    kk
    [incs,inits,ends]=setSweep(incs,inits,ends,[ny,nx]);
    TT=fsweep_2d_h(W,TT_old,inits(1),incs(1),ends(1),inits(2),incs(2),ends(2),h);
    if(max(max(abs(TT-TT_old)))<0.01)
        break;
    end
    TT_old=TT;
end
return;




