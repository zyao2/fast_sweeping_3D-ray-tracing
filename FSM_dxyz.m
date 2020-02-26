function TT=FSM_dxyz(Velocity,TT)
global dimsize 
maxSweeps=8;
NDims=3;

[incs,inits,ends]=initial_fs(NDims);
for k=1:maxSweeps-1
    [incs,inits,ends]=setSweep(incs,inits,ends, dimsize);
end
for k=1:maxSweeps
    [incs,inits,ends]=setSweep(incs,inits,ends, dimsize);
    TT=fsweep_3d_dxyz(Velocity,TT,incs,inits,ends);
end
      
