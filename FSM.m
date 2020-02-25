function TT=FSM(Velocity,TT)
global dimsize

maxSweeps=8;
NDims=3;
[incs,inits,ends]=initial_fs(NDims);

for k=1:maxSweeps
    [incs,inits,ends]=setSweep(incs,inits,ends, dimsize);
    TT=fsweep_3d(Velocity,TT,incs,inits,ends);
end

return;

      
