function u=fsweep_3d(W,u, incs,inits,ends)
global h t_inf
[NX,NY,NZ]=size(W);
for ix=inits(1):incs(1):ends(1)     
    for iy=inits(2):incs(2):ends(2) 
        for iz=inits(3):incs(3):ends(3)
            tt=u(ix,iy,iz);
            if (ix==1)
                t1=u(2,iy,iz);
            elseif(ix==NX)
                t1=u(NX-1,iy,iz);
            else
                t1=min(u(ix-1,iy,iz),u(ix+1,iy,iz));
            end
            
            if (iy==1)
                t2=u(ix,2,iz);
            elseif(iy==NY)
                t2=u(ix,NY-1,iz);
            else
                t2=min(u(ix,iy-1,iz),u(ix,iy+1,iz));
            end

            if (iz==1)
                t3=u(ix,iy,2);
            elseif(iz==NZ)
                t3=u(ix,iy,NZ-1);
            else
                t3=min(u(ix,iy,iz-1),u(ix,iy,iz+1));
            end 
            if(t1>t2)
                a=t1;
                t1=t2;
                t2=a;
            end
            if(t1>t3)
                a=t1;
                t1=t3;
                t3=a;
            end
            if(t2>t3)
                a=t2;
                t2=t3;
                t3=a;
            end 
            if(t1==t_inf || tt<t1 )
                continue;
            end
            hf=h/W(ix,iy,iz);
            ubar=t1+hf;
            if(ubar>t2)
                bb=2*hf*hf-(t1-t2)^2;
                ubar=0.5*(t1+t2+sqrt(bb));
               if(ubar>t3 && t3>t1+t2)
                   bb=2*(t1*t2-t1*t1-t2*t2+t1*t3+t2*t3-t3*t3);
                   ubar=(t1+t2+t3+sqrt(bb))/3+hf*hf;
               end
            end
            if(u(ix,iy,iz)>ubar)
                u(ix,iy,iz)=ubar;
            end
        end
    end
end

