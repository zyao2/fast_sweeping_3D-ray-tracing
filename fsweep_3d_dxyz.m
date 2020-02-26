function u=fsweep_3d(W,u, incs,inits,ends)
global dxyz t_inf
[NX,NY,NZ]=size(W);
dt=zeros(3,1);
t=zeros(3,1);

for ix=inits(1):incs(1):ends(1)     
    for iy=inits(2):incs(2):ends(2) 
        for iz=inits(3):incs(3):ends(3)
            tt=u(ix,iy,iz);
            if (ix==1)
                t(1)=u(2,iy,iz);
            elseif(ix==NX)
                t(1)=u(NX-1,iy,iz);
            else
                t(1)=min(u(ix-1,iy,iz),u(ix+1,iy,iz));
            end
            
            if (iy==1)
                t(2)=u(ix,2,iz);
            elseif(iy==NY)
                t(2)=u(ix,NY-1,iz);
            else
                t(2)=min(u(ix,iy-1,iz),u(ix,iy+1,iz));
            end

            if (iz==1)
                t(3)=u(ix,iy,2);
            elseif(iz==NZ)
                t(3)=u(ix,iy,NZ-1);
            else
                t(3)=min(u(ix,iy,iz-1),u(ix,iy,iz+1));
            end 
            for k=1:3
             dt(k)=dxyz(k)/W(ix,iy,iz)*0;
            end
           indx=(1:3)';
            if(t(1)>t(2)+dt(2))
                a=t(1);
                k=indx(1);
                t(1)=t(2);
                indx(1)=indx(2);
                t(2)=a;
                indx(2)=k;
            end
            if(t(1)>t(3))
                a=t(1);
                k=indx(1);
                t(1)=t(3);
                indx(1)=indx(3);
                t(3)=a;
                indx(3)=k;
            end
            if(t(2)>t(3))
                a=t(2);
                k=indx(2);
                t(2)=t(3);
                indx(2)=indx(3);
                t(3)=a;
                indx(3)=k;
            end 
            if(t(1)==t_inf || tt<t(1) )
                continue;
            end
            
            hf=dxyz(indx(1))/W(ix,iy,iz);
            ubar=t(1)+hf;
            if(ubar>t(2))
                aa=dxyz(indx(1))^2+dxyz(indx(2))^2;
                bb=-2*(dxyz(indx(2))^2*t(1)+dxyz(indx(1))^2*t(2));
                cc=dxyz(indx(2))^2*t(1)^2+dxyz(indx(1))^2*t(2)^2-(dxyz(indx(1))*dxyz(indx(2))/W(ix,iy,iz))^2;
                bbb=bb^2-4*aa*cc;          
                ubar=(-bb+sqrt(bbb))/(2*aa);

               if(ubar>t(3)&& t(3)>t(1)+t(2))
                    aa=dxyz(indx(1))^2*dxyz(indx(3))^2+dxyz(indx(1))^2*dxyz(indx(2))^2+...
                        dxyz(indx(3))^2*dxyz(indx(2))^2;
                    bb=-2*(dxyz(indx(2))^2*dxyz(indx(3))^2*t(1)+dxyz(indx(1))^2*dxyz(indx(3))^2*t(2)+...
                        dxyz(indx(2))^2*dxyz(indx(1))^2*t(3));
                    cc=(dxyz(indx(2))*dxyz(indx(3))*t(1))^2+(dxyz(indx(1))*dxyz(indx(3))*t(2))^2+...
                        (dxyz(indx(1))*dxyz(indx(2))*t(3))^2;
                    cc=cc-(dxyz(indx(1))*dxyz(indx(2))*dxyz(indx(3))/W(ix,iy,iz))^2;
                    bbb=bb^2-4*aa*cc;
                   
                   ubar1=(-bb+sqrt(bbb))/(2*aa);
                                      
                   if(abs(ubar1-ubar)>1.e-4)
                       stop
                   end
                   
               end
            end
            if(u(ix,iy,iz)>ubar)
                u(ix,iy,iz)=ubar;
            end
        end
    end
end

