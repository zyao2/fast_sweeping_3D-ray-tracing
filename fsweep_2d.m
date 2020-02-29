function u=fsweep_2d(W,u,ny1,ndy,ny2,nx1,ndx,nx2,dy,dx)
global t_inf
[nodey,nodex]=size(W);
for i=ny1:ndy:ny2 %y     
    for j=nx1:ndx:nx2%column     
        if (i==1)
            uymin = u(2,j);
        elseif(i==nodey)
            uymin = u(nodey-1,j);
        else
            uymin = min(u(i-1,j),u(i+1,j));
        end
        if (j==1)
             uxmin = u(i,2);
        elseif(j==nodex)
             uxmin = u(i,nodex-1);
        else
            uxmin = min(u(i,j-1),u(i,j+1));
        end
        umin=min(uxmin,uymin);
        if(umin==t_inf || umin>u(i,j))
            continue;
        end
        if(uxmin<uymin)
            t2=uymin;
            ubar=uxmin+dx*W(i,j);
        else
            t2=uxmin;
            ubar=uymin+dy*W(i,j);
        end
        if(ubar>t2)
            a = (1/dy^2 + 1/dx^2);
            b = -2*(uymin/dy^2 + uxmin/dx^2);
            c = (uymin/dy)^2 + (uxmin/dx)^2 - (W(i,j))^2;                
            ubar=(-b+sqrt(b^2-4*a*c))/(2*a);
        end        
        if(u(i,j)>ubar)
            u(i,j)=ubar;
        end
    end
end
