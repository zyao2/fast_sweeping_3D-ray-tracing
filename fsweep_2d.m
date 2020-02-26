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
        sx=1;sy=1;
        if(uymin==t_inf)
            sy=0;
        end
        if(uxmin==t_inf)
            sx=0;
        end
        a = (sy/dy^2 + sx/dx^2);
        b = -2*(sy*uymin/dy^2 + sx*uxmin/dx^2);
        c = (sy*uymin/dy)^2 + (sx*uxmin/dx)^2 - (W(i,j))^2;
        bb=b^2-4*a*c;
        if (bb>0)                  
            ubar=(-b+sqrt(b^2-4*a*c))/(2*a);
        else
            if(uxmin<=uymin)
                ubar=uxmin+dx*W(i,j);
            else
                ubar = uymin+dy*W(i,j);
            end
        end
         
        if(u(i,j)>ubar)
            u(i,j)=ubar;
        end
    end
end
