
clear all
%close all
global NDims 
global dimsize dxyz t_inf
dxyz=[1,1,2];
NDims=3;
t_inf=99999;
nz=100;nx=100;ny=300;
Velocity=ones(nz,nx,ny);
v0=3;vz=0.1;
for k=1:nz
    vv=v0+(k-1)*vz;
    for i=1:nx
        for j=1:ny
            Velocity(k,i,j)=vv;
        end
    end
end
Velocity(1:30,:,:)=1;
Velocity(31:60,:,:)=4;
Velocity(61:end,:,:)=10;
%Velocity=Velocity*0+abs(rand(size(Velocity)))+1;
TT=zeros(size(Velocity))+t_inf;
isz=1; isx=10; isy=20;

TT(isz, isx, isy)=0;
for i=-5:5
    k1=i+isz;
    if(k1<1)
        continue;
    end
    aa=(k1-1)*dxyz(1);
    aa=aa*aa;
    for j=-5:5
        k2=j+isx;
        if(k2<1)
            continue;
        end
        bb=(k2-1)*dxyz(2);
        bb=bb*bb;
        for k=-5:5
            k3=j+isy;
            if(k3<1)
                continue;
            end
            cc=(k3-1)*dxyz(3);
            cc=cc*cc;  
            dis=sqrt(aa+bb+cc);
            %TT(k1,k2,k3)=dis*Velocity(k1,k2,k3);
        end
    end
end
dimsize=[nz,nx,ny];
tic
TT1=FSM_dxyz(Velocity,TT);
toc

TTT=squeeze(TT1(:,isx,:));

figure;hold on;
contour(TTT,50)
%return
tic
options.order = 1;

G0=grad(TTT, 1,2, options);
%Normalize the gradient to obtained \(G(x) = G_0(x)/\norm{G_0(x)}\), 
%in order to have unit speed geodesic curve (parameterized by arc length).

G = G0 ./ repmat( sqrt( sum(G0.^2, 3) ), [1 1 2])+1.e-15;
 G(end,end)=0;

%Initialize the path with the ending point.
x0=[isz;isy];
x1=[1;180];
ray = x1;
tic
tau = 0.1;
[m,n]=size(TTT);

         
%get path

ray = x1;
for i=1:5000
    ind=floor(ray(:,end));
    aa=[G(ind(1),ind(2),1);G(ind(1),ind(2),2)];
    aa=ray(:,end) - tau*aa;%Geval(G, ray(:,end));
    if(aa(1)<1)
        aa(1)=1;
    end
    if(aa(2)<1)
        aa(2)=1;
    end
    ray(:,end+1) = aa;%ray(:,end) - tau*Geval(G, ray(:,end));
    if norm(ray(:,end)-x0)<1
        break;
    end
end
ray(:,end+1) = x0;

toc

hh = plot(ray(2,:),ray(1,:)); set(hh, 'LineWidth', 2);
hh = plot(x0(2),x0(1), '.r'); set(hh, 'MarkerSize', 25);
hh = plot(x1(2),x1(1), '.b'); set(hh, 'MarkerSize', 25);
axis ij;