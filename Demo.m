
clear all
%close all
global NDims h
global dimsize
NDims=3; h=1;
global t_inf
t_inf=99999.;
nz=100;nx=300;ny=300;
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
%Velocity=Velocity*0+1;
TT=zeros(size(Velocity))+t_inf;
isz=1; isx=10; isy=20;

TT(isz, isx, isy)=0;
dimsize=[nz,nx,ny];
tic
TT1=FSM(Velocity,TT);
toc

TTT=squeeze(TT1(:,isx,:));

figure;hold on;
contour(TTT,50)
%return
tic
options.order = 2;
G0 = grad(TTT, options);

%Normalize the gradient to obtained \(G(x) = G_0(x)/\norm{G_0(x)}\), 
%in order to have unit speed geodesic curve (parameterized by arc length).

G = G0 ./ repmat( sqrt( sum(G0.^2, 3) ), [1 1 2])+1.e-15;
 G(end,end)=0;

%Initialize the path with the ending point.
x0=[isz;isy];
x1=[1;150];
ray = x1;
tic
tau = 0.1;
[m,n]=size(TTT);
Geval = @(G,x)[interp2(1:m,1:n,G(:,:,1),x(2),x(1)); ...
             interp2(1:m,1:n,G(:,:,2),x(2),x(1)) ];
         
%get path
for kk=1:1
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
end
toc

hh = plot(ray(2,:),ray(1,:)); set(hh, 'LineWidth', 2);
hh = plot(x0(2),x0(1), '.r'); set(hh, 'MarkerSize', 25);
hh = plot(x1(2),x1(1), '.b'); set(hh, 'MarkerSize', 25);
axis ij;


