clear all
%close all
%(y,x)
global t_inf
dx=3;dy=1;
t_inf=99999;
m=1000;
n = 500;
Velocity = zeros(m,n); % constant velocity

Velocity(1:end,:)=1;
for i=1:n
    k=floor((i-1)*0.2)+200;
    Velocity(k:end,i)=2;
end

for i=1:n
    k=floor((i-1)*0.5)+400;
    Velocity(k:end,i)=4;
end

v0=3;
vz=0.1;
for i=1:m
    v=v0+(i-1)*vz;
   %Velocity(i,:)=v;
end
%Velocity(300:end,:)=0.01;
%Velocity(100:200,350:400)=120;


Velocity=1./Velocity;
%Velocity=Velocity*0+1;

sy=0;sx=71; %source location
ry=0;rx=1207; %reciver location
%inital time map
T0 = zeros(m,n)+t_inf;
%source location
isy=floor(sy/dy);
isx=floor(sx/dx);
%near source time approcimation
for i=-6:6
    k1=i+isy;
    if(k1<0)
        continue;
    end
    aa=(k1)*dy-sy;
    aa=aa*aa;
    for j=-5:5
        k2=j+isx;
        if(k2<0)
            continue;
        end
        bb=(k2)*dx-sx;
         bb=bb*bb;
        dis=sqrt(aa+bb);
        T0(k1+1,k2+1)=dis*Velocity(k1+1,k2+1);
    end
end


tic;
TT_new=eikonal_time2d(Velocity,T0,dy,dx);
%TT_new=TT_new(1:2:end,:);dx=2*dx;
toc
TTT=TT_new;
figure;hold on;
contour(TT_new,50)
%return
tic
G0 = grad_dxy(TTT, dy,dx);
% note that the grid size has bee scaled to one in grad_dxy

%Normalize the gradient to obtained \(G(x) = G_0(x)/\norm{G_0(x)}\), 
%in order to have unit speed geodesic curve (parameterized by arc length).
G = G0 ./ repmat( sqrt( sum(G0.^2, 3) ), [1 1 2]);

%Initialize the path with the ending point. !! noticed that the gradient
%has bee scaled with grid size
x0=[sy/dy;sx/dx];
x1=[ry/dy;rx/dx];
ray = x1;
tic
tau = 1;
Geval = @(G,x)[interp2(1:n,1:n,G(:,:,1),x(2),x(1)); ...
             interp2(1:n,1:n,G(:,:,2),x(2),x(1)) ];
         
%get path
for kk=1:1
ray = x1;
icon=0;
for i=1:5*n/tau
    ind=floor(ray(:,end));
    ind=ind+1;    % noticed that index (1:m,1:n)
    aa=[G(ind(1),ind(2),1);G(ind(1),ind(2),2)];
    aa=ray(:,end) - tau*aa;%Geval(G, ray(:,end)); 
    ray(:,end+1) = aa;%ray(:,end) - tau*Geval(G, ray(:,end));
    if (ray(1,end)<0 ||ray(2,end)<0)
        break;
    end
end
ray(:,end) = x0;
end
%ray=ray-1;
toc

h = plot(ray(2,:),ray(1,:)); set(h, 'LineWidth', 2);
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;
