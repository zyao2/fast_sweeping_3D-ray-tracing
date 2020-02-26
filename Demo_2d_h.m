
clear all
%close all
%(ny,nx)
global t_inf
t_inf=99999;
n = 1000;
h=1;
Velocity = zeros(n); % constant velocity
Velocity(1:200,:)=1;
Velocity(201:300,:)=2;
Velocity(301:400,:)=10;
Velocity(401:500,:)=15;
Velocity(501:end,:)=20;
Velocity(100:200,700:800)=120;


v0=3;
vz=0.1;
for i=1:n
    v=v0+(i-1)*vz;
    Velocity(i,:)=v;
end
%Velocity(300:end,:)=0.01;
%Velocity(100:200,700:800)=120;


Velocity=1./Velocity;
%Velocity=Velocity*0+1;


isy=1;isx=1; %source location

%inital time map
T0 = zeros(n)+t_inf;
T0(isy,isx) = 0;

%near source
for i=-5:5
    k1=i+isy;
    if(k1<1)
        continue;
    end
    aa=(k1-1)*h;
    aa=aa*aa;
    for j=-5:5
        k2=j+isx;
        if(k2<1)
            continue;
        end
        bb=(k2-1)*h;
        bb=bb*bb;
        dis=sqrt(aa+bb);
        %T0(k1,k2)=dis*Velocity(k1,k2);
    end
end



tic;
TT_new=eikonal_time2d_h(Velocity,T0,h);

toc

TTT=TT_new;
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
x0=[isy;isx];
x1=[1;900];
ray = x1;
tic
tau = 5;
Geval = @(G,x)[interp2(1:n,1:n,G(:,:,1),x(2),x(1)); ...
             interp2(1:n,1:n,G(:,:,2),x(2),x(1)) ];
         
%get path
for kk=1:1
ray = x1;
for i=1:1.5*n/tau
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

h = plot(ray(2,:),ray(1,:)); set(h, 'LineWidth', 2);
h = plot(x0(2),x0(1), '.r'); set(h, 'MarkerSize', 25);
h = plot(x1(2),x1(1), '.b'); set(h, 'MarkerSize', 25);
axis ij;

