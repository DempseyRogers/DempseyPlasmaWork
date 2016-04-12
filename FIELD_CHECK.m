clc
clear all
clf 
% Proposed Poloidal Field


[x,y]=meshgrid(1:.05:1.5,0:.05:1);

conx=2.5;
cony=.5;
dy=2*x-conx;
dx=-y+cony;

quiver(x,y,dx,dy);

j=ones(size(x)).*1.25;
k=ones(size(x)).*.5;
hold on 
plot(x,k,'r*') 
plot(j,y,'r*')
axis('square')
%%
%Curren Poloidal Field
clc
clear all
clf

[x,y]=meshgrid(1:.05:1.5,0:.05:1);

conx=2.5;
cony=.5;
dy=cos(2*pi*(x-.5));
dx=cos(pi*(y));

quiver(x,y,dx,dy);

j=ones(size(x)).*1.25;
k=ones(size(x)).*.5;
hold on 
plot(x,k,'r*') 
plot(j,y,'r*')
axis('square')



%%


clc
clear all
clf
% close all



syms x 
a=1;b=1;c=1;

f=a*exp(-(x-b)^2/(2*c^2));
ezplot(f)

%%
clc
clf
clear all
syms R Z


rho=1;
eta=.1;   
nu=.1;

%Domain
R_IN=1;
R_OUT=1.5;
Z_DOWN=0;
Z_UP=1;

%Spatial Step Sizes
dr=.01;
dz=.01;
dt=.000001 ;

%Controls
PoloidialRingAmp=1%12.18;%good
G_AMP=10;
Bz_AMP=PoloidialRingAmp;
Br_AMP=PoloidialRingAmp;
EXTERNAL_Const_Z=2.69;%good
VrI=0;
VpI=.27;
VzI=0;

%Tollerance
Conf_Toll=.7;   

%Time Parameters
SIM_DUR=20;
RESOLUTION_TIME=1;

% eta=1e-10;
EXTERNAL_Const_Zt(1)=EXTERNAL_Const_Z;


Bp=G_AMP*exp(-(R-R_IN-(R_OUT-R_IN)/2).^2/(.25^2)).*exp(-(Z-Z_DOWN-(Z_UP-Z_DOWN)/2).^2/(.25^2)) ;

ezsurf(Bp)




