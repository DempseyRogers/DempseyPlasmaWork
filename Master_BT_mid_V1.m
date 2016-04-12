tic
clc
clear all
clf
close all

%Unitless Constants
rho=1;
eta=.18;  % .03-.18
nu=.76;

%Domain
R_IN=.5;
R_OUT=1.5;
Z_DOWN=0;
Z_UP=1;

%Spatial Step Sizes
dr=.01;
dz=.01;
dt=.00001 ;

%Controls
PoloidialRingAmp=13.1;
% G_AMP=6.58;
G_AMP=13.6;
Bz_AMP=PoloidialRingAmp;
Br_AMP=PoloidialRingAmp;
EXTERNAL_Const_Z=0;%13.6;%good


%otped
VrI= 0;
VpI= 0;%.27;%.27;
VzI= 0;%3.3333e-07;

%Tollerance
Conf_Toll=.5;   

%Time Parameters
SIM_DUR=20;
RESOLUTION_TIME=1;


% for i=2:100
% 
%     VzI(i)=VzI(i-1)+.01   ; 
%     
[time,  Bp ] = MidPt_MHD(R_IN,   R_OUT,   Z_DOWN,   Z_UP,   G_AMP,  ... Bz_AMP,   Br_AMP,  
    PoloidialRingAmp,PoloidialRingAmp,EXTERNAL_Const_Z, ...
    Conf_Toll,   dr,   dz,   dt,   SIM_DUR,   RESOLUTION_TIME,   rho,   eta,   nu, VpI,VpI,VzI);
% 
% VzI(i)
%                     
% percent=i
% if time(i)-time(i-1)<-5e4
% break
% end
% 
%  end
% vpa(VzI)
% vpa(time)

%constant Zfield 
%3.35 -->converged to hundreths for best time 


toc
% 
% plot(VzI(2:size(time,2)), time(2:size(time,2)))
% title('Vertical Velocity Vs. time')
% xlabel('Dimensionless')
% ylabel('Time (s)')
% 
% 
% close all
% for loops=1:RESOLUTION_TIME:size(Bp,3)
% clf
% subplot(1,2,1)
% contour(Bp(:,:,loops))
% axis([1 size(Bp,1) 1 size(Bp,2)])
% axis([(R_IN)/dr (R_OUT)/dr (R_IN)/dr (R_OUT)/dr])
% 
% subplot(1,2,2)
% caxis([bpmin,bpmax]);
% surface(Bp(:,:,loops))
% axis([1 size(Bp,1) 1 size(Bp,2)])
% axis('square')
% title(loops)
% drawnow 
% end


