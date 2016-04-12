function [Confinement_Time,  Bp ] = MidPt_MHD(R_IN,   R_OUT,   Z_DOWN,   Z_UP,   G_AMP,   Bz_AMP,   Br_AMP,   EXTERNAL_Const_Z, ...
                        Conf_Toll,   dr,   dz,   dt,   SIM_DUR,   RESOLUTION_TIME,   rho,   eta,   nu, VrI,VpI,VzI)
%dt MUST BE 1/100 SPATIAL 
%   Detailed explanation goes here
% Constants
% rho=1;
% eta=.1;
% nu=.1;
 c=1;

tic
% MHD EQUATIONS
%dB/dt=curl(VxB)+ata*Laplacian(B)
%               and
%dV/dt=Jxv+rho*nue*Laplacian(V)-(V dot Grad(v))
%               Where:
% J=curl(B)=current density
% rho= mass density                 p
% ata= resistive constant           n
% nue= viscous constant             v



tf = SIM_DUR;

% Dimensions
r1=R_IN;
r2=R_OUT;
z1=Z_DOWN;
z2=Z_UP;
r=r1:dr:r2;
z=z1:dz:z2;
t=0:dt:tf;

% Initial
Rmatrix=zeros(size(z,2),size(r,2));
Zmatrix=zeros(size(z,2),size(r,2));


for i= 1:size(z,2)
   for j= 1:size(r,2)    
       Rmatrix(i,j)=r(j);
       Zmatrix(i,j)=z(i);
   end
end
Bp(:,:,1)=G_AMP*exp(-(Rmatrix-R_IN-(R_OUT-R_IN)/2).^2/(.25^2)).*exp(-(Zmatrix-Z_DOWN-(Z_UP-Z_DOWN)/2).^2/(.25^2)) ;
 Bz(:,:,1)=Bz_AMP*cos(2*pi*(Rmatrix-.5)) ;
 Br(:,:,1)=Br_AMP*cos(pi*Zmatrix) ; %ITER 13.6T 

Vr(:,:,1)=VrI*Rmatrix      ;
Vp(:,:,1)=VpI*Rmatrix      ;
Vz(:,:,1)=VzI*Rmatrix      ;
 
constant=EXTERNAL_Const_Z;

initial_max_Bp=max(max(Bp(:,:,1)));
 
counter=0;
for l=1:size(t,2)
    
    
    counter=counter+1;
     
    [Jr, Jp, Jz]=MHDcurl_V2(Br(:,:,l),Bp(:,:,l),Bz(:,:,l),r,z);
    [curlVBr, curlVBp, curlVBz]=MHDcurl_V2(Vp(:,:,l).*Bz(:,:,l)-Vz(:,:,l).*Bp(:,:,l) ...
                                       ,Vz(:,:,l).*Br(:,:,l)-Vr(:,:,l).*Bz(:,:,l) ...
                                       ,Vr(:,:,l).*Bp(:,:,l)-Vp(:,:,l).*Br(:,:,l),r,z);
                                   
        Brmid=(curlVBr+eta*MHDlaplacian(Br(:,:, l), r, z))*dt/2+Br(:,:, l);
        Bpmid=(curlVBp+eta*MHDlaplacian(Bp(:,:, l), r, z))*dt/2+Bp(:,:, l);
        Bzmid=(curlVBz+eta*MHDlaplacian(Bz(:,:, l), r, z))*dt/2+Bz(:,:, l);
        Vrmid=((Jp.*Bz(:,:,l)-Jz.*Bp(:,:,l))/(c*rho) ...
                   +nu*MHDlaplacian(Vr(:,:,l), r, z)-Vr(:,:,l).*MHDgradR_V2 ...
                   (Vr(:,:,l),r)-Vz(:,:,l).*MHDgradZ_V2(Vr(:,:,l),z))*dt/2+Vr(:,:,l);
        Vpmid=((Jz.*Br(:,:,l)-Jr.*Bz(:,:,l))/(c*rho) ...
                    +nu*MHDlaplacian(Vp(:,:,l), r, z)-Vr(:,:,l).*MHDgradR_V2 ...
                    (Vp(:,:,l),r)-Vz(:,:,l).*MHDgradZ_V2(Vp(:,:,l),z))*dt/2+Vp(:,:,l);
        Vzmid=((Jr.*Bp(:,:,l)-Jp.*Br(:,:,l))/(c*rho) ...
                    +nu*MHDlaplacian(Vz(:,:,l), r, z)-Vr(:,:,l).*MHDgradR_V2 ...
                    (Vz(:,:,l),r)-Vz(:,:,l).*MHDgradZ_V2(Vz(:,:,l),z))*dt/2+Vz(:,:,l);

    [Jr, Jp, Jz]=MHDcurl_V2(Brmid,Bpmid,Bzmid,r,z);
    [curlVBr, curlVBp, curlVBz]=MHDcurl_V2(Vpmid.*Bzmid-Vzmid.*Bpmid ...
                                       ,Vzmid.*Brmid-Vrmid.*Bzmid ...
                                       ,Vrmid.*Bpmid-Vpmid.*Brmid,r,z);
                                           
        Br(:,:,l+1)=(curlVBr+eta*MHDlaplacian(Brmid, r, z))*dt+Br(:,:, l);
        Bp(:,:,l+1)=(curlVBp+eta*MHDlaplacian(Bpmid, r, z))*dt+Bp(:,:, l);
        Bz(:,:,l+1)=(curlVBz+eta*MHDlaplacian(Bzmid, r, z))*dt+Bz(:,:, l);
        Vr(:,:,l+1)=((Jp.*Bzmid-Jz.*Bpmid)/(c*rho) ...
                   +nu*MHDlaplacian(Vrmid, r, z)-Vrmid.*MHDgradR_V2 ...
                   (Vrmid,r)-Vzmid.*MHDgradZ_V2(Vrmid,z))*dt+Vr(:,:,l);
        Vp(:,:,l+1)=((Jz.*Brmid-Jr.*Bzmid)/(c*rho) ...
                    +nu*MHDlaplacian(Vpmid, r, z)-Vrmid.*MHDgradR_V2 ...
                    (Vpmid,r)-Vzmid.*MHDgradZ_V2(Vpmid,z))*dt+Vp(:,:,l);
        Vz(:,:,l+1)=((Jr.*Bpmid-Jp.*Brmid)/(c*rho) ...
                    +nu*MHDlaplacian(Vzmid, r, z)-Vrmid.*MHDgradR_V2 ...
                    (Vzmid,r)-Vzmid.*MHDgradZ_V2(Vzmid,z))*dt+Vz(:,:,l);

                
                %Boundary Conditions
                Bz(:,1:2,l+1)=constant;
                Bp(:,1:2,l+1)=0;
                Br(:,1:2,l+1)=0;
                
                
                for i=1:size(Br,1)
                   for j=3:size(Br,2)
                    if( (Rmatrix(i,j)-1)^2+(Zmatrix(i,j)-.5)^2 > (R_OUT-R_IN)/2)
                        
                Bz(i,j,l+1)=constant;
                Bp(i,j,l+1)=0;
                Br(i,j,l+1)=0; 
                
%                     elseif ( (Rmatrix(i,j)-1)^2+(Zmatrix(i,j)-.5)^2 > dr+(R_OUT-R_IN)/2)
%                 
%                 Bz(i,j,:)=NaN;
%                 Bp(i,j,:)=NaN;
%                 Br(i,j,:)=NaN; 
%                 
%                 Vr(:,:,1)=NaN;
%                 Vp(:,:,1)=NaN;
%                 Vz(:,:,1)=NaN;

                     end
                   end 
                end
                 if max(max(Bp(:,:,l)))<initial_max_Bp*Conf_Toll
                    break
                 end
                
                 if rem(counter,1000)== 0
                     counter
                 end
                 if rem(counter, 10000)== 0
                     clc
                 end
        
                 if counter > 30000
                     break  

                 end
                 
                 
                 
end
bpmax= max(max(max(Bp(:,:,:))));
bpmin=min(min(min(Bp(:,:,:))));


Confinement_Time=l*dt;


% for loops=1:RESOLUTION_TIME:size(Bp,3)
% clf
% subplot(1,2,1)
% contour(Bp(:,:,loops))
% axis([0 (1/dz+1)/2 0 1/dz+1 ])
% 
% subplot(1,2,2)
% caxis([bpmin,bpmax]);
% surface(Bp(:,:,loops))
%  axis([0 (z2/dz+1)/2 0 z2/dz+1 ])
% % axis('square')
% title(loops)
% drawnow 
% end



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


titlenumber=1

  close all
  FIVES=10;
  bpmax=max(max(max(Bp)));
  bpmin=min(min(min(Bp)));
  BaseName='NO_CONSTANT_Z';
for loops=1:FIVES:size(Bp,3)
    if loops == 1
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    end 
titlename=['Bp(T)'];  
subplot(1,2,1)
contour(Bp(:,:,loops))
axis([0 (Z_UP/dz+1) 0 Z_UP/dz+1 ])
title(titlename)
xlabel('Radial');
ylabel('Vertical');
subplot(1,2,2)
% caxis([bpmin,bpmax]);
surface(Bp(:,:,loops))
axis([0 (Z_UP/dz+1) 0 Z_UP/dz+1 bpmin bpmax ])
%    axis([0 (z2/dz+1) 0 z2/dz+1 -10^-2 10^-2 ])
 
title(titlename)
xlabel('Radial');
ylabel('Vertical');
zlabel('Bp (T)');
%   axis([z1/dz+1,z2/dz+1,r1/dr+1, r2/dr+1])
%   filename=eval(sprintf('BP_GAUSIAN_TEST',loops));
   drawnow 
    FileName=[BaseName,num2str(titlenumber),'.png'];
   saveas(gcf,FileName);
   titlenumber=titlenumber+1;
%       pause(1/60)
clf
%      pause(.01)
end



























end
