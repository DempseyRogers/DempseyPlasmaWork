clc
clear all
clf
close all

%Unitless Constants
rho=1;
eta=.1;  % .03-.18
nu=.1;
c=1;
%Domain
R_IN=100;
R_OUT=102;
Z_DOWN=100;
Z_UP=101;

%Spatial Step Sizes
dr=.1;
dz=.1;
dt=.0001 ;
tic
tf = 2;

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
Bp(:,:,1)=5 *ones(size(Rmatrix))   ;
Br(:,:,1)=zeros(size(Rmatrix))   ;
Bz(:,:,1)=0*ones(size(Rmatrix))  ;

Vr(:,:,1)=0*Rmatrix      ;
Vp(:,:,1)=0*Rmatrix      ;
Vz(:,:,1)=0*Rmatrix      ;
 
blz2=size(Bp,2);
blz1=size(Bp,1);
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

%                 Bp(:,2,l+1)=Bp(:,3,l+1);
%                 Bp(:,1,l+1)=Bp(:,3,l+1);


                 Bz(:,1:2,:)=4;
                 Bz(1:2,:,:)=4;
                 Bz(blz1-1:blz1,:,:)=4;
                 Bz(:,blz2-1:blz2,:)=4;
%                 Bp(:,blz2-1,l+1)=Bp(:,blz2-2,l+1);
%                 Bp(:,blz2,l+1)=Bp(:,blz2-2,l+1);
                 Br(:,1:2,:)=0;
                 Br(1:2,:,:)=0;
                 Br(blz1-1:blz1,:,:)=0;
                 Br(:,blz2-1:blz2,:)=0;
                
%                 Bp(2,:,l+1)=Bp(3,:,l+1);
%                 Bp(1,:,l+1)=Bp(3,:,l+1);
%                 Bp(:,2,l+1)=Bp(:,3,l+1);
%                 Bp(:,1,l+1)=Bp(:,3,l+1);
%                 Bz(:,2,l+1)=Bz(:,3,l+1);
%                 Bz(:,1,l+1)=Bz(:,3,l+1);
%                 Br(:,2,l+1)=Br(:,3,l+1);
%                 Br(:,1,l+1)=Br(:,3,l+1);                
                
%                 Bp(blz1-1,:,l+1)=Bp(blz1-2,:,l+1);
%                 Bp(blz1,:,l+1)=Bp(blz1-2,:,l+1);
%                 Bp(:,blz2-1,l+1)=Bp(:,blz2-2,l+1);
%                 Bp(:,blz2,l+1)=Bp(:,blz2-2,l+1);
%                 Bz(:,blz2-1,l+1)=Bz(:,blz2-2,l+1);
%                 Bz(:,blz2,l+1)=Bz(:,blz2-2,l+1);                
%                 Br(:,blz2-1,l+1)=Br(:,blz2-2,l+1);
%                 Br(:,blz2,l+1)=Br(:,blz2-2,l+1);                 
                 
                 
end







  close all

  FIVES=10;
  bpmax=max(max(max(Bp)));
  bpmin=min(min(min(Bp)));
  BaseName='ALFVEN_rz';
for loops=1:FIVES:size(Bp,3)
    if loops == 1
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
    end 
titlename=['Bp(:,:,',num2str(loops),')'];  
subplot(1,2,1)
contour(Bp(:,:,loops))
axis([0 ((R_OUT-R_IN)/dr+1) 0 (Z_UP-Z_DOWN)/dz+1 ])
title(titlename)
xlabel('R (.01 (.5m-1.5m))');
ylabel('Z (.01 (0m-1m))');
subplot(1,2,2)
% caxis([bpmin,bpmax]);
surf(Bp(:,:,loops))
axis([0 ((R_OUT-R_IN)/dr+1) 0 (Z_UP-Z_DOWN)/dz+1 bpmin bpmax ])
%    axis([0 (z2/dz+1) 0 z2/dz+1 -10^-2 10^-2 ])
 
title(titlename)
xlabel('R (.01 (.5,-1.5m))');
ylabel('Z (.01 (0m-1m))');
zlabel('Br Magnitutde (T)');
%   axis([z1/dz+1,z2/dz+1,r1/dr+1, r2/dr+1])
%   filename=eval(sprintf('BP_GAUSIAN_TEST',loops));
   drawnow 
%     FileName=[BaseName,num2str(loops),'.png'];
%    saveas(gcf,FileName);
   
%       pause(1/60)
clf
%      pause(.01)
end


























