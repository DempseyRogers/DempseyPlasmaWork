%%  First attempt at MHD system of equations.  Initially started with 
%   Eulers meathod, Implimenting MHDcurl, MHDlaplacian, MHDgradR, and 
%   MHDgradZ 
%   Dempsey Rogers 6-4-15
% clc
clear all
clf
% close all
% MHD EQUATIONS
%dB/dt=curl(VxB)+ata*Laplacian(B)
%               and
%dV/dt=Jxv+rho*nue*Laplacian(V)-(V dot Grad(v))
%               Where:
% J=curl(B)=current density
% rho= mass density                 p
% ata= resistive constant           n
% nue= viscous constant             v

% Constants
rho=1;
eta=1;
nu=1;
c=1;

% Stepsize
dr=.1;
dz=.1;
dt=.001;

tf = 1;

% Dimensions
r=100:dr:102;
z=0:dz:4;
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

%Br(:,:,1)=.13*sin(pi*Zmatrix) .*cos(pi) ;
%Bp(:,:,1)=zeros(size(Rmatrix)) ;
%Bz(:,:,1)=13.6*ones(size(Rmatrix))  ; %ITER 13.6T 

%Bp(:,:,1)=.01*sin(pi*Zmatrix) .*cos(pi) ;
%Br(:,:,1)=zeros(size(Rmatrix)) ;
%Bz(:,:,1)=4*ones(size(Rmatrix))  ; 

Bp(:,:,1)=.01*sin(pi*Zmatrix)    ;
Br(:,:,1)=zeros(size(Rmatrix))   ;
Bz(:,:,1)=4*ones(size(Rmatrix))  ;

Vr(:,:,1)=0*Rmatrix      ;
Vp(:,:,1)=0*Rmatrix      ;
Vz(:,:,1)=0*Rmatrix      ;
 
%for l = 1 : 1400
for l=1:size(t,2)
     
    [Jr, Jp, Jz]=MHDcurl(Br(:,:,l),Bp(:,:,l),Bz(:,:,l),r,z);
    [curlVBr, curlVBp, curlVBz]=MHDcurl(Vp(:,:,l).*Bz(:,:,l)-Vz(:,:,l).*Bp(:,:,l) ...
                                       ,Vz(:,:,l).*Br(:,:,l)-Vr(:,:,l).*Bz(:,:,l) ...
                                       ,Vr(:,:,l).*Bp(:,:,l)-Vp(:,:,l).*Br(:,:,l),r,z);
                                   
        Brmid=(curlVBr+eta*MHDlaplacian(Br(:,:, l), r, z))*dt/2+Br(:,:, l);
        Bpmid=(curlVBp+eta*MHDlaplacian(Bp(:,:, l), r, z))*dt/2+Bp(:,:, l);
        Bzmid=(curlVBz+eta*MHDlaplacian(Bz(:,:, l), r, z))*dt/2+Bz(:,:, l);
        Vrmid=((Jp.*Bz(:,:,l)-Jz.*Bp(:,:,l))/(c*rho) ...
                   +nu*MHDlaplacian(Vr(:,:,l), r, z)-Vr(:,:,l).*MHDgradR ...
                   (Vr(:,:,l),r)-Vz(:,:,l).*MHDgradZ(Vr(:,:,l),z))*dt/2+Vr(:,:,l);
        Vpmid=((Jz.*Br(:,:,l)-Jr.*Bz(:,:,l))/(c*rho) ...
                    +nu*MHDlaplacian(Vp(:,:,l), r, z)-Vr(:,:,l).*MHDgradR ...
                    (Vp(:,:,l),r)-Vz(:,:,l).*MHDgradZ(Vp(:,:,l),z))*dt/2+Vp(:,:,l);
        Vzmid=((Jr.*Bp(:,:,l)-Jp.*Br(:,:,l))/(c*rho) ...
                    +nu*MHDlaplacian(Vz(:,:,l), r, z)-Vr(:,:,l).*MHDgradR ...
                    (Vz(:,:,l),r)-Vz(:,:,l).*MHDgradZ(Vz(:,:,l),z))*dt/2+Vz(:,:,l);

    [Jr, Jp, Jz]=MHDcurl(Brmid,Bpmid,Bzmid,r,z);
    [curlVBr, curlVBp, curlVBz]=MHDcurl(Vpmid.*Bzmid-Vzmid.*Bpmid ...
                                       ,Vzmid.*Brmid-Vrmid.*Bzmid ...
                                       ,Vrmid.*Bpmid-Vpmid.*Brmid,r,z);
                                           
        Br(:,:,l+1)=(curlVBr+eta*MHDlaplacian(Brmid, r, z))*dt+Br(:,:, l);
        Bp(:,:,l+1)=(curlVBp+eta*MHDlaplacian(Bpmid, r, z))*dt+Bp(:,:, l);
        Bz(:,:,l+1)=(curlVBz+eta*MHDlaplacian(Bzmid, r, z))*dt+Bz(:,:, l);
        Vr(:,:,l+1)=((Jp.*Bzmid-Jz.*Bpmid)/(c*rho) ...
                   +nu*MHDlaplacian(Vrmid, r, z)-Vrmid.*MHDgradR ...
                   (Vrmid,r)-Vzmid.*MHDgradZ(Vrmid,z))*dt+Vr(:,:,l);
        Vp(:,:,l+1)=((Jz.*Brmid-Jr.*Bzmid)/(c*rho) ...
                    +nu*MHDlaplacian(Vpmid, r, z)-Vrmid.*MHDgradR ...
                    (Vpmid,r)-Vzmid.*MHDgradZ(Vpmid,z))*dt+Vp(:,:,l);
        Vz(:,:,l+1)=((Jr.*Bpmid-Jp.*Brmid)/(c*rho) ...
                    +nu*MHDlaplacian(Vzmid, r, z)-Vrmid.*MHDgradR ...
                    (Vzmid,r)-Vzmid.*MHDgradZ(Vzmid,z))*dt+Vz(:,:,l);
        
%Br(1,:,l+1) = zeros(1,20);
%Br(2,:,l+1) = zeros(1,20);
%Br(40,:,l+1) = zeros(1,20);
%Br(41,:,l+1) = zeros(1,20);
%Bp(2,:,l+1) = zeros(1,21);
%Bp(40,:,l+1) = zeros(1,21);

Bp(1,:,l+1) = zeros(1,21);
Bp(2,:,l+1) = zeros(1,21);
Bp(40,:,l+1) = zeros(1,21);
Bp(41,:,l+1) = zeros(1,21);
Bp(:,1,l+1) = zeros(1,41);
Bp(:,2,l+1) = zeros(1,41);
Bp(:,20,l+1) = zeros(1,41);
Bp(:,21,l+1) = zeros(1,41);

Br(1,:,l+1) = zeros(1,21);
Br(2,:,l+1) = zeros(1,21);
Br(40,:,l+1) = zeros(1,21);
Br(41,:,l+1) = zeros(1,21);
Br(:,1,l+1) = zeros(1,41);
Br(:,2,l+1) = zeros(1,41);
Br(:,20,l+1) = zeros(1,41);
Br(:,21,l+1) = zeros(1,41);

Bz(1,:,l+1) = 4*ones(1,21) ;
Bz(2,:,l+1) = 4*ones(1,21) ;
Bz(40,:,l+1) = 4*ones(1,21) ;
Bz(41,:,l+1) = 4*ones(1,21) ;
Bz(:,1,l+1) = 4*ones(41,1) ;
Bz(:,2,l+1) = 4*ones(41,1) ;
Bz(:,20,l+1) = 4*ones(41,1) ;
Bz(:,21,l+1) = 4*ones(41,1) ;

Vz(1,:,l+1) = zeros(1,21);
Vz(2,:,l+1) = zeros(1,21);
Vz(40,:,l+1) = zeros(1,21);
Vz(41,:,l+1) = zeros(1,21);
Vr(:,1,l+1) = zeros(1,41);
Vr(:,2,l+1) = zeros(1,41);
Vr(:,20,l+1) = zeros(1,41);
Vr(:,21,l+1) = zeros(1,41);




%Bz(:,1,l+1) = zeros(1,41);
%Br(2,:,l+1) = zeros(1,20);
%Bz(:,20,l+1) = zeros(1,41);
%Bz(:,20,l+1) = zeros(1,41);

%Bp(:,1,l+1) = zeros(1,41);
%Br(2,:,l+1) = zeros(1,20);
%Bz(:,20,l+1) = zeros(1,41);
%Bp(:,20,l+1) = zeros(1,41);


end
  surf(Bp(:,:,l+1))
% subplot(2,3,1)
%   surf(Vr(:,:,80))
% subplot(2,3,2)
%   surf(Vp(:,:,80))
% subplot(2,3,3)
%   surf(Vz(:,:,80))  
%   
% subplot(2,3,4)
%   surf(Br(:,:,80))
% subplot(2,3,5)
%   surf(Bp(:,:,80))
% subplot(2,3,6)
%   surf(Bz(:,:,80))  


% subplot(3,3,1)
%   surf(Vr(:,:,2))
%    title('Vr(1)')
% subplot(3,3,2)
%   surf(Vr(:,:,round(size(Vr,3)/2)))
%    title('Vr(mid)')
% subplot(3,3,3)
%   surf(Vr(:,:,size(Vr,3)))
%    title('Vr(end)')
% 
% subplot(3,3,4)
%   surf(Vp(:,:,2))
%    title('Vp(1)')
% subplot(3,3,5)
%   surf(Vp(:,:,round(size(Vr,3)/2)))
%    title('Vp(mid)')
% subplot(3,3,6)
%   surf(Vp(:,:,size(Vr,3)))
%    title('Vr(end)')
% 
% subplot(3,3,7)
%   surf(Vz(:,:,2))
%    title('Vz(1)')
% subplot(3,3,8)
%   surf(Vz(:,:,round(size(Vr,3)/2)))
%    title('Vz(mid)')
% subplot(3,3,9)
%   surf(Vz(:,:,size(Vr,3)))
%    title('Vz(end)')

