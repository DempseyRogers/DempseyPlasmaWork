%%  First attempt at MHD system of equations.  Initially started with 
%   Eulers meathod, Implimenting MHDcurl, MHDlaplacian, MHDgradR, and 
%   MHDgradZ 
%   Dempsey Rogers 6-4-15
 clc
clear all
clf
 close all
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

% Constants
rho=1;
eta=.1;
nu=.1;
c=1;

% Stepsize
dr=.01;
dz=.01;
dt=.0001;

tf = 2;

% Dimensions
r1=1;
r2=1.5;
z1=0;
z2=1;
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

%Br(:,:,1)=.13*sin(pi*Zmatrix) .*cos(pi) ;
%Bp(:,:,1)=zeros(size(Rmatrix)) ;
%Bz(:,:,1)=13.6*ones(size(Rmatrix))  ; %ITER 13.6T 

Bp(:,:,1)=13.6*exp(-(Rmatrix-1.25).^2/.25^2).*exp(-(Zmatrix-.5).^2/.25^2) ;
Bz(:,:,1)=13.6*cos(pi*Rmatrix) ;
Br(:,:,1)=13.6*cos(pi*Zmatrix) ; %ITER 13.6T 

Vr(:,:,1)=0*Rmatrix      ;
Vp(:,:,1)=0*Rmatrix      ;
Vz(:,:,1)=0*Rmatrix      ;
 
constant=1;


%for l = 1 : 1400
for l=1:size(t,2)
    
    
    
    
    
     
    [Jr, Jp, Jz]=MHDcurl_V2(Br(:,:,l),Bp(:,:,l),Bz(:,:,l),r,z);
    [curlVBr, curlVBp, curlVBz]=MHDcurl_V2(Vp(:,:,l).*Bz(:,:,l)-Vz(:,:,l).*Bp(:,:,l) ...
                                       ,Vz(:,:,l).*Br(:,:,l)-Vr(:,:,l).*Bz(:,:,l) ...
                                       ,Vr(:,:,l).*Bp(:,:,l)-Vp(:,:,l).*Br(:,:,l),r,z);
    
        Br(:,:,l+1)=(curlVBr+eta*MHDlaplacian(Br(:,:, l), r, z))*dt+Br(:,:, l);
        Bp(:,:,l+1)=(curlVBp+eta*MHDlaplacian(Bp(:,:, l), r, z))*dt+Bp(:,:, l);
        Bz(:,:,l+1)=(curlVBz+eta*MHDlaplacian(Bz(:,:, l), r, z))*dt+Bz(:,:, l);
        Vr(:,:,l+1)=((Jp.*Bz(:,:,l)-Jz.*Bp(:,:,l))/(c*rho) ...
                   +nu*MHDlaplacian(Vr(:,:,l), r, z)-Vr(:,:,l).*MHDgradR_V2 ...
                   (Vr(:,:,l),r)-Vz(:,:,l).*MHDgradZ_V2(Vr(:,:,l),z))*dt+Vr(:,:,l);
        Vp(:,:,l+1)=((Jz.*Br(:,:,l)-Jr.*Bz(:,:,l))/(c*rho) ...
                    +nu*MHDlaplacian(Vp(:,:,l), r, z)-Vr(:,:,l).*MHDgradR_V2 ...
                    (Vp(:,:,l),r)-Vz(:,:,l).*MHDgradZ_V2(Vp(:,:,l),z))*dt+Vp(:,:,l);
        Vz(:,:,l+1)=((Jr.*Bp(:,:,l)-Jp.*Br(:,:,l))/(c*rho) ...
                    +nu*MHDlaplacian(Vz(:,:,l), r, z)-Vr(:,:,l).*MHDgradR_V2 ...
                    (Vz(:,:,l),r)-Vz(:,:,l).*MHDgradZ_V2(Vz(:,:,l),z))*dt+Vz(:,:,l);
        
                
                
%                 Bz(:,1:2,l+1)=constant;
%                 Bp(:,1:2,l+1)=0;
%                 Br(:,1:2,l+1)=0;
%                 
%                 
%                 for i=1:size(Br,1)
%                    for j=3:size(Br,2)
%                     if( (Rmatrix(i,j)-1)^2+(Zmatrix(i,j)-.5)^2 > .25)
%                         
%                 Bz(i,j,l+1)=constant;
%                 Bp(i,j,l+1)=0;
%                 Br(i,j,l+1)=0;                      
%                         
%                         
%                     end
%                    end 
%                 end
%                 
                
                
                
%Br(1,:,l+1) = zeros(1,20);
%Br(2,:,l+1) = zeros(1,20);
%Br(40,:,l+1) = zeros(1,20);
%Br(41,:,l+1) = zeros(1,20);
 
% Bp(1,:,l+1) = zeros(1,size(Bp,2));
% Bp(41,:,l+1) = zeros(1,size(Bp,2));
%Br(2,:,l+1) = zeros(1,20);
%Bp(40,:,l+1) = zeros(1,20);


%Bz(:,1,l+1) = zeros(1,41);
%Br(2,:,l+1) = zeros(1,20);
%Bz(:,20,l+1) = zeros(1,41);
%Bz(:,20,l+1) = zeros(1,41);

%Bp(:,1,l+1) = zeros(1,41);
%Br(2,:,l+1) = zeros(1,20);
%Bz(:,20,l+1) = zeros(1,41);
%Bp(:,20,l+1) = zeros(1,41);


end


% for loops=1:size(Bp,3)
%   plot(Bp(:,10,loops))
%   axis([0 60 -8e-3 8e-3])
%   title(loops)
%   drawnow
%   
% end
toc

% BaseName='BP_GAUSIAN_TEST_';
% close all
% FIVES=100;
% for loops=1:size(Bp,3)
%     subplot(1,2,1)
%    contour(Bp(:,:,loops))
%     axis([0 (z2/dz+1)/2 0 z2/dz+1 ])
% subplot(1,2,2)
%    surface(Bp(:,:,loops))
% %  axis([0 (z2/dz+1) 0 z2/dz+1 0e-3 10e-3])
%    axis([0 (z2/dz+1)/2 0 z2/dz+1 ])
%    title(loops)
% %  axis([z1/dz+1,z2/dz+1,r1/dr+1, r2/dr+1])
% %  filename=eval(sprintf('BP_GAUSIAN_TEST',loops));
%    drawnow 
% 
%    FIVES=FIVES+1;
%    if FIVES==100
%        drawnow 
%    FileName=[BaseName,num2str(loops),'.pdf'];
%    saveas(gcf,FileName);
%    FIVES=1
%    end
%   
%   %   pause(.01)
% end
bpmin=min(min(min(Bp(:,:,:))));
bpmax= max(max(max(Bp(:,:,:))));

BaseName='BP_GAUSIAN_TEST_P';
clf
close all

FIVES=100;
for loops=1:100:size(Bp,3)
     subplot(1,2,1)
%    contour(Bp(:,:,loops))
contour(Bp(:,:,loops),'LevelStep',.5*10^-4);
     axis([0 (z2/dz+1)/2 0 z2/dz+1 ])
 subplot(1,2,2)
 
     caxis([bpmin,bpmax]);
 
 
   surface(Bp(:,:,loops))
 axis([0 (z2/dz+1)/2 0 z2/dz+1 ])
%    axis([0 (z2/dz+1) 0 z2/dz+1 -10^-2 10^-2 ])
   title(loops)
%  axis([z1/dz+1,z2/dz+1,r1/dr+1, r2/dr+1])
%  filename=eval(sprintf('BP_GAUSIAN_TEST',loops));
   drawnow 


%        drawnow 
%    FileName=[BaseName,num2str(loops),'.pdf'];
%    saveas(gcf,FileName);
   
%       pause(1/60)
 clf
%      pause(.01)
end





% BaseName='MyDataFile_';
% for k=1:10
% FileName=[BaseName,num2str(k)]
% end



% 
% 
% 
% for loops=1:size(Bp,3)
%     
%     surf(Bp(:,:,loops))
%     axis([0 z2/dz+1 0 z2/dz+1 0e-3 10e-3])
%     title(loops)  
% % y = x.^n;
% % plot(x,y)
% % drawnow
% frame = getframe(1);
% im = frame2im(frame);
% [A,map] = rgb2ind(im,256); 
% 	if loops == 1;
% 		imwrite(A,map,'Test1','gif','LoopCount',Inf,'DelayTime',1);
% 	else
% 		imwrite(A,map,'Test1','gif','WriteMode','append','DelayTime',1);
% 	end
% end




% subplot(2,3,1)
%   surf(Vr(:,:,80))
% subplot(2,3,2)
%   surf(Vp(:,:,80))
% subplot(2,3,3)
%   surf(Vz(:,:,80))  
%   
%   
%   subplot(2,3,4)
%   surf(Br(:,:,80))
% subplot(2,3,5)
%   surf(Bp(:,:,80))
% subplot(2,3,6)
%   surf(Bz(:,:,80))  
%   
%         
% subplot(3,3,1)
% surf(Vr(:,:,2))
% title('Vr(1)')
% subplot(3,3,2)
% surf(Vr(:,:,round(size(Vr,3)/2)))
% title('Vr(mid)')
% subplot(3,3,3)
% surf(Vr(:,:,size(Vr,3)))
% title('Vr(end)')
% 
% subplot(3,3,4)
% surf(Vp(:,:,2))
% title('Vp(1)')
% subplot(3,3,5)
% surf(Vp(:,:,round(size(Vr,3)/2)))
% title('Vp(mid)')
% subplot(3,3,6)
% surf(Vp(:,:,size(Vr,3)))
% title('Vr(end)')
% 
% subplot(3,3,7)
% surf(Vz(:,:,2))
% title('Vz(1)')
% subplot(3,3,8)
% surf(Vz(:,:,round(size(Vr,3)/2)))
% title('Vz(mid)')
% subplot(3,3,9)
% surf(Vz(:,:,size(Vr,3)))
% title('Vz(end)')
% 
