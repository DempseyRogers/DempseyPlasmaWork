function [ divRZ,divRZE ] = MHDdivergence(Rcomp, PHIcomp, Zcomp, r, z  )
%MHDgradient Takes a vector, r phi z, to a multivaried function. This
%function opporates under the axisymmetric assumption, dphi/dt=0.
%4th order
R=Rcomp;
P=PHIcomp;
Z=Zcomp;

dr=r(2)-r(1);
dz=z(2)-z(1);

% r Matrix ---------------------------------------------------------------
s = size(R);
rmatrix=zeros(s(1),size(r,2));
for i = 1 : s(1)
    rmatrix(i,:) = r;
end

%INTERIOR 
% r Divergance ---------------------------------------------------------------
vrmat=rmatrix.*R;
dfdr=((1/12).*circshift(vrmat,[0,-2])-(2/3).*circshift(vrmat,[0,-1]) ...
      -(1/12).*circshift(vrmat,[0,2])+(2/3).*circshift(vrmat,[0,1]))/dr;
Rcomponent=1./rmatrix.*dfdr;  
% z Divergance ---------------------------------------------------------------

Zcomponent=((1/12).*circshift(Z,[-2,0])-(2/3).*circshift(Z,[-1,0]) ...
      -(1/12).*circshift(Z,[2,0])+(2/3).*circshift(Z,[1,0]))/dz;

  Zcomponent=-1*Zcomponent;
  
divRZ=Rcomponent+Zcomponent; 










% Edges
dfdrE=((-1).*circshift(vrmat,[0,1])+circshift(vrmat,[0,0]))/dr;
RcomponentE=1./rmatrix.*dfdrE;  
% z Divergance ---------------------------------------------------------------

ZcomponentE=((-1).*circshift(Z,[1,0])+circshift(Z,[0,0]))/dz;

%  ZcomponentE=-1*ZcomponentE;
  
divRZE=RcomponentE+ZcomponentE; 

divRZ(size(divRZ,1)-1,:)=divRZE(size(divRZE,1)-1,:);
divRZ(size(divRZ,1),:)=divRZE(size(divRZE,1),:);


divRZ(:,size(divRZ,2)-1)=divRZE(:,size(divRZE,2)-1);
divRZ(:,size(divRZ,2))=divRZE(:,size(divRZE,2));




dfdrEE=((-1).*circshift(vrmat,[0,0])+circshift(vrmat,[0,-1]))/dr;
RcomponentEE=1./rmatrix.*dfdrEE;  
% z Divergance ---------------------------------------------------------------

ZcomponentEE=((-1).*circshift(Z,[-1,0])+circshift(Z,[0,0]))/dz;

ZcomponentEE=-1*ZcomponentEE;
  
divRZEE=RcomponentEE+ZcomponentEE; 

divRZ(1,:)=divRZEE(1,:);
divRZ(2,:)=divRZEE(2,:);


divRZ(:,1)=divRZEE(:,1);
divRZ(:,2)=divRZEE(:,2);

% fix up corner
divRZ(size(divRZ,1),1)=(Z(size(divRZ,1),1) - Z(size(divRZ,1)-1,1))/dz+...
    (Z(size(divRZ,1),2) - Z(size(divRZ,1),1))/dr;

divRZ(size(divRZ,1),2)=(Z(size(divRZ,1),2) - Z(size(divRZ,1)-1,2))/dz+...
    (Z(size(divRZ,1),3) - Z(size(divRZ,1),1))/2/dr;

end
% R=Rcomp;
% P=PHIcomp;
% Z=Zcomp;
% 
% dr=r(2)-r(1);
% dz=z(2)-z(1);
% 
% % r Matrix ---------------------------------------------------------------
% s = size(R);
% rmatrix=zeros(s(1),size(r,2));
% for i = 1 : s(1)
%     rmatrix(i,:) = r;
% end
% 
% 
% % r Divergance ---------------------------------------------------------------
% vrmat=rmatrix.*R;
% dfdr=(-(1/2).*circshift(vrmat,[0,-1])+(1/2).*circshift(vrmat,[0,1]))/dr;
% 
% 
% % z Divergance ---------------------------------------------------------------
% Rcomponent=1./r.*dfdr;
% Zcomponent=(-(1/2).*circshift(Z,[0,-1])+(1/2).*circshift(Z,[0,1]))/dz;
% 
% divRZ=Rcomponent+Zcomponent; 
% 
% end

