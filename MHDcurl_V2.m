function [ curlR, curlPHI, curlZ ] = MHDcurl_V2( Rcomp, PHIcomp, Zcomp, r, z  )
%%Given the R,PHI, and Z components of a matrix with r and z, MHDcurl, will output the
%R PHI and Z components of Curl in cylindrical coordinates. Code was built
%under axisymmetric assumption d/dphi=0.
%
%Meathods: 4th order innards and second order edges.


R=Rcomp;
P=PHIcomp;
Z=Zcomp;

dr=r(2)-r(1);
dz=z(2)-z(1);



% Matrix Dimension Eval
col=size(R,2);
row=size(R,1);

% r Matrix ---------------------------------------------------------------
s = size(P);
rmatrix=zeros(s(1),size(r,2));
for i = 1 : s(1)
    rmatrix(i,:) = r;
end

% r curl ---------------------------------------------------------------
curlR=((1/12).*circshift(P,[-2,0])-(2/3).*circshift(P,[-1,0]) ...
      -(1/12).*circshift(P,[2,0])+(2/3).*circshift(P,[1,0]))/dz;

% z curl ---------------------------------------------------------------  
  rphi=rmatrix.*P;
zin=((1/12).*circshift(rphi,[0,-2])-(2/3).*circshift(rphi,[0,-1]) ...
      -(1/12).*circshift(rphi,[0,2])+(2/3).*circshift(rphi,[0,1]))/dr;
curlZ=1./rmatrix.*zin;
curlZ=-1*curlZ;


% phi curl ---------------------------------------------------------------
drdz=-((1/12).*circshift(R,[-2,0])-(2/3).*circshift(R,[-1,0]) ...
      -(1/12).*circshift(R,[2,0])+(2/3).*circshift(R,[1,0]))/dz;
dzdr=-((1/12).*circshift(Z,[0,-2])-(2/3).*circshift(Z,[0,-1]) ...
      -(1/12).*circshift(Z,[0,2])+(2/3).*circshift(Z,[0,1]))/dr;
curlPHI=drdz-dzdr;



% --------------------------------------------------------------
% EDGES
% dPdz(height,:)=-(-1.5*P(height,:)+2*P(height-1,:)-.5*P(height-2,:))/dz;
% dPdz(row,:) = -((-11/6)*P(row,:)+3*P(row-1,:) -1.5* P(row-2,:)+(1/3)* P(row-3,:))/dz;
% dPdz(row-1,:)=-((-11/6)*P(row-1,:)+3*P(row-2,:) -1.5* P(row-3,:)+(1/3)* P(row-4,:))/dz;
%GP
dPdz(row-1,:) = -((1/12)*P(row,:)+(-2/3)*P(row,:) +(2/3)* P(row-2,:)+(-1/12)* P(row-3,:))/dz;
dPdz(row,:) = -((1/12)*P(row,:)+(-2/3)*P(row,:) +(2/3)* P(row-1,:)+(-1/12)* P(row-2,:))/dz;

curlR(row-1,:)=-dPdz(row-1,:);
curlR(row,:)=-dPdz(row,:);
%  zinEE(:,col) =    (-1/r(col) *1*((-11/6)*rphi(:,col)+3*rphi(:,col-1) ...
%                     -1.5* rphi(:,col-2)+(1/3)* rphi(:,col-3))/dr);            
% zinEE(:,col-1) = (-1/r(col-1)) *((-11/6)*rphi(:,col-1)+3*rphi(:,col-2) ...
%                    -1.5* rphi(:,col-3)+(1/3)* rphi(:,col-4))/dr;
%GP          
zinEE(:,col-1) = -(-1/r(col-1)) *((1/12)*rphi(:,col)+(-2/3)*rphi(:,col) ...
                   +(2/3)* rphi(:,col-2)+(-1/12)* rphi(:,col-3))/dr;
zinEE(:,col) = -(-1/r(col)) *((1/12)*rphi(:,col)+(-2/3)*rphi(:,col) ...
                   +(2/3)* rphi(:,col-1)+(-1/12)* rphi(:,col-2))/dr;
               
curlZ(:,col-1)=zinEE(:,col-1);
curlZ(:,col)=zinEE(:,col);
%GP
dpdz(2,:) = ((1/12)*P(1,:)+(-2/3)*P(1,:) +(2/3)* P(3,:)-(1/12)* P(4,:))/dz;
dpdz(1,:) = ((1/12)*P(1,:)+(-2/3)*P(1,:) +(2/3)* P(2,:)-(1/12)* P(3,:))/dz;
curlR(1,:)=-dpdz(1,:);
curlR(2,:)=-dpdz(2,:);
% zinEE(:,1) = 1/r(1) * (-.5*rphi(:,3)+2*rphi(:,2)-1.5*rphi(:,1) ) / dr;
% zinEE(:,2) = 1/r(2) *  (rphi(:,3)-rphi(:,1) ) / 2 / dr;
%GP
 zinEE(:,2) = 1/r(2) *  ((1/12)*rphi(:,1)-(2/3)*rphi(:,1)+(2/3)*rphi(:,3)-(1/12)*rphi(:,4)) / dr;
 zinEE(:,1) = 1/r(1) *  ((1/12)*rphi(:,1)-(2/3)*rphi(:,1)+(2/3)*rphi(:,2)-(1/12)*rphi(:,3)) / dr;
% zinEE(:,1) =    -(-1/r(1) *1*((-11/6)*rphi(:,1)+3*rphi(:,2) ...
%                     -1.5* rphi(:,3)+(1/3)* rphi(:,4))/dr);            
% zinEE(:,2) = -(-1/r(2)) *((-11/6)*rphi(:,2)+3*rphi(:,3) ...
%                    -1.5* rphi(:,4)+(1/3)* rphi(:,5))/dr;
zinEE(:,2) = -(-1/r(2)) *((1/12)*rphi(:,1)+(-2/3)*rphi(:,1) ...
                   +(2/3)* rphi(:,3)-(1/12)* rphi(:,4))/dr;
zinEE(:,2) = -(-1/r(1)) *((1/12)*rphi(:,1)+(-2/3)*rphi(:,1) ...
                   +(2/3)* rphi(:,2)-(1/12)* rphi(:,3))/dr;               

curlZ(:,2)=zinEE(:,2);
curlZ(:,1)=zinEE(:,1);
% drdzEE(1,:)= ((-11/6)*R(1,:)+3*R(2,:) ...
%                      -1.5*R(3,:)+(1/3)*R(4,:))/dz;
% drdzEE(2,:)=((-11/6)*R(2,:)+3*R(3,:) ...
%                      -1.5*R(4,:)+(1/3)*R(5,:))/dz;
%GP
drdzEE(2,:)=((1/12)*R(1,:)+(-2/3)*R(1,:) ...
                     +(2/3)*R(3,:)+(-1/12)*R(4,:))/dz;
drdzEE(1,:)=((1/12)*R(1,:)+(-2/3)*R(1,:) ...
                     +(2/3)*R(2,:)+(-1/12)*R(3,:))/dz;
% dzdrEE(:,1)= ((-11/6)*Z(:,1)+3*Z(:,2) ...
%                      -1.5*Z(:,3)+(1/3)*Z(:,4))/dr;
% dzdrEE(:,2)=(Z(:,3)-Z(:,1))/2/dr;
%GP
dzdrEE(:,2)= ((1/12)*Z(:,1)+(-2/3)*Z(:,1) ...
                     +(2/3)*Z(:,3)+(-1/12)*Z(:,4))/dr;
dzdrEE(:,1)= ((1/12)*Z(:,1)+(-2/3)*Z(:,1) ...
                     +(2/3)*Z(:,2)+(-1/12)*Z(:,3))/dr;                 
% curlPHIEE=-drdzEE-dzdrEE;
curlPHI(2,:)=drdzEE(2,:)-dzdr(2,:);
curlPHI(:,2)=drdz(:,2)-dzdrEE(:,2);
curlPHI(:,1)=drdz(:,1)-dzdrEE(:,1);
curlPHI(1,:)=drdzEE(1,:)-dzdr(1,:);
% drdz1(row-1,:)= (R(row,:) - R(row-2,:) )/2/dz;
% drdz1(row,:)= -(-1.5*R(row,:)+2*R(row-1,:)-0.5*R(row-2,:))/dz;
%GP



%#######################################

%#######################################

%#######################################

drdz1(row-1,:)= -((1/12)*R(row,:) +(-2/3)*R(row,:)+(2/3)*R(row-2,:) +(-1/12)*R(row-3,:) )/dz;
drdz1(row,:)= -((1/12)*R(row,:) +(-2/3)*R(row,:)+(2/3)*R(row-1,:) +(-1/12)*R(row-2,:) )/dz;
% dzdr1(:,col)=(Z(:,col)-Z(:,col-1))/dr;
% dzdr1(:,col-1)=(Z(:,col)-Z(:,col-2))/2/dr; 
%GP
dzdr1(:,col-1)=-((1/12)*Z(:,col)+(-2/3)*Z(:,col)+(2/3)*Z(:,col-2)+(-1/12)*Z(:,col-3))/dr; 
dzdr1(:,col)=-((1/12)*Z(:,col)+(-2/3)*Z(:,col)+(2/3)*Z(:,col-1)+(-1/12)*Z(:,col-2))/dr; 



%#######################################

%####################################### drdz?

%####################################### dzdr?



curlPHI(row,:)=drdz1(row,:)-dzdr(row,:);
curlPHI(row-1,:)=drdz1(row-1,:)-dzdr(row-1,:);
curlPHI(:,col) =drdz(:,col) - dzdr1(:,col);
curlPHI(:,col-1) = drdz(:,col-1) - dzdr1(:,col-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (MIN, MIN)
% Curl R
% curlR(2,1) = -((-11/6)*P(2,1)+3*P(3,1) -1.5* P(4,1)+(1/3)* P(5,1))/dz;
% curlR(2,2) = -((-11/6)*P(2,2)+3*P(3,2) -1.5* P(4,2)+(1/3)* P(5,2))/dz;
% curlR(1,1) = -((-11/6)*P(1,1)+3*P(2,1) -1.5* P(3,1)+(1/3)* P(4,1))/dz;
% curlR(1,2) = -((-11/6)*P(1,2)+3*P(2,2) -1.5* P(3,2)+(1/3)* P(4,2))/dz;
%GP
curlR(2,2) = -((1/12)*P(1,2)+(-2/3)*P(1,2) +(2/3)* P(3,2)+(-1/12)* P(4,2))/dz;
curlR(2,1) = -((1/12)*P(1,1)+(-2/3)*P(1,1) +(2/3)* P(3,1)+(-1/12)* P(4,1))/dz;
curlR(1,2) = -((1/12)*P(1,2)+(-2/3)*P(1,2) +(2/3)* P(2,2)+(-1/12)* P(3,2))/dz;
curlR(1,1) = -((1/12)*P(1,1)+(-2/3)*P(1,1) +(2/3)* P(2,1)+(-1/12)* P(3,1))/dz;
% Curl Z
% dzindr(1,1) = ((-11/6)*rphi(1,1)+3*rphi(1,2) -1.5* rphi(1,3)+(1/3)* rphi(1,4))/dr;
% dzindr(2,1) = ((-11/6)*rphi(2,1)+3*rphi(2,2) -1.5* rphi(2,3)+(1/3)* rphi(2,4))/dr;
% dzindr(1,2) = ((-11/6)*rphi(1,2)+3*rphi(1,3) -1.5* rphi(1,4)+(1/3)* rphi(1,5))/dr;
% dzindr(2,2) = ((-11/6)*rphi(2,2)+3*rphi(2,3) -1.5* rphi(2,4)+(1/3)* rphi(2,5))/dr;
%GP
dzindr(2,2) = -((1/12)*rphi(2,1)+(-2/3)*rphi(2,1) +(2/3)* rphi(2,3)+(-1/12)* rphi(2,4))/dr;
dzindr(2,1) = -((1/12)*rphi(2,1)+(-2/3)*rphi(2,1) +(2/3)* rphi(2,2)+(-1/12)* rphi(2,3))/dr;
dzindr(1,2) = -((1/12)*rphi(1,1)+(-2/3)*rphi(1,1) +(2/3)* rphi(1,3)+(-1/12)* rphi(1,4))/dr;
dzindr(1,1) = -((1/12)*rphi(1,1)+(-2/3)*rphi(1,1) +(2/3)* rphi(1,2)+(-1/12)* rphi(1,3))/dr;
% dzindr(1,2) = (rphi(1,3) - rphi(1,1))/2/dr;%
% dzindr(2,2) = (rphi(2,3) - rphi(2,1))/2/dr;%
curlZ(1,1) = -1/rmatrix(1,1)*dzindr(1,1);
curlZ(2,1) = -1/rmatrix(2,1)*dzindr(2,1);
curlZ(1,2) = -1/rmatrix(1,2)*dzindr(1,2);
curlZ(2,2) = -1/rmatrix(2,2)*dzindr(2,2);
% Curl PHI
% dzdr(1,1)=(-1.5*Z(1,1)+2*Z(1,2)-0.5*Z(1,3))/dr;
% dzdr(1,1)=((-11/6)*Z(1,1)+3*Z(1,2) -1.5* Z(1,3)+(1/3)*Z(1,4))/dr;
% dzdr(2,1)=((-11/6)*Z(2,1)+3*Z(2,2) -1.5* Z(2,3)+(1/3)*Z(2,4))/dr;
% dzdr(1,2)=((-11/6)*Z(1,2)+3*Z(1,3) -1.5* Z(1,4)+(1/3)*Z(1,5))/dr;
% dzdr(2,2)=((-11/6)*Z(2,2)+3*Z(2,3) -1.5* Z(2,4)+(1/3)*Z(2,5))/dr;
%GP
dzdr(2,2)=-((1/12)*Z(2,1)+(-2/3)*Z(2,1) +(2/3)* Z(2,3)+(-1/12)*Z(2,4))/dr;
dzdr(2,1)=-((1/12)*Z(2,1)+(-2/3)*Z(2,1) +(2/3)* Z(2,2)+(-1/12)*Z(2,3))/dr;
dzdr(1,2)=-((1/12)*Z(1,1)+(-2/3)*Z(1,1) +(2/3)* Z(1,3)+(-1/12)*Z(1,4))/dr;
dzdr(1,1)=-((1/12)*Z(1,1)+(-2/3)*Z(1,1) +(2/3)* Z(1,2)+(-1/12)*Z(1,3))/dr;
% drdz(1,1)=(-1.5*R(1,1)+2*R(2,1)-.5*R(3,1))/dz;
% drdz(1,1)=((-11/6)*R(1,1)+3*R(2,1) -1.5* R(3,1)+(1/3)* R(4,1))/dz;
% drdz(1,2)=((-11/6)*R(1,2)+3*R(2,2) -1.5* R(3,2)+(1/3)* R(4,2))/dz;
% drdz(2,1)=((-11/6)*R(2,1)+3*R(3,1) -1.5* R(4,1)+(1/3)* R(5,1))/dz;
% drdz(2,2)=((-11/6)*R(2,2)+3*R(3,2) -1.5* R(4,2)+(1/3)* R(5,2))/dz;
%GP

%######################################

%######################################
drdz(2,2)=-((1/12)*R(1,2)+(-2/3)*R(1,2)+ (2/3)* R(3,2)+(-1/12)* R(4,2))/dz;
drdz(2,1)=-((1/12)*R(1,1)+(-2/3)*R(1,1)+ (2/3)* R(3,1)+(-1/12)* R(4,1))/dz;
drdz(1,2)=-((1/12)*R(1,2)+(-2/3)*R(1,2)+ (2/3)* R(2,2)+(-1/12)* R(3,2))/dz;
drdz(1,1)=-((1/12)*R(1,1)+(-2/3)*R(1,1)+ (2/3)* R(2,1)+(-1/12)* R(3,1))/dz;

%######################################

%######################################
curlPHI(1,2)=-drdz(1,2)-dzdr(1,2);
curlPHI(2,2)=-drdz(2,2)-dzdr(2,2);
curlPHI(1,1)=-drdz(1,1)-dzdr(1,1);
curlPHI(2,1)=-drdz(2,1)-dzdr(2,1);
% (MIN,MAX)
% Curl R
% curlR(1,width)=-((P(2,width) - P(1,width)))/dz;
% curlR(1,col)=-((-11/6)*P(1,col)+3*P(2,col) -1.5* P(3,col)+(1/3)* P(4,col))/dz;%
% curlR(1,col-1)=-((-11/6)*P(1,col-1)+3*P(2,col-1) -1.5* P(3,col-1)+(1/3)* P(4,col-1))/dz;
% curlR(2,col)=-((-11/6)*P(2,col)+3*P(3,col) -1.5* P(4,col)+(1/3)* P(5,col))/dz;
% curlR(2,col-1)=-((-11/6)*P(2,col-1)+3*P(3,col-1) -1.5* P(4,col-1)+(1/3)* P(5,col-1))/dz;
%GP



%#######################################

%#######################################

%#######################################

curlR(2,col-1)=-((1/12)*P(1,col-1)+(-2/3)*P(1,col-1) +(2/3)* P(3,col-1)+(-1/12)* P(4,col-1))/dz;
curlR(2,col)=-((1/12)*P(1,col)+(-2/3)*P(1,col) +(2/3)* P(3,col)+(-1/12)* P(4,col))/dz;
curlR(1,col-1)=-((1/12)*P(1,col-1)+(-2/3)*P(1,col-1) +(2/3)* P(2,col-1)+(-1/12)* P(3,col-1))/dz;
curlR(1,col)=-((1/12)*P(1,col)+(-2/3)*P(1,col) +(2/3)* P(2,col)+(-1/12)* P(3,col))/dz;



%#######################################

%#######################################

%#######################################
% Curl Z
% dzindr(1,width) = -(-0.5*rphi(1,width-2) + 2*rphi(1,width-1)-1.5*rphi(1,width))/dr;
% dzindr(1,col)=-((-11/6)*rphi(1,col)+3*rphi(1,col-1) -1.5* rphi(1,col-2)+(1/3)* rphi(1,col-3))/dr;
% dzindr(2,col) = -((-11/6)*rphi(2,col)+3*rphi(2,col-1) -1.5* rphi(2,col-2)+(1/3)* rphi(2,col-3))/dr;
% dzindr(1,col-1) =-((-11/6)*rphi(1,col-1)+3*rphi(1,col-2) -1.5* rphi(1,col-3)+(1/3)* rphi(1,col-4))/dr;
% dzindr(2,col-1) =-((-11/6)*rphi(2,col-1)+3*rphi(2,col-2) -1.5* rphi(2,col-3)+(1/3)* rphi(2,col-4))/dr;
%GP
dzindr(2,col-1) =-((1/12)*rphi(2,col)+(-2/3)*rphi(2,col) +(2/3)* rphi(2,col-2)+(-1/12)* rphi(2,col-3))/dr;
dzindr(2,col) =-((1/12)*rphi(2,col)+(-2/3)*rphi(2,col) +(2/3)* rphi(2,col-1)+(-1/12)* rphi(2,col-2))/dr;
dzindr(1,col-1) =-((1/12)*rphi(1,col)+(-2/3)*rphi(1,col) +(2/3)* rphi(1,col-2)+(-1/12)* rphi(1,col-3))/dr;
dzindr(1,col) =-((1/12)*rphi(1,col)+(-2/3)*rphi(1,col) +(2/3)* rphi(1,col-1)+(-1/12)* rphi(1,col-2))/dr;

curlZ(1,col) = -1/rmatrix(1,col)*dzindr(1,col);
curlZ(2,col) = -1/rmatrix(2,col)*dzindr(2,col);
curlZ(1,col-1)= -1/rmatrix(1,col-1)*dzindr(1,col-1);
curlZ(2,col-1) = -1/rmatrix(2,col-1)*dzindr(2,col-1);
% Curl PHI
% dzdr(1,col)=-((-11/6)*Z(1,col)+3*Z(1,col-1) -1.5* Z(1,col-2)+(1/3)*Z(1,col-3))/dr;
% dzdr(2,col)=-((-11/6)*Z(2,col)+3*Z(2,col-1) -1.5* Z(2,col-2)+(1/3)*Z(2,col-3))/dr;
% dzdr(1,col-1)=-((-11/6)*Z(1,col-1)+3*Z(1,col-2) -1.5* Z(1,col-3)+(1/3)*Z(1,col-4))/dr;
% dzdr(2,col-1)=-((-11/6)*Z(2,col-1)+3*Z(2,col-2) -1.5* Z(2,col-3)+(1/3)*Z(2,col-4))/dr;
%GP
dzdr(2,col-1)=-((1/12)*Z(2,col)+(-2/3)*Z(2,col) +(2/3)* Z(2,col-2)+(-1/12)*Z(2,col-3))/dr;
dzdr(2,col)=-((1/12)*Z(2,col)+(-2/3)*Z(2,col) +(2/3)* Z(2,col-1)+(-1/12)*Z(2,col-2))/dr;
dzdr(1,col-1)=-((1/12)*Z(1,col)+(-2/3)*Z(1,col) +(2/3)* Z(1,col-2)+(-1/12)*Z(1,col-3))/dr;
dzdr(1,col)=-((1/12)*Z(1,col)+(-2/3)*Z(1,col) +(2/3)* Z(1,col-1)+(-1/12)*Z(1,col-2))/dr;
% drdz(1,col)=((-11/6)*R(1,col)+3*R(2,col) -1.5* R(3,col)+(1/3)* R(4,col))/dz;
% drdz(1,col-1)=((-11/6)*R(1,col-1)+3*R(2,col-1) -1.5* R(3,col-1)+(1/3)* R(4,col-1))/dz;
% drdz(2,col)=((-11/6)*R(2,col)+3*R(3,col) -1.5* R(4,col)+(1/3)* R(5,col))/dz;
% drdz(2,col-1)=((-11/6)*R(2,col-1)+3*R(3,col-1) -1.5* R(4,col-1)+(1/3)* R(5,col-1))/dz;
%GP



%######################################
%######################################
drdz(2,col-1)=-((1/12)*R(1,col-1)+(-2/3)*R(1,col-1)+(2/3)* R(3,col-1)+(-1/12)* R(4,col-1))/dz;
drdz(2,col)=-((1/12)*R(1,col)+(-2/3)*R(1,col)+(2/3)* R(3,col)+(-1/12)* R(4,col))/dz;
drdz(1,col-1)=-((1/12)*R(1,col-1)+(-2/3)*R(1,col-1)+(2/3)* R(2,col-1)+(-1/12)* R(3,col-1))/dz;
drdz(1,col)=-((1/12)*R(1,col)+(-2/3)*R(1,col)+(2/3)* R(2,col)+(-1/12)* R(3,col))/dz;

%######################################
%######################################

curlPHI(1,col-1)=-drdz(1,col-1)-dzdr(1,col-1);
curlPHI(2,col-1)=-drdz(2,col-1)-dzdr(2,col-1);
curlPHI(1,col)=-drdz(1,col)-dzdr(1,col);
curlPHI(2,col)=-drdz(2,col)-dzdr(2,col);
% (MAX,MAX)
% Curl R
% % curlR(height,width)=-1*((P(height,width) - P(height-1,width)))/dz;
% curlR(row,col) = ((-11/6)*P(row,col)+3*P(row-1,col) -1.5*...
%                       P(row-2,col)+(1/3)* P(row-3,col))/dz;
% curlR(row-1,col-1)=((-11/6)*P(row-1,col-1)+3*P(row-2,col-1) -1.5*...
%                       P(row-3,col-1)+(1/3)* P(row-4,col-1))/dz;
% curlR(row,col-1) = ((-11/6)*P(row,col-1)+3*P(row-1,col-1) -1.5*...
%                       P(row-2,col-1)+(1/3)* P(row-3,col-1))/dz;
% curlR(row-1,col)=((-11/6)*P(row,col)+3*P(row-1,col) -1.5*...
%                       P(row-2,col)+(1/3)* P(row-3,col))/dz;                  
%GP
curlR(row-1,col-1)=((1/12)*P(row,col-1)+(-2/3)*P(row,col-1) +(2/3)*P(row-2,col-1)+(-1/12)* P(row-3,col-1))/dz;
curlR(row-1,col)=((1/12)*P(row,col)+(-2/3)*P(row,col) +(2/3)*P(row-2,col)+(-1/12)* P(row-3,col))/dz;
curlR(row,col-1)=((1/12)*P(row,col-1)+(-2/3)*P(row,col-1) +(2/3)*P(row-1,col-1)+(-1/12)* P(row-2,col-1))/dz;
curlR(row,col)=((1/12)*P(row,col)+(-2/3)*P(row,col) +(2/3)*P(row-1,col)+(-1/12)* P(row-2,col))/dz;
% % Curl Z
% dzindr(row,col)    =-((-11/6)*rphi(row,col)+3*rphi(row,col-1)-1.5* ...
%                      rphi(row,col-2)+(1/3)*rphi(row,col-3))/dr;
% dzindr(row-1,col)  =-((-11/6)*rphi(row-1,col)+3*rphi(row-1,col-1)-1.5* ...
%                      rphi(row-1,col-2)+(1/3)*rphi(row-1,col-3))/dr;
% dzindr(row,col-1)  =-((-11/6)*rphi(row,col-1)+3*rphi(row,col-2)-1.5* ...
%                      rphi(row,col-3)+(1/3)*rphi(row,col-4))/dr;
% dzindr(row-1,col-1)=-((-11/6)*rphi(row-1,col-1)+3*rphi(row-1,col-2)-1.5* ...
%                      rphi(row-1,col-3)+(1/3)*rphi(row-1,col-4))/dr;
%GP                 
dzindr(row-1,col-1)=-((1/12)*rphi(row-1,col)+(-2/3)*rphi(row-1,col)+(2/3)*rphi(row-1,col-2)+(-1/12)*rphi(row-1,col-3))/dr;
dzindr(row-1,col)=-((1/12)*rphi(row-1,col)+(-2/3)*rphi(row-1,col)+(2/3)*rphi(row-1,col-1)+(-1/12)*rphi(row-1,col-2))/dr;
dzindr(row,col-1)=-((1/12)*rphi(row,col)+(-2/3)*rphi(row,col)+(2/3)*rphi(row,col-2)+(-1/12)*rphi(row,col-3))/dr;
dzindr(row,col)=-((1/12)*rphi(row,col)+(-2/3)*rphi(row,col)+(2/3)*rphi(row,col-1)+(-1/12)*rphi(row,col-2))/dr;

curlZ(row,col)= 1./rmatrix(row,col).*dzindr(row,col);
curlZ(row,col-1)= 1./rmatrix(row,col-1).*dzindr(row,col-1);
curlZ(row-1,col)= 1./rmatrix(row-1,col).*dzindr(row-1,col);
curlZ(row-1,col-1)= 1./rmatrix(row-1,col-1).*dzindr(row-1,col-1);
% % Curl PHI
% dzdr(row,col)=-((-11/6)*Z(row,col)+3*Z(row,col-1) -1.5* Z(row,col-2)+(1/3)*Z(row,col-3))/dr;
% dzdr(row-1,col)=-((-11/6)*Z(row-1,col)+3*Z(row-1,col-1) -1.5* Z(row-1,col-2)+(1/3)*Z(row-1,col-3))/dr;
% dzdr(row,col-1)=-((-11/6)*Z(row,col-1)+3*Z(row,col-2) -1.5* Z(row,col-3)+(1/3)*Z(row,col-4))/dr;
% dzdr(row-1,col-1)=-((-11/6)*Z(row-1,col-1)+3*Z(row-1,col-2) -1.5* Z(row-1,col-3)+(1/3)*Z(row-1,col-4))/dr;
%GP
dzdr(row-1,col-1)=-((1/12)*Z(row-1,col)+(-2/3)*Z(row-1,col)+(2/3)* Z(row-1,col-2)+(-1/12)*Z(row-1,col-3))/dr;
dzdr(row-1,col)=-((1/12)*Z(row-1,col)+(-2/3)*Z(row-1,col)+(2/3)* Z(row-1,col-1)+(-1/12)*Z(row-1,col-2))/dr;
dzdr(row,col-1)=-((1/12)*Z(row,col)+(-2/3)*Z(row,col)+(2/3)* Z(row,col-2)+(-1/12)*Z(row,col-3))/dr;
dzdr(row,col)=-((1/12)*Z(row,col)+(-2/3)*Z(row,col)+(2/3)* Z(row,col-1)+(-1/12)*Z(row,col-2))/dr; 
% drdz(row,col)=  -((-11/6)*R(row,col)+3*  R(row-1,col) -1.5*   R(row-2,col)+(1/3)*   R(row-3,col))/dz;
% drdz(row,col-1)=-((-11/6)*R(row,col-1)+3*R(row-1,col-1) -1.5* R(row-2,col-1)+(1/3)* R(row-3,col-1))/dz;
% drdz(row-1,col)=  -((-11/6)*R(row-1,col)+3*  R(row-2,col) -1.5*   R(row-3,col)+(1/3)*   R(row-4,col))/dz;
% drdz(row-1,col-1)=-((-11/6)*R(row-1,col-1)+3*R(row-2,col-1) -1.5* R(row-3,col-1)+(1/3)* R(row-4,col-1))/dz;
%GP
drdz(row-1,col-1)=-((1/12)*R(row,col-1)+(-2/3)*R(row,col-1)+(2/3)* R(row-2,col-1)+(-1/12)* R(row-3,col-1))/dz;
drdz(row-1,col)=-((1/12)*R(row,col)+(-2/3)*R(row,col)+(2/3)* R(row-2,col)+(-1/12)* R(row-3,col))/dz;
drdz(row,col-1)=-((1/12)*R(row,col-1)+(-2/3)*R(row,col-1)+(2/3)* R(row-1,col-1)+(-1/12)* R(row-2,col-1))/dz;
drdz(row,col)=-((1/12)*R(row,col)+(-2/3)*R(row,col)+(2/3)* R(row-1,col)+(-1/12)* R(row-2,col))/dz;

curlPHI(row-1,col)=drdz(row-1,col)-dzdr(row-1,col);
curlPHI(row-1,col-1)=drdz(row-1,col-1)-dzdr(row-1,col-1);
curlPHI(row,col)=drdz(row,col)-dzdr(row,col);
curlPHI(row,col-1)=drdz(row,col-1)-dzdr(row,col-1);
% % (MAX,MIN)
% curlR(row,1) = ((-11/6)*P(row,1)+3*P(row-1,1) -1.5*...
%                       P(row-2,1)+(1/3)* P(row-3,1))/dz;
% curlR(row-1,2)=((-11/6)*P(row-1,2)+3*P(row-2,2) -1.5*...
%                       P(row-3,2)+(1/3)* P(row-4,2))/dz;
% curlR(row,2) = ((-11/6)*P(row,2)+3*P(row-1,2) -1.5*...
%                       P(row-2,2)+(1/3)* P(row-3,2))/dz;
% curlR(row-1,1)=((-11/6)*P(row,1)+3*P(row-1,1) -1.5*...
%                       P(row-2,1)+(1/3)* P(row-3,1))/dz; 
%GP


%######################################
%######################################
%######################################


curlR(row-1,2) = ((1/12)*P(row,2)+(-2/3)*P(row,2) +(2/3)*P(row-2,2)+(-1/12)* P(row-3,2))/dz;
curlR(row-1,1) = ((1/12)*P(row,1)+(-2/3)*P(row,1) +(2/3)*P(row-2,1)+(-1/12)* P(row-3,1))/dz;
curlR(row,2) = ((1/12)*P(row,2)+(-2/3)*P(row,2) +(2/3)*P(row-1,2)+(-1/12)* P(row-2,2))/dz;
curlR(row,1) = ((1/12)*P(row,1)+(-2/3)*P(row,1) +(2/3)*P(row-1,1)+(-1/12)* P(row-2,1))/dz;



% Curl Z
% dzindr(row,1)    =((-11/6)*rphi(row,1)+3*rphi(row,2)-1.5* ...
%                      rphi(row,3)+(1/3)*rphi(row,4))/dr;
% dzindr(row-1,1)  =((-11/6)*rphi(row-1,1)+3*rphi(row-1,2)-1.5* ...
%                      rphi(row-1,3)+(1/3)*rphi(row-1,4))/dr;
% dzindr(row,2)  =((-11/6)*rphi(row,2)+3*rphi(row,3)-1.5* ...
%                      rphi(row,4)+(1/3)*rphi(row,5))/dr;
% dzindr(row-1,2)=((-11/6)*rphi(row-1,2)+3*rphi(row-1,3)-1.5* ...
%                      rphi(row-1,4)+(1/3)*rphi(row-1,5))/dr;
%GP                 
dzindr(row-1,2)  =((1/12)*rphi(row-1,1)+(-2/3)*rphi(row-1,1)+(2/3)*rphi(row-1,3)+(-1/12)*rphi(row-1,4))/dr;
dzindr(row-1,1)  =((1/12)*rphi(row-1,1)+(-2/3)*rphi(row-1,1)+(2/3)*rphi(row-1,2)+(-1/12)*rphi(row-1,3))/dr;
dzindr(row,2)  =((1/12)*rphi(row,1)+(-2/3)*rphi(row,1)+(2/3)*rphi(row,3)+(-1/12)*rphi(row,4))/dr;
dzindr(row,1)  =((1/12)*rphi(row,1)+(-2/3)*rphi(row,1)+(2/3)*rphi(row,2)+(-1/12)*rphi(row,3))/dr;
                
curlZ(row,1)= 1./rmatrix(row,1).*dzindr(row,1);
curlZ(row,2)= 1./rmatrix(row,2).*dzindr(row,2);
curlZ(row-1,1)= 1./rmatrix(row-1,1).*dzindr(row-1,1);
curlZ(row-1,2)= 1./rmatrix(row-1,2).*dzindr(row-1,2);
% Curl PHI
% dzdr(row,1)=((-11/6)*Z(row,1)+3*Z(row,2) -1.5* Z(row,3)+(1/3)*Z(row,4))/dr;
% dzdr(row-1,1)=((-11/6)*Z(row-1,1)+3*Z(row-1,2) -1.5* Z(row-1,3)+(1/3)*Z(row-1,4))/dr;
% dzdr(row,2)=((-11/6)*Z(row,2)+3*Z(row,3) -1.5* Z(row,4)+(1/3)*Z(row,5))/dr;
% dzdr(row-1,2)=((-11/6)*Z(row-1,2)+3*Z(row-1,3) -1.5* Z(row-1,4)+(1/3)*Z(row-1,5))/dr;
%GP
dzdr(row-1,2)=((1/12)*Z(row-1,1)+(-2/3)*Z(row-1,1)+(2/3)*Z(row-1,3)+(-1/12)*Z(row-1,4))/dr;
dzdr(row-1,1)=((1/12)*Z(row-1,1)+(-2/3)*Z(row-1,1)+(2/3)*Z(row-1,2)+(-1/12)*Z(row-1,3))/dr;
dzdr(row,2)=((1/12)*Z(row,1)+(-2/3)*Z(row,1)+(2/3)*Z(row,3)+(-1/12)*Z(row,4))/dr;
dzdr(row,1)=((1/12)*Z(row,1)+(-2/3)*Z(row,1)+(2/3)*Z(row,2)+(-1/12)*Z(row,3))/dr;
% drdz(row,1)=  -((-11/6)*R(row,1)+3*  R(row-1,1) -1.5*   R(row-2,1)+(1/3)*   R(row-3,1))/dz;
% drdz(row,2)=-((-11/6)*R(row,2)+3*R(row-1,2) -1.5* R(row-2,2)+(1/3)* R(row-3,2))/dz;
% drdz(row-1,1)=  -((-11/6)*R(row-1,1)+3*  R(row-2,1) -1.5*   R(row-3,1)+(1/3)*   R(row-4,1))/dz;
% drdz(row-1,2)=-((-11/6)*R(row-1,2)+3*R(row-2,2) -1.5* R(row-3,2)+(1/3)* R(row-4,2))/dz;
%GP
drdz(row-1,2)=-((1/12)*R(row,2)+(-2/3)*R(row,2)+(2/3)* R(row-2,2)+(-1/12)* R(row-3,2))/dz;
drdz(row-1,1)=-((1/12)*R(row,1)+(-2/3)*R(row,1)+(2/3)* R(row-2,1)+(-1/12)* R(row-3,1))/dz;
drdz(row,2)=-((1/12)*R(row,2)+(-2/3)*R(row,2)+(2/3)* R(row-1,2)+(-1/12)* R(row-2,2))/dz;
drdz(row,1)=-((1/12)*R(row,1)+(-2/3)*R(row,1)+(2/3)* R(row-1,1)+(-1/12)* R(row-2,1))/dz;


curlPHI(row-1,1)=drdz(row-1,1)-dzdr(row-1,1);
curlPHI(row-1,2)=drdz(row-1,2)-dzdr(row-1,2);
curlPHI(row,1)=drdz(row,1)-dzdr(row,1);
curlPHI(row,2)=drdz(row,2)-dzdr(row,2);

end




