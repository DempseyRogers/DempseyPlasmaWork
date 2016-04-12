function [ ddz ] = MHDgradZ( zcomp ,z )
%Given the Z component of a function and z, MHDgradZ, will output the
%derivative of zcomp with respect to z,(d(zcomp)/dz).
% 
%Meathods: 4th order inneards and second order edges.
   

Z=zcomp;
% Width=size(Z,2);
row=size(Z,1);

dz=z(2)-z(1);


ddz=-((1/12).*circshift(Z,[-2,0])-(2/3).*circshift(Z,[-1,0]) ...
     -(1/12).*circshift(Z,[2, 0])+(2/3).*circshift(Z,[1, 0]))/dz;

  
  ddz(1,:)=((-11/6)*Z(1,:)+3*Z(2,:)-1.5*Z(3,:)+(1/3)*Z(4,:))/dz;
  ddz(2,:)=((-11/6)*Z(2,:)+3*Z(3,:)-1.5*Z(4,:)+(1/3)*Z(5,:))/dz;
  ddz(row,:)=-((-11/6)*Z(row,:)+3*Z(row-1,:)-1.5*Z(row-2,:)+(1/3)*Z(row-3,:))/dz;
  ddz(row-1,:)=-((-11/6)*Z(row-1,:)+3*Z(row-2,:)-1.5*Z(row-3,:)+(1/3)*Z(row-4,:))/dz;
  %GP
%     ddz(2,:)=((1/12)*Z(2,:)+(-2/3)*Z(2,:)+(2/3)*Z(3,:)+(-1/12)*Z(4,:))/dz;
%     ddz(1,:)=((1/12)*Z(1,:)+(-2/3)*Z(1,:)+(2/3)*Z(2,:)+(-1/12)*Z(3,:))/dz;
%     ddz(row-1,:)=-((1/12)*Z(row-1,:)+(-2/3)*Z(row-1,:)+(2/3)*Z(row-2,:)+(-1/12)*Z(row-3,:))/dz;
%     ddz(row,:)=-((1/12)*Z(row,:)+(-2/3)*Z(row,:)+(2/3)*Z(row-1,:)+(-1/12)*Z(row-2,:))/dz;
%   
    
    
    
end

