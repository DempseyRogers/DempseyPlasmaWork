function [ ddx ] = MHDgradR( rcomp ,r )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

R=rcomp;

dr=r(2)-r(1);

rows=size(rcomp,1);
cols=size(rcomp,2);

ddx= ((1/12)*circshift(R,[0,-2])+ ...
     (-2/3)*circshift(R,[0,-1])+ ...
      (2/3)*circshift(R,[0,1])+ ...
    (-1/12)*circshift(R,[0,2]))/dr;





ddx(2,:)=((1/12)*R(2,:)+ ...
     (-2/3)*R(2,:)+ ...
      (2/3)*R(3,:)+ ...
    (-1/12)*R(4,:))/dr;

ddx(1,:)=((1/12)*R(1,:)+ ...
     (-2/3)*R(1,:)+ ...
      (2/3)*R(2,:)+ ...
    (-1/12)*R(3,:))/dr;




ddx(rows-1,:)=((1/12)*R(rows-1,:)+ ...
     (-2/3)*R(rows-1,:)+ ...
      (2/3)*R(rows-2,:)+ ...
    (-1/12)*R(rows-3,:))/dr;

ddx(rows,:)=((1/12)*R(rows,:)+ ...
     (-2/3)*R(rows,:)+ ...
      (2/3)*R(rows-1,:)+ ...
    (-1/12)*R(rows-2,:))/dr;






end

