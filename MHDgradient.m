function [ rcomp, zcomp ] = MHDgradient(fxn, r, z )
%MHDgradient Takes a multivaried function, r phi z, to vector. This
%function opporates under the axisymmetric assumption, dphi/dt=0.

dr=r(2)-r(1);
dz=z(2)-z(1);

rcomp=((1/12).*circshift(fxn,[0,-2])-(2/3).*circshift(fxn,[0,-1]) ...
      -(1/12).*circshift(fxn,[0,2])+(2/3).*circshift(fxn,[0,1]))/dr;

zcomp=((1/12).*circshift(fxn,[-2,0])-(2/3).*circshift(fxn,[-1,0]) ...
      -(1/12).*circshift(fxn,[2,0])+(2/3).*circshift(fxn,[1,0]))/dz;
  zcomp=-1*zcomp;



end