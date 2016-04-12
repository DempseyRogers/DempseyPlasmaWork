function [ ddr ] = MHDgradR_V2( rcomp ,r )
%Given the R component of a function and r, MHDgradR, will output the
%derivative of zcomp with respect to r,(d(rcomp)/dr). 
% 
%Meathods: 4th order inneards and second order edges.

R=rcomp;
% row=size(R,2);
col=size(R,2);
% Height=size(R,1);

dr=r(2)-r(1);


ddr=-((1/12).*circshift(R,[0,-2])-(2/3).*circshift(R,[0,-1]) ...
     -(1/12).*circshift(R,[0, 2])+(2/3).*circshift(R,[0, 1]))/dr;

  
%   ddr(:,1)= ((-11/6)*R(:,1)+3*R(:,1) -1.5* R(:,3)+(1/3)* R(:,4))/dr;
%   ddr(:,2)=((-11/6)*R(:,2)+3*R(:,3) -1.5* R(:,4)+(1/3)* R(:,5))/dr;
%   ddr(:,col)=-((-11/6)*R(:,col)+3*R(:,col-1) -1.5* R(:,col-2)+(1/3)* R(:,col-3))/dr;
%   ddr(:,col-1)=-((-11/6)*R(:,col-1)+3*R(:,col-2) -1.5* R(:,col-3)+(1/3)* R(:,col-4))/dr;
   %GP  
%     ddr(:,2)=((1/12)*R(:,2)+(-2/3)*R(:,2)+(2/3)* R(:,3)+(-1/12)* R(:,4))/dr;
%     ddr(:,1)=((1/12)*R(:,1)+(-2/3)*R(:,1)+(2/3)* R(:,2)+(-1/12)* R(:,3))/dr;
%     ddr(:,col-1)=-((1/12)*R(:,col-1)+(-2/3)*R(:,col-1)+(2/3)* R(:,col-2)+(1/12)* R(:,col-3))/dr;
%     ddr(:,col)=-((1/12)*R(:,col)+(-2/3)*R(:,col)+(2/3)* R(:,col-1)+(1/12)* R(:,col-2))/dr;
  
    ddr(:,2)=((1/12)*R(:,1)+(-2/3)*R(:,1)+(2/3)* R(:,3)+(-1/12)* R(:,4))/dr;
    ddr(:,1)=((1/12)*R(:,1)+(-2/3)*R(:,1)+(2/3)* R(:,2)+(-1/12)* R(:,3))/dr;
    
    ddr(:,col-1)=-((1/12)*R(:,col)+(-2/3)*R(:,col)+(2/3)* R(:,col-2)+(-1/12)* R(:,col-3))/dr;
    ddr(:,col)=-((1/12)*R(:,col)+(-2/3)*R(:,col)+(2/3)* R(:,col-1)+(-1/12)* R(:,col-2))/dr;  

 
    
    
end

