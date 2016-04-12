function [ laplacianRZ ] = MHDlaplacian( fxn, r, z )
%%Given a multivarried function with r and z, MHDlaplacian, will output the
%the functions laplacian opporator in cylindrical coordinates. Code was built
%under axisymmetric assumption d/dphi=0.
% 
%Meathods: 4th order inneards and second order edges.



%BEGINING OF MHDlaplacian

dr=r(2)-r(1);
dz = z(2)-z(1);

% Matrix dimension Eval
col=size(fxn,2);
row=size(fxn,1);

    
    rmatrix=zeros(row, col);
    for i = 1 : row
        rmatrix(i,:) = r;
    end
% r Laplacian ---------------------------------------------------------------
dfdr=-((1/12).*circshift(fxn,[0,-2])-(2/3).*circshift(fxn,[0,-1]) ...
      -(1/12).*circshift(fxn,[0,2])+(2/3).*circshift(fxn,[0,1]))/dr;
  

d2fdr2=(-(1/12).*circshift(fxn,[0,-2])+(4/3).*circshift(fxn,[0,-1])...
        -(5/2).*circshift(fxn,[0,0]) ...
        -(1/12).*circshift(fxn,[0,2]) +(4/3).*circshift(fxn,[0,1]))/dr^2;
  
Rcomponent=1./rmatrix.*(dfdr+rmatrix.*d2fdr2);

% z Laplacian --------------------------------------------------------------
Zcomponent=(-(1/12).*circshift(fxn,[-2,0])+(4/3).*circshift(fxn,[-1,0])...
           -(5/2).*circshift(fxn,[0,0]) ...
            -(1/12).*circshift(fxn,[2,0]) +(4/3).*circshift(fxn,[1,0]))/dz^2;

% Zcomponent = (circshift(fxn,[-1,0])-2*fxn+circshift(fxn,[1,0]))/dz^2;       
























% Edge Work

% r Laplacian ---------------------------------------------------------------
% dfdrE(:,1)=((-3/2)*fxn(:,1)+ ...
%                (2)*fxn(:,2)+ ...
%             (-1/2)*fxn(:,3))/dr;        
% dfdrE(:,2)=((-1/2)*fxn(:,1)+ ...
%              ( 1/2)*fxn(:,3))/dr;
% 
% d2fdr2E(:,1)=((2)*fxn(:,1)+ ...
%              (-5)*fxn(:,2)+ ...  
%               (4)*fxn(:,3)+ ...
%              (-1)*fxn(:,4))/dr^2;
% d2fdr2E(:,2)=(fxn(:,1)+ ...
%              (-2)*fxn(:,2)+ ...  
%              fxn(:,3))/dz^2; 
         %GP         
dfdrE(:,2)=((1/12)*fxn(:,1)+( -2/3)*fxn(:,1)+(2/3)*fxn(:,3)+( -1/12)*fxn(:,4))/dr;
dfdrE(:,1)=((1/12)*fxn(:,1)+( -2/3)*fxn(:,1)+(2/3)*fxn(:,2)+( -1/12)*fxn(:,3))/dr;
         
d2fdr2E(:,2)=((-1/12)*fxn(:,1)+(4/3)*fxn(:,1)+(-5/2)*fxn(:,2)+(4/3)*fxn(:,3)+(-1/12)*fxn(:,4))/dz^2; 
d2fdr2E(:,1)=((-1/12)*fxn(:,1)+(4/3)*fxn(:,1)+(-5/2)*fxn(:,1)+(4/3)*fxn(:,2)+(-1/12)*fxn(:,3))/dz^2; 
         
         
         
         
         
         
% dfdrE(:,col)=-((-3/2)*fxn(:,col)+ ...
%                    (2)*fxn(:,col-1)+ ...
%                 (-1/2)*fxn(:,col-2))/dr;              
% dfdrE(:,col-1)=-((-1/2)*fxn(:,col)+ ...
%                    ( 1/2)*fxn(:,col-2))/dr;
%          
% d2fdr2E(:,col)=((2)*fxn(:,col)+ ...
%                   (-5)*fxn(:,col-1)+ ...  
%                    (4)*fxn(:,col-2)+ ...
%                   (-1)*fxn(:,col-3))/dr^2;              
% d2fdr2E(:,col-1)=(fxn(:,col)+ ...
%                    (-2)*fxn(:,col-1)+ ...  
%                    fxn(:,col-2))/dz^2; 
%GP               
dfdrE(:,col-1)=((1/12)*fxn(:,col)+( -2/3)*fxn(:,col)+(2/3)*fxn(:,col-2)+( -1/12)*fxn(:,col-3))/dr;
dfdrE(:,col)=((1/12)*fxn(:,col)+( -2/3)*fxn(:,col)+(2/3)*fxn(:,col-1)+( -1/12)*fxn(:,col-2))/dr;  

d2fdr2E(:,col-1)=((-1/12)*fxn(:,col)+(4/3)*fxn(:,col)+(-5/2)*fxn(:,col-1)+(4/3)*fxn(:,col-2)+(-1/12)*fxn(:,col-3))/dz^2; 
d2fdr2E(:,col)=((-1/12)*fxn(:,col)+(4/3)*fxn(:,col)+(-5/2)*fxn(:,col)+(4/3)*fxn(:,col-1)+(-1/12)*fxn(:,col-2))/dz^2; 
               
               
% Rcomponent(:,col)=1./rmatrix(:,col).*(dfdrE(:,col)+rmatrix(:,col).*d2fdr2E(:,col));
% Rcomponent(:,col-1)=1./rmatrix(:,col-1).*(dfdrE(:,col-1)+rmatrix(:,col-1).*d2fdr2E(:,col-1));
% Rcomponent(:,1)=1./rmatrix(:,1).*(dfdrE(:,1)+rmatrix(:,1).*d2fdr2E(:,1));
% Rcomponent(:,2)=1./rmatrix(:,2).*(dfdrE(:,2)+rmatrix(:,2).*d2fdr2E(:,2));
              
               
Rcomponent(:,col)=1./rmatrix(:,col).*(dfdrE(:,col)+rmatrix(:,col).*d2fdr2E(:,col));
Rcomponent(:,col-1)=1./rmatrix(:,col-1).*(dfdrE(:,col-1)+rmatrix(:,col-1).*d2fdr2E(:,col-1));
Rcomponent(:,1)=1./rmatrix(:,1).*(dfdrE(:,1)+rmatrix(:,1).*d2fdr2E(:,1));
Rcomponent(:,2)=1./rmatrix(:,2).*(dfdrE(:,2)+rmatrix(:,2).*d2fdr2E(:,2));

% Edge Work 

% z Laplacaian
% First Row 

       

% % Last Row                  
% Zcomponent(row,:)=((2)*fxn(row,:)+ ...
%                   (-5)*fxn(row-1,:)+ ...  
%                    (4)*fxn(row-2,:)+ ...
%                   (-1)*fxn(row-3,:))/dz^2;
%               
% Zcomponent(row-1,:)=(fxn(row,:)+ ...
%                        (-2)*fxn(row-1,:)+ ...  
%                        fxn(row-2,:))/dz^2; 
%GP
Zcomponent(2,:)=((-1/12)*fxn(1,:)+ (4/3)*fxn(1,:)+(-5/2)*fxn(2,:)+ (4/3)*fxn(3,:)+ (-1/12)*fxn(4,:))/dz^2;
Zcomponent(1,:)=((-1/12)*fxn(1,:)+ (4/3)*fxn(1,:)+(-5/2)*fxn(1,:)+ (4/3)*fxn(2,:)+ (-1/12)*fxn(3,:))/dz^2;
Zcomponent(row-1,:)=((-1/12)*fxn(row,:)+ (4/3)*fxn(row,:)+(-5/2)*fxn(row-1,:)+ (4/3)*fxn(row-2,:)+ (-1/12)*fxn(row-3,:))/dz^2;
Zcomponent(row,:)=((-1/12)*fxn(row,:)+ (4/3)*fxn(row,:)+(-5/2)*fxn(row,:)+ (4/3)*fxn(row-1,:)+ (-1/12)*fxn(row-2,:))/dz^2;


% surf(laplacianRZ)
% (row,col)
% (min,min) 

% r Laplacian ---------------------------------------------------------------
% dfdr(2,1)=((-3/2)*fxn(2,3)+2*fxn(2,2)-(1/2)*fxn(2,1))/dr; 
% dfdr(1,2)=(-.5*fxn(1,3)+.5*fxn(1,1))/dr;
% dfdr(2,2)=(-.5*fxn(2,3)+.5*fxn(2,1))/dr;
%GP
dfdr(2,2)=((1/12)*fxn(2,1)+(-2/3)*fxn(2,1)+(2/3)*fxn(2,3)+(-1/12)*fxn(2,4))/dr;
dfdr(2,1)=((1/12)*fxn(2,1)+(-2/3)*fxn(2,1)+(2/3)*fxn(2,2)+(-1/12)*fxn(2,3))/dr;
dfdr(1,2)=((1/12)*fxn(1,1)+(-2/3)*fxn(1,1)+(2/3)*fxn(1,3)+(-1/12)*fxn(1,4))/dr;
dfdr(1,1)=((1/12)*fxn(1,1)+(-2/3)*fxn(1,1)+(2/3)*fxn(1,2)+(-1/12)*fxn(1,3))/dr;


% d2fdr2(1,1)=-((2)*fxn(1,4)-5*fxn(1,3)+(4)*fxn(1,2)-fxn(1,1))/(dr^2);    
% d2fdr2(2,1)=-((2)*fxn(2,4)-5*fxn(2,3)+(4)*fxn(2,2)-fxn(2,1))/(dr^2); 
% d2fdr2(1,2) = (fxn(1,1)-2*fxn(1,2) ...
%                            +fxn(1,3) )/dz^2;
% d2fdr2(2,2) = (fxn(2,1)-2*fxn(2,2) ...
%                            +fxn(2,3) )/dz^2;
%GP
d2fdr2(2,2)=-((-1/12)*fxn(2,1)-(4/3)*fxn(2,1)+(-5/2)*fxn(2,2)+(4/3)*fxn(2,3)-(-1/12)*fxn(2,4))/(dr^2);  
d2fdr2(2,1)=-((-1/12)*fxn(2,1)-(4/3)*fxn(2,1)+(-5/2)*fxn(2,1)+(4/3)*fxn(2,2)-(-1/12)*fxn(2,3))/(dr^2);  
d2fdr2(1,2)=-((-1/12)*fxn(1,1)-(4/3)*fxn(1,1)+(-5/2)*fxn(1,2)+(4/3)*fxn(1,3)-(-1/12)*fxn(1,4))/(dr^2);  
d2fdr2(1,1)=-((-1/12)*fxn(1,1)-(4/3)*fxn(1,1)+(-5/2)*fxn(1,1)+(4/3)*fxn(1,2)-(-1/12)*fxn(1,3))/(dr^2);  




% (min, max)

% dfdr(1,col)=((-3/2)*fxn(1,col)+2*fxn(1,col-1)-(1/2)*fxn(1,col-2))/dr;
% dfdr(2,col)=((-3/2)*fxn(2,col)+2*fxn(2,col-1)-(1/2)*fxn(2,col-2))/dr; 
% dfdr(1,col-1)=(.5*fxn(1,col)-.5*fxn(1,col-2))/dr;
% dfdr(2,col-1)=(.5*fxn(2,col)-.5*fxn(2,col-2))/dr;
%  d2fdr2(1,col)=-((2)*fxn(1,col-3)-5*fxn(1,col-2)+(4)*fxn(1,col-1)-fxn(1,col))/(dr^2);    
%  d2fdr2(2,col)=-((2)*fxn(2,col-3)-5*fxn(2,col-2)+(4)*fxn(2,col-1)-fxn(2,col))/(dr^2);         
%GP
dfdr(2,col-1)=((1/12)*fxn(2,col)+(-2/3)*fxn(2,col)+(2/3)*fxn(2,col-2)+(-1/12)*fxn(2,col-3))/dr;
dfdr(2,col)=((1/12)*fxn(2,col)+(-2/3)*fxn(2,col)+(2/3)*fxn(2,col-1)+(-1/12)*fxn(2,col-2))/dr;
dfdr(1,col-1)=((1/12)*fxn(1,col)+(-2/3)*fxn(1,col)+(2/3)*fxn(1,col-2)+(-1/12)*fxn(1,col-3))/dr;
dfdr(1,col)=((1/12)*fxn(1,col)+(-2/3)*fxn(1,col)+(2/3)*fxn(1,col-1)+(-1/12)*fxn(1,col-2))/dr;



% d2fdr2(1,col)=-(-1*fxn(1,col-3)+4*fxn(1,col-2)-5*fxn(1,col-1)+2*fxn(1,col))/(dr^2);    
% d2fdr2(2,col)=-(-1*fxn(2,col-3)+4*fxn(2,col-2)-5*fxn(2,col-1)+2*fxn(2,col))/(dr^2);         
% d2fdr2(1,col-1) = (fxn(1,col) - 2*fxn(1,col-1) ...
%                            +fxn(1,col-2) )/dz^2;
% d2fdr2(2,col-1) = (fxn(2,col) - 2*fxn(2,col-1) ...
%                            +fxn(2,col-2) )/dz^2;
%GP     
d2fdr2(2,col-1)=((-1/12)*fxn(2,col)+(4/3)*fxn(2,col)+(-5/2)*fxn(2,col-1)+(4/3)*fxn(2,col-2)+(-1/12)*fxn(2,col-3))/(dr^2);
d2fdr2(2,col)=((-1/12)*fxn(2,col)+(4/3)*fxn(2,col)+(-5/2)*fxn(2,col)+(4/3)*fxn(2,col-1)+(-1/12)*fxn(2,col-2))/(dr^2);                   
d2fdr2(1,col-1)=((-1/12)*fxn(1,col)+(4/3)*fxn(1,col)+(-5/2)*fxn(1,col-1)+(4/3)*fxn(1,col-2)+(-1/12)*fxn(1,col-3))/(dr^2);           
d2fdr2(1,col)=((-1/12)*fxn(1,col)+(4/3)*fxn(1,col)+(-5/2)*fxn(1,col)+(4/3)*fxn(1,col-1)+(-1/12)*fxn(1,col-2))/(dr^2);                       
                   
% (max,min)

%  Zcomponent(row,1)=  (-1*fxn(row-3,1)+4*fxn(row-2,1)-5* ...
%                         fxn(row-1,1)+2*fxn(row,1))/(dz^2);    
% 
% Zcomponent(row,2)=   (-1*fxn(row-3,2)+4*fxn(row-2,2)-5* ...
%                         fxn(row-1,2)+2*fxn(row,2))/(dz^2); 
%                     
% Zcomponent(row-1,1) = (fxn(row,1) - 2*fxn(row-1,1) ...
%                            +fxn(row-2,1) )/dz^2;
% 
% Zcomponent(row-1,2) = (fxn(row,2) - 2*fxn(row-1,2)+fxn(row-2,2) )/dz^2;
%GP
Zcomponent(row-1,2) = ((-1/12)*fxn(row,2) +(4/3)*fxn(row,2) +(-5/2)*fxn(row-1,2)+(4/3)*fxn(row-2,2)+(-1/12)*fxn(row-3,2) )/dz^2;
Zcomponent(row-1,1) = ((-1/12)*fxn(row,1) +(4/3)*fxn(row,1) +(-5/2)*fxn(row-1,1)+(4/3)*fxn(row-2,1)+(-1/12)*fxn(row-3,1) )/dz^2;
Zcomponent(row,2) = ((-1/12)*fxn(row,2) +(4/3)*fxn(row,2) +(-5/2)*fxn(row,2)+(4/3)*fxn(row-1,2)+(-1/12)*fxn(row-2,2) )/dz^2;
Zcomponent(row,1) = ((-1/12)*fxn(row,1) +(4/3)*fxn(row,1) +(-5/2)*fxn(row,1)+(4/3)*fxn(row-1,1)+(-1/12)*fxn(row-2,1) )/dz^2;

% % (max,max)
% plot 6 corners
%  Zcomponent(row,col)= (-1*fxn(row-3,col)+4*fxn(row-2,col)-5* ...
%                             fxn(row-1,col)+2*fxn(row,col))/(dz^2);    
% 
% Zcomponent(row,col-1)=(-1*fxn(row-3,col-1)+4*fxn(row-2,col-1)-5* ...
%                            fxn(row-1,col-1)+2*fxn(row,col-1))/(dz^2); 
%                 
% Zcomponent(row-1,col) = (fxn(row,col) - 2*fxn(row-1,col) ...
%                            +fxn(row-2,col) )/dz^2;
% 
% Zcomponent(row-1,col-1) = (fxn(row,col-1) - 2*fxn(row-1,col-1) ...
%                            +fxn(row-2,col-1) )/dz^2;
 %GP                       
Zcomponent(row-1,col-1) = ((-1/12)*fxn(row,col-1) +(4/3)*fxn(row,col-1) +(-5/2)*fxn(row-1,col-1)+(4/3)*fxn(row-2,col-1)+(-1/12)*fxn(row-3,col-1) )/dz^2;
Zcomponent(row-1,col) = ((-1/12)*fxn(row,col) +(4/3)*fxn(row,col) +(-5/2)*fxn(row-1,col)+(4/3)*fxn(row-2,col)+(-1/12)*fxn(row-3,col) )/dz^2;                       
Zcomponent(row,col-1) = ((-1/12)*fxn(row,col-1) +(4/3)*fxn(row,col-1) +(-5/2)*fxn(row,col-1)+(4/3)*fxn(row-1,col-1)+(-1/12)*fxn(row-2,col-1) )/dz^2;                   
Zcomponent(row,col) = ((-1/12)*fxn(row,col) +(4/3)*fxn(row,col) +(-5/2)*fxn(row,col)+(4/3)*fxn(row-1,col)+(-1/12)*fxn(row-2,col) )/dz^2;





laplacianRZ=Rcomponent+Zcomponent;
end
