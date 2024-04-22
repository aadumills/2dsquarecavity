

%code by Abraham Adu-Mills
%remove comment from line 109 to end to plot results
clear
clc
dpdx=-45;      %pressure gradient

dz=0.005;   %dz
dy=0.005;   %dy
b=dy/dz;    % beta
L=1.5;      %Length of domain
h=1;        %height of domain
mu=0.4;     %dynamic viscosity
coeff=1/(2*(1+b^2));
ystar=L/L;
zstar=h/L;

JM=(ystar/dy)+1;  %number of points in the y direction
KM=round((zstar/dz))+1;  %number of points in the z direction
y=linspace(0,L,JM);
z=linspace(0,h,KM);
ustar=zeros(KM,JM,2);    %data for u* using gauss-seidel iteration method
ustarpsor=zeros(KM,JM,2);
uexact=zeros(KM,JM);
% D=linspace(0,0,199);
m=1;
%exact solution
A=(16*L^2)/(mu*pi^3);
B=-dpdx;
while m<500
    C=(-1)^((m-1)/2);
    
    D=1-cosh(m*pi*z(2:KM-1)./(2*L))./cosh(m*pi*h/(2*L));
    E=cos(m*pi*y(2:JM-1)./(2*L))./m^3;
  
        uexact(2:KM-1,2:JM-1)=uexact(2:KM-1,2:JM-1)+A.*B.*C.*D'.*E;
     
    
    m=m+2;
end

uexact(1,:,:)=uexact(2,:,:);
uexact(:,1,:)=uexact(:,2,:);
uexact(end,:,:)=0;
uexact(:,end,:)=0;


%%%%%%%%%%%%%%%%%%%%%%%



omg=1;


ustar(2:KM-1,2:JM-1,1)=0.1;  %initial guess for gauss-seidel iteration method
ustar(2:KM-1,2:JM-1,2)=0.15; 
ustarpsor(2:KM-1,2:JM-1,1)=0.1;  %initial guess for point successive over relaxation method 
ustarpsor(2:KM-1,2:JM-1,2)=0.15; 
       


it=0;  %number of iterations
%Apply point Gauss-Seidel Iteration Method
while any(max(abs((ustar(:,:,2)-ustar(:,:,1))))>=1e-6)
ustar(:,:,1)=ustar(:,:,2);
ustar(2:KM-1,2:JM-1,2) = coeff*(ustar(2:KM-1,3:JM,1)+ustar(2:KM-1,1:JM-2,2)+(b^2)*(ustar(3:KM,2:JM-1,1)+ustar(1:KM-2,2:JM-1,2))+(dy^2));

%apply boundary conditions for gauss-seidel


ustar(1,:,:)=ustar(2,:,:);
ustar(:,1,:)=ustar(:,2,:);
ustar(end,:,:)=0;
ustar(:,end,:)=0;


it=it+1; %increase number of iteration

end
it2=0;
tic
while any(max(abs((ustarpsor(:,:,2)-ustarpsor(:,:,1))))>=1e-06)
ustarpsor(:,:,1)=ustarpsor(:,:,2);

ustarpsor(2:KM-1,2:JM-1,2) =(1-omg)*ustarpsor(2:KM-1,2:JM-1,1)+ omg*coeff*(ustarpsor(2:KM-1,3:JM,1)+ustarpsor(2:KM-1,1:JM-2,2)+(b^2)*(ustarpsor(3:KM,2:JM-1,1)+ustarpsor(1:KM-2,2:JM-1,2))+(dy^2));


%apply boundary conditions for psor
ustarpsor(1,:,:)=ustarpsor(2,:,:);
ustarpsor(:,1,:)=ustarpsor(:,2,:);
ustarpsor(end,:,:)=0;
ustarpsor(:,end,:)=0;



it2=it2+1; %increase number of iteration

end
toc

ugs=ustar(:,:,2)*L^2*(-dpdx)/mu;
upsor=ustarpsor(:,:,2)*L^2*(-dpdx)/mu;
error_gs=abs((ugs(2:KM-1,2:JM-1)-uexact(2:KM-1,2:JM-1))./uexact(2:KM-1,2:JM-1))*100;
error_psor=abs((upsor(2:KM-1,2:JM-1)-uexact(2:KM-1,2:JM-1))./uexact(2:KM-1,2:JM-1))*100;



% %plot the error distribution for gauss-seidel method
% figure(1)
% hold on
% surf(y(2:JM-1),z(2:KM-1),error_gs)
%  title('Error distribution (Point Gauss-Seidel Method')
% xlabel('y'),zlabel('error%'),ylabel('z')
% grid on
% mesh(y(2:JM-1),z(2:KM-1),error_gs)
%  % title('Error distribution (Point Gauss-Seidel Method')
% xlabel('y'),zlabel('error%'),ylabel('z')
% colorbar
% hold off
% 
% %plot the gauss-seidel solution
figure(2)
grid on
hold on
mesh(y,z,ugs)
colorbar
 % title('Two dimensional flow field (Point Gauss-Seidel Method')
xlabel('y'),zlabel('velocity (ms^-1)'),ylabel('z')

% 
% mesh(y,z,ugs)
%  title('Two dimensional flow field (Point Gauss-Seidel Method')
% xlabel('y'),zlabel('velocity (ms^-1)'),ylabel('z')
% colorbar
% hold off

%plot the solution point succesive over relaxation method

% surf(y,z,upsor)
%  title('Two dimensional flow field (PSOR Method: \omega =1)')
% xlabel('y'),zlabel('velocity (ms^-1)'),ylabel('z')
% figure(3)
% hold on
% grid on
%  mesh(y,z,upsor)
%  % title('Two dimensional flow field (PSOR Method: \omega =1)')
% xlabel('y'),zlabel('velocity (ms^-1)'),ylabel('z')
% colorbar
% hold off

%plot the error for PSOR method
% figure(4)
% hold on
% grid on
% mesh(y(2:JM-1),z(2:KM-1),error_psor)
% colorbar
% % title('Error distribution (PSOR Method)')
% xlabel('y'),zlabel('error%'),ylabel('z')

% 
% mesh(y(2:JM-1),z(2:KM-1),error_psor)
%  title('Error distribution (PSOR Method)')
% xlabel('y'),zlabel('error%'),ylabel('z')
% colorbar
% hold off
% 
% 
