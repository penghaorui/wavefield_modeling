%implementation of wavefield model of elastic media with finite-difference time domain approach 
%grid parameter
dx = 10;
dy = 10;
nx = 150;
ny = 150;

%material parameters
lambda=1.6875e10*ones(2*ny,2*nx);
mu=1.6875e10*ones(2*ny,2*nx);
den=2500*ones(2*ny,2*nx);
den(1:ny,:) = 3500;

%PML thickness
xthick = 20;
ythick = 20;

%damping factor
D0=log(1/0.001)*3*sqrt((lambda(1,1)+2*mu(1,1))/den(1,1)) /2/xthick;
Dx = zeros(2*ny,2*nx);
Dy = zeros(2*ny,2*nx);

for i = 1:2*xthick
    Dx(:,[i,2*nx-(i-1)])=D0 *((2*xthick-(i-1))/(2*xthick))^2;
    Dy([i,2*ny-(i-1)],:)=D0 *((2*ythick-(i-1))/(2*ythick))^2;
end

% time step
dt = 0.001;

% time of delay
t0 = 0.045;

%times of reiterate
T=170;

%peak frequency of Ricker wavelet
f0 = 30;

%order of accuracy
M=4;

if M==1;
    a(1)=1;
end

if M==2
    a(1)=9/8;
    a(2)=-1/24;
end

if M>2    
    for i=1:M
        a(i)=(-1)^(i+1)/(2*i-1);
        for k=1:i-1
            a(i)=a(i)*(2*k-1)^2;
        end
        for k=i+1:M
            a(i)=a(i)*(2*k-1)^2;
        end
        for k=1:i-1
            a(i)=a(i)/abs((2*i-1)^2-(2*k-1)^2);
        end
        for k=i+1:M
            a(i)=a(i)/abs((2*i-1)^2-(2*k-1)^2);
        end
    end
end



vx1_x= zeros(ny,nx);
vx1_y= zeros(ny,nx);
vx2_x= zeros(ny,nx);
vx2_y= zeros(ny,nx);
vx1 =zeros(ny,nx);
vx2 =zeros(ny,nx);

vy1_x= zeros(ny,nx);
vy1_y =zeros(ny,nx);
vy2_x =zeros(ny,nx);
vy2_y =zeros(ny,nx);
vy1 =zeros(ny,nx);
vy2 =zeros(ny,nx);

txx1_x =zeros(ny,nx);
txx1_y =zeros(ny,nx);
txx2_x =zeros(ny,nx);
txx2_y =zeros(ny,nx);
txx1 =zeros(ny,nx);
txx2 =zeros(ny,nx);

tyy1_x =zeros(ny,nx);
tyy1_y =zeros(ny,nx);
tyy2_x =zeros(ny,nx);
tyy2_y =zeros(ny,nx);
tyy1 =zeros(ny,nx);
tyy2 =zeros(ny,nx);

txy1_x =zeros(ny,nx);
txy1_y =zeros(ny,nx);
txy2_x =zeros(ny,nx);
txy2_y =zeros(ny,nx);
txy1 =zeros(ny,nx);
txy2 =zeros(ny,nx);

% iterate through time
for i = 1: T
    
    %source
    vx1_x(ny/2,nx/2) = 0.5*(1-2*pi^2*f0^2*((i*dt-0.045)^2))*exp(-pi^2*f0^2*(i*dt-0.045)^2); 
    vx1_y(ny/2,nx/2) = 0.5*(1-2*pi^2*f0^2*((i*dt-0.045)^2))*exp(-pi^2*f0^2*(i*dt-0.045)^2);
    vx1(ny/2,nx/2) = vx1_x(ny/2,nx/2) + vx1_y(ny/2,nx/2);
    
    for y = M+1:ny-M
        for x = M+1:nx-M 
            
            vx2_x(y,x) = vx1_x(y,x) * (1 - 0.5*dt*Dx(2*(y-1)+1, 2*(x-1)+1)) / (1 + 0.5*dt*Dx(2*(y-1)+1, 2*(x-1)+1)); 
            vx2_y(y,x) = vx1_y(y,x) * (1 - 0.5*dt*Dy(2*(y-1)+1, 2*(x-1)+1)) / (1 + 0.5*dt*Dy(2*(y-1)+1, 2*(x-1)+1));
            
            vy2_x(y,x) = vy1_x(y,x) * (1 - 0.5*dt*Dx(2*y, 2*x)) / (1 + 0.5*dt*Dx(2*y, 2*x)); 
            vy2_y(y,x) = vy1_y(y,x) * (1 - 0.5*dt*Dy(2*y, 2*x)) / (1 + 0.5*dt*Dy(2*y, 2*x));
            
            for m = 1:M
                vx2_x(y,x) = vx2_x(y,x) + dt/dx/den(2*(y-1)+1,2*(x-1)+1) / (1 + 0.5*dt*Dx(2*(y-1)+1, 2*(x-1)+1)) * a(m)*(txx1(y,x+m-1)- txx1(y,x-m));               
                vx2_y(y,x) = vx2_y(y,x) + dt/dy/den(2*(y-1)+1,2*(x-1)+1) / (1 + 0.5*dt*Dy(2*(y-1)+1, 2*(x-1)+1)) * a(m)*(txy1(y+m-1,x)- txy1(y-m,x));      
                
                vy2_x(y,x) = vy2_x(y,x) + dt/dx/den(2*y,2*x) / (1 + 0.5*dt*Dx(2*y, 2*x)) * a(m)*(txy1(y,x+m)- txy1(y,x-(m-1)));
                vy2_y(y,x) = vy2_y(y,x) + dt/dy/den(2*y,2*x) / (1 + 0.5*dt*Dy(2*y, 2*x)) * a(m)*(tyy1(y+m,x)- tyy1(y-(m-1),x)); 
            end
            vx2(y,x) = vx2_x(y,x) + vx2_y(y,x);
            vy2(y,x) = vy2_x(y,x) + vy2_y(y,x); 
        end
    end
    
    for y = M+1:ny-M
        for x = M+1:nx-M
            
            txx2_x(y,x) = txx1_x(y,x) * (1 - 0.5*dt*Dx(2*(y-1)+1, 2*x)) / (1 + 0.5*dt*Dx(2*(y-1)+1, 2*x));
            txx2_y(y,x) = txx1_y(y,x) * (1 - 0.5*dt*Dy(2*(y-1)+1, 2*x)) / (1 + 0.5*dt*Dy(2*(y-1)+1, 2*x));
            
            tyy2_x(y,x) = tyy1_x(y,x) * (1 - 0.5*dt*Dx(2*(y-1)+1, 2*x)) / (1 + 0.5*dt*Dx(2*(y-1)+1, 2*x));
            tyy2_y(y,x) = tyy1_y(y,x) * (1 - 0.5*dt*Dy(2*(y-1)+1, 2*x)) / (1 + 0.5*dt*Dy(2*(y-1)+1, 2*x));
            
            txy2_x(y,x) = txy1_x(y,x) * (1 - 0.5*dt*Dx(2*y, 2*(x-1)+1)) / (1 + 0.5*dt*Dx(2*y, 2*(x-1)+1));
            txy2_y(y,x) = txy1_y(y,x) * (1 - 0.5*dt*Dy(2*y, 2*(x-1)+1)) / (1 + 0.5*dt*Dy(2*y, 2*(x-1)+1));
            
            for m = 1:M
                txx2_x(y,x) =txx2_x(y,x) + dt/dx*(lambda(2*(y-1)+1,2*x) + 2*mu(2*(y-1)+1,2*x)) / (1 + 0.5*dt*Dx(2*(y-1)+1, 2*x))* a(m)*(vx2(y,x+m)- vx2(y,x-(m-1)));
                txx2_y(y,x) =txx2_y(y,x) + dt/dy*lambda(2*(y-1)+1,2*x) / (1 + 0.5*dt*Dy(2*(y-1)+1, 2*x)) * a(m)*(vy2(y+m-1,x)- vy2(y-m,x));  
                       
                tyy2_x(y,x) =tyy2_x(y,x) + dt/dx*lambda(2*(y-1)+1,2*x) / (1 + 0.5*dt*Dx(2*(y-1)+1, 2*x)) * a(m)*(vx2(y,x+m)- vx2(y,x-(m-1)));
                tyy2_y(y,x) =tyy2_y(y,x) + dt/dy*(lambda(2*(y-1)+1,2*x) + 2*mu(2*(y-1)+1,2*x)) / (1 + 0.5*dt*Dy(2*(y-1)+1, 2*x))* a(m)*(vy2(y+m-1,x)- vy2(y-m,x));
                       
                txy2_x(y,x) =txy2_x(y,x) + dt/dx*mu(2*y,2*(x-1)+1) / (1 + 0.5*dt*Dx(2*y,2*(x-1)+1)) * a(m)*(vy2(y,x+m-1)- vy2(y,x-m));
                txy2_y(y,x) =txy2_y(y,x) + dt/dy*mu(2*y,2*(x-1)+1) / (1 + 0.5*dt*Dy(2*y,2*(x-1)+1)) * a(m)*(vx2(y+m,x)- vx2(y-(m-1),x)); 
            end
            txx2(y,x) = txx2_x(y,x) + txx2_y(y,x);
            tyy2(y,x) = tyy2_x(y,x) + tyy2_y(y,x);
            txy2(y,x) = txy2_x(y,x) + txy2_y(y,x);
        end
    end
    
    vx1_x = vx2_x;
    vx1_y = vx2_y;
    vx1 = vx2; 
    
    vy1_x = vy2_x;
    vy1_y = vy2_y;
    vy1 = vy2;
    
    txx1_x = txx2_x;
    txx1_y = txx2_y;
    txx1 = txx2;
    
    tyy1_x= tyy2_x;
    tyy1_y = tyy2_y;
    tyy1 = tyy2;
    
    txy1_x = txy2_x;
    txy1_y = txy2_y;
    txy1 = txy2;
end
imagesc(vx1);figure(gcf);
hold on
