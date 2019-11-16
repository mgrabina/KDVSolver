%lie trotter
%Kdv Dos_solitones

clear all
clc
set(gca,'FontSize',8)
set(gca,'LineWidth',2)

N = 256;
x = linspace(-10,10,N);
delta_x = x(2) - x(1);
delta_k = 2*pi/(N*delta_x);

k = [0:delta_k:(N/2-1)*delta_k,0,-(N/2-1)*delta_k:delta_k:-delta_k];
c_1=13;
c_2 =3;

u = 1/2*c_1*(sech(sqrt(c_1)*(x+8)/2)).^2 + 1/2*c_2*(sech(sqrt(c_2)*(x+1)/2)).^2;

delta_t = 0.4/N^2;
t=0;
plot(x,u,'LineWidth',2)
axis([-10 10 0 10])
xlabel('x')
ylabel('u')
text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',14)
drawnow

tmax = 1.5; nplt = floor((tmax/100)/delta_t); nmax = round(tmax/delta_t);
udata = u.'; tdata = 0;
U1 = fft(u);
U2 = fft(u);
for n = 1:nmax-40000
    t = n*delta_t;
    
    % lineal
    U1 = U1.*exp(1i*k.^3*delta_t);
    
    % no lineal
    U1 = U1  - (3i*k*delta_t).*fft((real(ifft(U1))).^2);
    
    % no lineal
    U2 = U2 - (3i*k*delta_t).*fft((real(ifft(U2))).^2);
    % lineal
    U2 = U2.*exp(1i*k.^3*delta_t);
    
    U = -U1/6 + U2*(2/3);
    
    if mod(n,nplt) == 0
        u = real(ifft(U));
        udata = [udata u.']; tdata = [tdata t];
        if mod(n,4*nplt) == 0
            plot(x,u,'LineWidth',2)
            axis([-10 10 0 10])
            xlabel('x')
            ylabel('u')
            text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',10)
            drawnow
        end
    end
end

figure

waterfall(x,tdata(1:4:end),udata(:,1:4:end)')
xlabel x, ylabel t, axis([-10 10 0 tmax 0 10]), grid off
zlabel u