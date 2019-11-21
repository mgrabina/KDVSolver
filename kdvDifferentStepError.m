
%Kdv Dos_solitones
% compares step1 = 0.4/256^² with step 2 = (0.4/256^²) * 1 / stepDifference
% in order '4' againts order 6
function kdvDifferentStepError(stepDifference)
kdvOrder = 4;
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
delta_t2 = 0.4/(N^2*stepDifference);

t=0;
plot(x,u,'LineWidth',2)
axis([-10 10 0 10])
xlabel('x')
ylabel('u')
text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',14)
drawnow

tmax = 1.5; nplt = floor((tmax/100)/delta_t); nmax = round(tmax/delta_t);
udata = u.'; tdata = 0;
for i=1:1:kdvOrder
    U{i} = fft(u);
end

for i=1:1:kdvOrder
    U2{i} = fft(u);
end


for i=1:1:6
    U3{i} = fft(u);
end

% orders 2 , 4 , 6
orders = {[-1/6,2/3], [1/90,-2/9,0,32/45], [-1/1680,1/15,-27/80,0,0,27/35]};
map = [2, 2, 4, 4, 6, 6];
order = orders{1,kdvOrder/2};
order2 = orders{1,3};
time = 1;
for n = 1:nmax-40000
    
    t = n*delta_t;
    % step 1
    for i = 1:1:kdvOrder
        U{i} = calculateOrder(delta_t/ceil(i/2),k,U{i},map(i),i);
    end
    retU = 0;
    for i=1:1:kdvOrder
        retU = retU + 2 * order(i) * U{i};
    end
    
    
    % step 1 order 6
    retU3 = 0;
    for i = 1:1:6
        U3{i} = calculateOrder(delta_t/ceil(i/2),k,U3{i},map(i),i);
    end
    for i=1:1:6
        retU3 = retU3 + 2 * order2(i) * U3{i};
    end
    
    % step 2
    retU2 = 0;
    for z = 1:stepDifference
        for i = 1:1:kdvOrder
            U2{i} = calculateOrder(delta_t2/ceil(i/2),k,U2{i},map(i),i);
        end
        retU2 = 0;
        for i=1:1:kdvOrder
            retU2 = retU2 + 2 * order(i) * U2{i};
        end
    end
    
    if mod(n,nplt) == 0
        u = real(ifft(retU));
        u2 = real(ifft(retU2));
        u3 = real(ifft(retU3));
        udata = [udata u.']; tdata = [tdata t];
        if mod(n,4*nplt) == 0
            subplot(1,2,1)
            t2(time) = t;
            error1(time) = mean(abs(u3 - u));
            error2(time) = mean(abs(u3 - u2));
            plot(x,u,'LineWidth',2)
            axis([-10 10 0 10])
            xlabel('x')
            ylabel('u')
            text(6,9,['t = ',num2str(t,'%1.2f')],'FontSize',10)
            subplot(1,2,2)
            plot(t2,error1,"b*");
            hold on;
            plot(t2,error2,"b-");
            xlabel('time[s]')
            legend('Decremented step')
            ylabel('Mean Error')
            hold off;
            time = time + 1;
            drawnow
        end
    end
end

figure

waterfall(x,tdata(1:4:end),udata(:,1:4:end)')
xlabel x, ylabel t, axis([-10 10 0 tmax 0 10]), grid off
zlabel u
end

function ret=linear(delta_t,k,U)
ret = U.*exp(1i*k.^3*delta_t);
end

function ret=nonlinear(delta_t,k,U)
ret = U - (3i*k*delta_t).*fft((real(ifft(U))).^2);
end

function ret=calculateOrder(delta_t,k,U,order,index)
ret = U;
if(mod(index,2)==0)
    for i=1:1:order/2
        ret = nonlinear(delta_t,k,ret);
        ret = linear(delta_t,k,ret);
    end
else
    for i=1:1:order/2
        ret = linear(delta_t,k,ret);
        ret = nonlinear(delta_t,k,ret);
    end
end
end




