clear;clc;
%% 测试信号
t=linspace(0,1,1000);
%s=2*sin(2*pi*20*t)+1.5*cos(sin(pi*5*t));
dt=t(2)-t(1);
s=2.*cos(2.*pi.*(80.*t))+(1+0.5.*cos(2.*t)).*exp(-t/10).*cos(10.*pi.*(8.*t+6.*t.^2)...
    +0.3.*cos(t))+(2+0.2.*cos(t)).*sin(10*pi.*(5.*t+0.3.*cos(6.*t)));
[stTFR,f]= ST(s, t);
df=f(2)-f(1);
magnitude_s=abs(stTFR);
x=zeros(1,length(s));
for t_index=1:length(t )
    for f_index=1:length(f )
        for tau=1:length(t )
            x(t_index)=stTFR(f_index,tau)*exp(1i*2*pi*f(f_index)...
                *t(t_index))*dt+x(t_index);
        end        
    end
    x(t_index)=x(t_index)*df;
end
real_x=real(x);
figure(1);
p1=plot(t,s,'b-');
hold on;
p2=plot(t,real_x,'r--');
legend([p1, p2], {'原始信号', '重构信号'}, ...
       'Location', 'northeast', ...  % 右上角
       'FontSize', 12);
xlabel('时间');
ylabel('幅值');
title('信号对比图');
grid on;
figure(2);
imagesc(t,f,magnitude_s');
axis xy;
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('ST');
