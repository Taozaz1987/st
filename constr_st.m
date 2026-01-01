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
real_x = real(sum(stTFR, 1) * df);
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
