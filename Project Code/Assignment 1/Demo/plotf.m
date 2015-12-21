clear all
godunov
eoscheme

plot(x,ugod,x,ueo,x,u0g,x,u0e);
title('At t=5');
axis([0,5,0,1])
legend('Godunov','Einquist Osher','Godunov at t=0','EO at t=0');