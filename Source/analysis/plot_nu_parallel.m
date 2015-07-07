function [] = plot_nu_parallel(rho,t,muc,delta,Delta,nupar,c),  
figure(98431); hold on;
x = t.*Delta.^nupar;
y = rho.*t.^delta;  

plot(x,y,[c]);
end