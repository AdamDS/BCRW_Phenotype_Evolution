critical = 0.337;
D = mutability'-critical;
lD = log10(abs(D));
lAP = log10(AVG_POPS);
lSP = log10(STD_POPS);
lower = find(D(4:end)<0);
upper = find(D>0);
[mlower,blower,smlower,sblower] = linear_fit(lD(lower),lAP(lower),[],1);
[mupper,bupper,smupper,sbupper] = linear_fit(lD(upper),lAP(upper),[],1);
[mlower,blower,smlower,sblower] = linear_fit(lD(lower),lSP(lower),[],1);
[mupper,bupper,smupper,sbupper] = linear_fit(lD(upper),lSP(upper),[],1);
figure, plot(lD,lAP);