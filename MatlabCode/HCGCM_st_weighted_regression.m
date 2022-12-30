function[beta,r2]=HCGCM_st_weighted_regression(SST,ATOM,w)



w(w<0.01)=nan;

SST=SST.*sqrt(w);
ATOM=ATOM.*sqrt(w);

SST=SST(:);
ATOM=ATOM(:);
w=w(:);

index=logical((~isnan(ATOM)).*(~isnan(SST)));
SST=SST(index);
ATOM=ATOM(index);
w=w(index);



if length(SST)<50
    beta=nan(1,3);
    r2=nan;
else
    [b,bint,r]=regress(ATOM,[SST,sqrt(w)]);
    beta=[b(1),bint(1,:)];
    r2=1-sum(r.^2)/sum(ATOM.^2);
end



