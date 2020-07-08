function [yPred, tPred, m] = fitGrowth2(Y,t)
%define some other params

yPred=0;
tPred=0;
m=0;

if min(Y)>0
    yStar=2*min(Y);
    yPrime=0.6*max(Y);

    tStar=t(sum(Y<yStar));
    tPrime=t(max(1,sum(Y<yPrime)));

    %xform Y into log space
    yStar=log10(yStar);
    yPrime=log10(yPrime);

    %calculate slope in log space
    m=(yPrime-yStar)/(tPrime-tStar);

    %generate predicted curve
    tPred=t(t>tStar);
    tPred=tPred(tPred<tPrime);

    tAdj=tPred-tStar;
    yPred=yStar+m.*tAdj;

    yPred=10.^yPred;
end

end


