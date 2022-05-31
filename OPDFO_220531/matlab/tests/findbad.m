req.maxdim = 20;
req.maxcon = 20;
req.type='ubln';
plist=secup(req);
np=length(plist);
for ip=1:np
    p=macup(plist{ip});
    tic
    [a1, b1,c1, d1] = cobyla(p);
    t1=toc
    tic
    [a2, b2, c2, d2] = cobyla(p);
    t2=toc
    if t1/t2 >= 5*d1.funcCount/d2.funcCount
        p
        break
    end
end
