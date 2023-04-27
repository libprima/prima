function [tlist, plist] = timing(maxdim)

req.mindim = 1;
req.maxdim = maxdim;
req.maxcon = 500*maxdim;
req.type='ubln';
plist = secup(req);
tlist= NaN(length(plist), 1);

for ip = 1 : length(plist)
    plist{ip}
    tic;
    cobyla(macup(plist{ip}));
    tlist(ip) = toc;
    tlist(ip)
end

save('timing.mat', 'plist', 'tlist');
