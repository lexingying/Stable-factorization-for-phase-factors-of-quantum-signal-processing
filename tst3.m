%kps = 10*[1 2 4 8 16 32];
%kps = [1000];

kps = [16 64 256 1024]

Ctimes = zeros(size(kps));
Caerrs = zeros(size(kps));
Cds = zeros(size(kps));

for ii=1:numel(kps)
    kp = kps(ii);
    fprintf(1, '---------kp: %d\n', kp);
    
    N = max(1000*kp,1024);    N = 2*ceil(N/2);
    ts = [0:(N-1)]'/N * 2*pi;
    zs = exp(i*ts);
    xs = cos(ts);

    wfn = @(x) 1-exp(-(5*x*kp).^2);
    aext = 1./xs;     bad = find(abs(xs)<1e-10);    aext(bad) = 0;
    aext = aext.*wfn(xs);
    aext = aext/max(abs(aext)) * 0.3;
    
    aux = fft(aext);
    tmp = abs(aux(1:end/2));
    d = max(find(abs(tmp)>1e-12*max(tmp)));
    if(mod(d,2)==0) d=d+1; end;
    ks = [0:(N/2-1) -N/2:-1]';    aux = aux.*(abs(ks)<=d);
    as = ifft(aux);    fprintf(1, 'as error %1.1d\n', norm(as-aext)/norm(aext));
    
    D = 2*d+1;
    L = 2*d+2;
    
    bs = 0.3*sin(d*ts);
    
    nullAlg = 'fft';
    gall
    
    Cds(ii) = d;
    Ctimes(ii) = Ctime;
    Caerrs(ii) = Caerr;
end
save tst3.mat kps Ctimes Caerrs Cds


FS = 20;
if(1)
    figure(1);
    loglog(kps,Cds,'b-+'); hold on;grid on;
    xlabel('\kappa');
    ylabel('d');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst3_d');
    
    figure(2);
    loglog(kps,Ctimes,'b-+'); hold on;grid on;
    xlabel('\kappa');
    ylabel('time(sec)');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst3_t');
    
    figure(3);
    loglog(kps,Caerrs,'b-+'); hold on;grid on;
    xlabel('\kappa');
    ylabel('error');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst3_e');
end


if(1)
    kp =10;
    fprintf(1, '---------kp: %d\n', kp);
    
    N = max(1000*kp,1024);    N = 2*ceil(N/2);
    ts = [0:(N-1)]'/N * 2*pi;
    zs = exp(i*ts);
    xs = cos(ts);
    
    wfn = @(x) 1-exp(-(5*x*kp).^2);
    aext = 1./xs;     bad = find(abs(xs)<1e-10);    aext(bad) = 0;
    aext = aext.*wfn(xs);
    aext = aext/max(abs(aext)) * 0.3;
    
    aux = fft(aext);
    tmp = abs(aux(1:end/2));
    d = max(find(abs(tmp)>1e-11*max(tmp)));
    if(mod(d,2)==0) d=d+1; end;
    ks = [0:(N/2-1) -N/2:-1]';    aux = aux.*(abs(ks)<=d);
    as = ifft(aux);    fprintf(1, 'as error %1.1d\n', norm(as-aext)/norm(aext));
    
    figure(4);
    plot(xs,as,'b-'); hold on;grid on;
    xlabel('x');
    ylabel('a(x)');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst3_a');

end