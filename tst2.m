%deltas = 0.005 * [16 8 4 2 1];
deltas = 0.005 * [16 8 4 2 1];
%deltas = 0.01;

Ctimes = zeros(size(deltas));
Caerrs = zeros(size(deltas));
Cds = zeros(size(deltas));

for ii = 1:numel(deltas)
    delta = deltas(ii);
    fprintf(1, '---------delta: %d\n', delta);
    
    k = 20/delta;
    N = max(100*k,1024);    N = 2*ceil(N/2);
    ts = [0:(N-1)]'/N * 2*pi;
    zs = exp(i*ts);
    xs = cos(ts);
    
    pos = -1 + 2*(xs.^2-delta^2)/(1-delta^2);
    aext = zeros(size(pos));
    gud = find((abs(pos)<=1));
    bad = find((abs(pos)>1));
    aext(gud) = cos(k*acos(pos(gud)));
    aext(bad) = chebyshevT(k,pos(bad));
    
    aext = aext/max(abs(aext)) * 0.3;
    
    aux = fft(aext);
    tmp = abs(aux(1:end/2));
    d = max(find(abs(tmp)>10e-12*max(tmp)));    %d = 2*k;
    if(mod(d,2)==1) d=d+1; end; %make even
    ks = [0:(N/2-1) -N/2:-1]';    aux = aux.*(abs(ks)<=d);
    as = ifft(aux);    fprintf(1, 'as error %1.1d\n', norm(as-aext,Inf)/norm(aext,Inf));
    
    D = 2*d+1;
    L = 2*d+2;
    
    bs = 0.4*sin(d*ts);
    
    nullAlg = 'itr';
    gall
    
    Cds(ii) = d;
    Ctimes(ii) = Ctime;
    Caerrs(ii) = Caerr;
end
save tst2.mat     deltas    Ctimes    Caerrs    Cds

FS = 20;
if(1)
    figure(1);
    loglog(1./deltas,Cds,'b-+'); hold on; grid on;set(gca, 'XLimSpec', 'Tight');
    xlabel('1/\Delta');
    ylabel('d');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst2_d');
    
    figure(2);
    loglog(1./deltas,Ctimes,'b-+'); hold on;grid on;set(gca, 'XLimSpec', 'Tight');
    xlabel('1/\Delta');
    ylabel('time(sec)');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst2_t');
    
    figure(3);
    loglog(1./deltas,Caerrs,'b-+'); hold on;grid on;set(gca, 'XLimSpec', 'Tight');
    xlabel('1/\Delta');
    ylabel('error');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst2_e');
end




if(1)
    delta = 0.005*16;
    fprintf(1, '---------delta: %d\n', delta);
    
    k = 20/delta;
    N = max(100*k,1024);    N = 2*ceil(N/2);
    ts = [0:(N-1)]'/N * 2*pi;
    zs = exp(i*ts);
    xs = cos(ts);
    
    pos = -1 + 2*(xs.^2-delta^2)/(1-delta^2);
    aext = zeros(size(pos));
    gud = find((abs(pos)<=1));
    bad = find((abs(pos)>1));
    aext(gud) = cos(k*acos(pos(gud)));
    aext(bad) = chebyshevT(k,pos(bad));
    
    aext = aext/max(abs(aext)) * 0.3;

    
    aux = fft(aext);
    tmp = abs(aux(1:end/2));
    d = max(find(abs(tmp)>10e-12*max(tmp)));    %d = 2*k;
    if(mod(d,2)==1) d=d+1; end; %make even
    ks = [0:(N/2-1) -N/2:-1]';    aux = aux.*(abs(ks)<=d);
    as = ifft(aux);    %fprintf(1, 'as error %1.1d\n', norm(as-aext,Inf)/norm(aext,Inf));
    
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
    print(gcf, '-depsc', 'tst2_a');

end


