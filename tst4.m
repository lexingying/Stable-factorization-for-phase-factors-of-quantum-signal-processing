betas = 100*[1 2 4 8 16];

Ctimes = zeros(size(betas));
Caerrs = zeros(size(betas));
Cds = zeros(size(betas));

for ii=1:numel(betas)
    beta = betas(ii);
    fprintf(1, '---------beta: %d\n', beta);
    
    N = max(200*beta,1024);    N = 2*ceil(N/2);
    ts = [0:(N-1)]'/N * 2*pi;
    zs = exp(i*ts);
    xs = cos(ts);
    
    aext = 1./(1+exp(beta*xs)) - 0.5;
    aext = aext/max(abs(aext)) * 0.3;
    
    %Hamiltonian simulation
    aux = fft(aext);
    tmp = abs(aux(1:end/2));
    d = max(find(abs(tmp)>1e-12*max(tmp)));
    if(mod(d,2)==0) d=d+1; end;
    ks = [0:(N/2-1) -N/2:-1]';    aux = aux.*(abs(ks)<=d);
    as = ifft(aux);    fprintf(1, 'as error %1.1d\n', norm(as-aext)/norm(aext));
    
    %L = 2*d+2;
    D = 2*d+1;
    L = 2*d+2;
    
    bs = 0.4*sin(d*ts);

    nullAlg = 'itr';
    gall

    Cds(ii) = d;
    Ctimes(ii) = Ctime;
    Caerrs(ii) = Caerr;
end
save tst4.mat betas Ctimes Caerrs Cds



FS = 20;
if(1)
    figure(1);
    loglog(betas,Cds,'b-+'); hold on;grid on;
    xlabel('\beta');
    ylabel('d');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst4_d');
    
    figure(2);
    loglog(betas,Ctimes,'b-+'); hold on;grid on;
    xlabel('\beta');
    ylabel('time(sec)');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst4_t');
    
    figure(3);
    loglog(betas,Caerrs,'b-+'); hold on;grid on;
    xlabel('\beta');
    ylabel('error');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst4_e');
end

if(1)
    beta = 100;
    fprintf(1, '---------beta: %d\n', beta);
    
    N = max(200*beta,1024);    N = 2*ceil(N/2);
    ts = [0:(N-1)]'/N * 2*pi;
    zs = exp(i*ts);
    xs = cos(ts);
    
    aext = 1./(1+exp(beta*xs)) - 0.5;
    aext = aext/max(abs(aext)) * 0.3;
    
    %Hamiltonian simulation
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
    print(gcf, '-depsc', 'tst4_a');

end

