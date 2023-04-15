taus = 1000*[1:5];

CtimeRes = zeros(size(taus));
CaerrRes = zeros(size(taus));
CdRes = zeros(size(taus));

CtimeIms = zeros(size(taus));
CaerrIms = zeros(size(taus));
CdIms = zeros(size(taus));

for ii=1:numel(taus)
    tau = taus(ii);
    fprintf(1, '---------tau: %d\n', tau);
    N = max(40*tau,1024);    N = 2*ceil(N/2);
    
    %Hamiltonian simulation (real part)
    ts = [0:(N-1)]'/N * 2*pi;
    zs = exp(i*ts);
    xs = cos(ts);
    
    aext = 0.3 * cos(tau*xs);
    
    aux = fft(aext);
    tmp = abs(aux(1:end/2));
    d = max(find(abs(tmp)>1e-12*max(tmp)));
    if(mod(d,2)==1) d=d+1; end; %make even
    ks = [0:(N/2-1) -N/2:-1]';    aux = aux.*(abs(ks)<=d);
    as = ifft(aux);    fprintf(1, 'as error %1.1d\n', norm(as-aext)/norm(aext));
    
    D = 2*d+1;
    L = 2*d+2;
    
    bs = 0.4*sin(d*ts)+ 0*sin((d-2)*ts);
    
    nullAlg = 'itr';
    gall;
    
    CdRes(ii) = d;
    CtimeRes(ii) = Ctime;
    CaerrRes(ii) = Caerr;
    
    %Hamiltonian simulation (imag part)
    ts = [0:(N-1)]'/N * 2*pi;
    zs = exp(i*ts);
    xs = cos(ts);
    
    aext = 0.3 * sin(tau*xs);
    
    aux = fft(aext);
    tmp = abs(aux(1:end/2));
    d = max(find(abs(tmp)>1e-12*max(tmp)));
    if(mod(d,2)==0) d=d+1; end; %make odd
    ks = [0:(N/2-1) -N/2:-1]';    aux = aux.*(abs(ks)<=d);
    as = ifft(aux);    fprintf(1, 'as error %1.1d\n', norm(as-aext)/norm(aext));
    
    D = 2*d+1;
    L = 2*d+2;
    
    bs = 0.4*sin(d*ts)+ 0*sin((d-2)*ts);
    
    nullAlg = 'itr';
    gall
    
    CdIms(ii) = d;
    CtimeIms(ii) = Ctime;
    CaerrIms(ii) = Caerr;
end
save tst1.mat taus    CtimeRes    CaerrRes    CdRes    CtimeIms    CaerrIms    CdIms


    
FS = 20;
if(1)
    figure(1);
    loglog(taus,CdRes,'r-+'); hold on;grid on;
    loglog(taus,CdIms,'b-+');
    legend('Re', 'Im');
    xlabel('\tau');
    ylabel('d');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst1_d');
    
    figure(2);
    loglog(taus,CtimeRes,'r-+'); hold on;grid on;
    loglog(taus,CtimeIms,'b-+');
    legend('Re', 'Im');
    xlabel('\tau');
    ylabel('time(sec)');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst1_t');
    
    figure(3);
    loglog(taus,CaerrRes,'r-+'); hold on;grid on;
    loglog(taus,CaerrIms,'b-+');
    legend('Re', 'Im');
    xlabel('\tau');
    ylabel('error');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst1_e');
end

if(1)
    tau = 25;
    fprintf(1, '---------tau: %d\n', tau);
    N = max(40*tau,1024);    N = 2*ceil(N/2);
    ts = [0:(N-1)]'/N * 2*pi;
    zs = exp(i*ts);
    xs = cos(ts);
    
    %Hamiltonian simulation (real part)
    aext = 0.3 * cos(tau*xs);
    
    aux = fft(aext);
    tmp = abs(aux(1:end/2));
    d = max(find(abs(tmp)>1e-12*max(tmp)));
    if(mod(d,2)==1) d=d+1; end; %make even
    ks = [0:(N/2-1) -N/2:-1]';    aux = aux.*(abs(ks)<=d);
    as = ifft(aux);    fprintf(1, 'as error %1.1d\n', norm(as-aext)/norm(aext));
    asRe = as;
    
    aext = 0.3 * sin(tau*xs);
    
    aux = fft(aext);
    tmp = abs(aux(1:end/2));
    d = max(find(abs(tmp)>1e-12*max(tmp)));
    if(mod(d,2)==0) d=d+1; end; %make odd
    ks = [0:(N/2-1) -N/2:-1]';    aux = aux.*(abs(ks)<=d);
    as = ifft(aux);    fprintf(1, 'as error %1.1d\n', norm(as-aext)/norm(aext));
    asIm = as;
    
    figure(4);
    plot(xs,asRe,'b-'); hold on;grid on;
    plot(xs,asIm,'r-');
    legend('Re', 'Im');
    xlabel('x');
    ylabel('a(x)');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst1_a');

end


