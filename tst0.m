if(0)
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
    
    %L = 2*d+2;
    D = 2*d+1;
    L = 2*d+2;
    
    bs = 0*sin(d*ts);
    
    gall
    
    M1 = M;
    
    bs = 0.3*sin(d*ts);
    
    gall
    
    M2 = M;

    
end

FS = 20;
if(1)
    s = svd(M1);

    figure(1);
    imagesc(M1); colorbar; axis equal; axis tight;
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst0_M1');
    
    figure(2);
    plot(s, '-+');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst0_s1');
    
    
    s = svd(M2);

    figure(3);
    imagesc(M2); colorbar; axis equal; axis tight;
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst0_M2');
    
    figure(4);
    plot(s, '-+');
    set(gca, 'FontSize', FS);
    bb=get(gca);
    set(bb.XLabel, 'FontSize', FS);
    set(bb.YLabel, 'FontSize', FS);
    set(bb.ZLabel, 'FontSize', FS);
    set(bb.Title, 'FontSize', FS);
    print(gcf, '-depsc', 'tst0_s2');
    
end