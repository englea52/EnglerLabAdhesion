
clear all
close all

maxTau = 400*5/6;
load('ComsolShears')
l=X(:,1);
t=X(:,2);
t=sort(t);


L=46.65;
l(X(:,1)>L)=[];
t(X(:,1)>L)=[];
t=t./max(t)*maxTau;

%Number of image rows and columns
n=1;
m = 61;

%Distance and shears at that distance
Tau=t(ceil(linspace(1,numel(l),m)));



%Values are from S survival curve
muEpi = (557);     %ln(tau50) for epithelial cells
sigEpi = 1/0.51;   %Scale parameter for epithelial cell distribution

%These values need to be measured experimentally for more cancer lines
%(close to MDA-MB231 survival curve)
muCan = log(103);        %Average ln(tau50) for cancer cells
sigCan = 0.42;      %Scale parameter for cancer cell distribution
sMu = log(70);   %Std Dev for average shear
sSig = 0;         %Std Dev for average scale parameter


 

cl = get(groot,'defaultAxesColorOrder');





%Average Cells/image
C = 80000/24*(29543*1434)*(1.6^2/10000^2)/m; %80k cells/cm2 x Area / n regions


c=1;
for h=1:2
    
    %Cancer cell percentage
    pCancer = (h-1)/1;
    
    %Generate random ln(tau50) and scale factor for cancer cells
    %mu = normrnd(muCan,sMu);
    mu = gamrnd(1.9652, 71.5722);   
    %sigma = normrnd(sigCan,sSig);
    sigma=gamrnd(1.5018, 2.0386);
    
    for i=1:m
%         r1 = random('LogLogistic',mu, sigma,poissrnd(pCancer*C*n),1);
%         r2 = random('LogLogistic',log(muEpi), 1/sigEpi,poissrnd((1-pCancer)*C*n),1);
        r1 = random('Weibull',mu, sigma,poissrnd(pCancer*C*n),1);
        r2 = random('Weibull',810.4, 1.4154,poissrnd((1-pCancer)*C*n),1);
        Pos3{i} = [r1;r2];
    end
    for i=1:m
        adh3(i) = sum(Pos3{i}>Tau(i))/numel(Pos3{i});
    end
    
    
    
    %Fitting predictions for percent cancer cells (p) and their tau50 (m)
    [xData, yData] = prepareCurveData( Tau, adh3' );
    a1 = num2str(muEpi);
    b1 = num2str(sigEpi);
    eqn = 'p*exp(-(x/m)^b)+(1-p)*exp(-(x/537.47)^0.7894)';
    ft = fittype( eqn, 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [1 0 0];
    opts.StartPoint = [2 100 0.5];
    opts.Upper = [15 537.47 1];
    
    try
        [f, gof] = fit( xData, yData, ft, opts );
        
        G(h,1) = gof.sse;
        G(h,2) = gof.rsquare;
        G(h,3) = gof.dfe;
        G(h,4) = gof.adjrsquare;
        G(h,5) = gof.rmse;
        

        
        %Percent errors
        err(h,1) = (f.p-pCancer)/pCancer*100;
        err(h,2) = (exp(f.m)-exp(mu))/exp(mu);
        err(h,3) = ((f.b-sigma)/sigma)*100;

        %Absolute errors
        aErr(h,1) = abs(f.p-pCancer)*100;
        %aErr(h,2) = abs(f.m*(-log(0.5))^(1/f.b)-exp(mu));
        aErr(h,2) = abs(f.m*(-log(0.5))^(1/f.b)-mu*(-log(0.5))^(1/sigma));
        aErr(h,3) = abs(f.b-sigma);
        
        %Predicted and simulated values
        P(h,:) = [f.p,pCancer];
        M(h,:) = [f.m*(-log(0.5))^(1/f.b),mu*(-log(0.5))^(1/sigma)];%[f.m,exp(mu)];
        S(h,:) = [f.b,sigma];
        p4 = f(Tau);
    catch e
        fprintf(e.message)
        err(h,1) = NaN;
        err(h,2) = NaN;
        err(h,3) = NaN;
        
        P(h,:) = [NaN,pCancer];
        M(h,:) = [NaN, exp(mu)];
        S(h,:) = [NaN, sigma];
        p4 = zeros(size(Tau));
    end
    
    
    if pCancer==0 || pCancer==0.5 || pCancer==1
        adh(:,c) = adh3;
        s=scatter(Tau, adh3,10,cl(c,:),'filled');
        p1(c) = fancyPlot({Tau},{p4},{'xlabel','Shear Stess (dynes/cm^2)'},...
            {'ylabel','Survival Fraction'},{'lineStyle','-'},{'ylim',[0,1.2]},...
            {'color',cl(c,:)},{'spline'},{'lineWidth',2},{'fontSize',18});
         set(gca, 'XScale', 'log');
        uistack(s,'top')
        c=c+1;
    end
    
end

l=legend(p1(1:c),{'0%%','50%%','100%%'},'Location','best');
l.EdgeColor = [1,1,1];
l.FontSize = 18;
set(l,'EdgeColor','none');
set(l,'color','none');

%%
ex =  G(:,4)<0.5 | M(:,1)>600;
idx = ~ex;
figure
hold on
c90m = quantile(abs((M(idx,1)-M(idx,2))),0.8);
line([0,max(M(idx,2))],[0,0],'color','k','lineWidth',1)
line([0,max(M(idx,2))],[c90m,c90m],'color','r','lineWidth',1)
line([0,max(M(idx,2))],[-c90m,-c90m],'color','r','lineWidth',1)
scatter(M(idx,2),M(idx,1)-M(idx,2),10,'filled')
xlabel('Simulated \tau_5_0 (dyn/cm^2)')
ylabel('Absolute Error (dyn/cm^2)')
ylim([-100 100])
xlim([0 max(M(idx,2))])
box on
set(gca,'fontsize',18,'FontName', 'Calibri')
set(gca,'color','none')
box on

% figure
% [~,I]=sort(B(:,2));
% hold on
% plot(0:0.1:max(B(I,2)),0:0.1:max(B(I,2)),'k','lineWidth',2)
% scatter((B(I,2)),(B(I,1)),10,'filled')
% xlabel('Simulated Scale Parameter (\sigma)')
% ylabel('Predicted Scale Parameter (\sigma)')
% set(gca,'fontsize',18,'FontName', 'Calibri')
% set(gca,'color','none')
% box on

%%
ex =  G(:,4)<0.5 | M(:,1)>600;
idx = ~ex;
figure
hold on
c90p = quantile(abs(((P(idx,1)-P(idx,2))*100)),0.8);
line([0,1],[0,0],'color','k','lineWidth',1)
line([0,1],[c90p,c90p],'color','r','lineWidth',1)
line([0,1],[-c90p,-c90p],'color','r','lineWidth',1)
scatter(P(idx,2),(P(idx,1)-P(idx,2))*100,10,'filled')
xlabel('Simulated Cancer Fraction')
ylabel('Absolute Error (Cancer Percentage)')
ylim([-50 50])
yticks([-50 -40 -30 -20 -10 0 10 20 30 40 50])
set(gca,'fontsize',18,'FontName', 'Calibri')
set(gca,'color','none')
box on

