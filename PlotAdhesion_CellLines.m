clear all
close all

cl = get(groot,'defaultAxesColorOrder');


FilePath{1} = 'Z:\Ben Yeoman\Parallel Plate\Variable Width Chamber\Version 2\MDAMB231\Data';
FilePath{2} = 'Z:\Ben Yeoman\Parallel Plate\Variable Width Chamber\Version 2\MDAMB468\Data'; 

FilePath{3} = 'Z:\Ben Yeoman\Parallel Plate\Variable Width Chamber\Version 2\BT549\Data';
FilePath{5} = 'Z:\Ben Yeoman\Parallel Plate\Variable Width Chamber\Version 2\BT20\Data';

FilePath{4} = 'Z:\Ben Yeoman\Parallel Plate\Variable Width Chamber\Version 2\SUM1315\Data';
FilePath{6} = 'Z:\Ben Yeoman\Parallel Plate\Variable Width Chamber\Version 2\MCF7\Data'; 

FilePath{7} = 'Z:\Ben Yeoman\Parallel Plate\Variable Width Chamber\Version 2\MCF10A\Data'; 
FilePath{8} = 'Z:\Ben Yeoman\Parallel Plate\Variable Width Chamber\Version 2\MCF10AT\Data';
FilePath{9} = 'Z:\Ben Yeoman\Parallel Plate\Variable Width Chamber\Version 2\MCF10A-DCIS\Data';

eqn = 'exp(-(x/m)^b)';
ft = fittype( eqn, 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0];
opts.StartPoint = [1 1];
opts.Upper = [Inf 800];

maxTau = 400*5/6;
nL = 61;
Len = 46.5;
load('Z:\Ben Yeoman\Matlab\ComsolShears')
t=t./max(t)*maxTau;
L=linspace(0,Len,nL+1);
t(l>max(L))=[];
l(l>max(L))=[];
Shear=t(ceil(linspace(1,numel(l),nL)));


for F=1:numel(FilePath) 
    strid = strfind(FilePath{F},'\');
    name{F} = FilePath{F}(strid(end-1)+1:strid(end)-1);
    clear tPre tPost 
    nc=1;
    a=dir([FilePath{F} '\Slide*']);
    matfile = natsortfiles({a.name});

    for i=1:numel(a)
            
        try
            Folder = FilePath{F};
            load([Folder '\' matfile{i}])

            tPre(:,nc)=Pre;
            tPost(:,nc)=Post;
            nc=nc+1;

            y1 = Post./Pre;
            [xData, yData] = prepareCurveData( shear(~isnan(y1)), y1(~isnan(y1)));
            [f1, g1] = fit( xData, yData, ft, opts );
            t50{F}(i) = f1.m*(-log(0.5))^(1/f1.b);
            GF{F}(i) = g1.rsquare;
            if max(shear)<300
                stp=1;
            end

        catch e
            t50{F}(i) = 0;
            GF{F}(i) = 0;
            disp(e.message)
        end
        

    end
    nPre{F} = tPre;
    nPost{F} = tPost;
    Fpre{:,F} = sum(tPre,2);
    Fpost{:,F} = sum(tPost,2);


       
    y = Fpost{:,F}./Fpre{:,F};
    x = Shear;
    [xData, yData] = prepareCurveData( x, y);
    [f, gof] = fit( xData, yData, ft, opts );
    ci = confint(f,0.95);

    pl{:,F} = f(Shear);
    shr{:,F} = Shear;
    tau50(F,1) = f.m*(-log(0.5))^(1/f.b);
    ci_t50(:,F) = exp(ci(:,2))';
    sig(F,1) = f.b;
    nCells(F,1) = sum(Fpre{:,F});
    R(F,1) = gof.adjrsquare;
    ci = confint(f,0.95);
    
end

%%
cl = get(groot,'defaultAxesColorOrder');
for i=1:numel(FilePath)
    LegNames{i} = name{i};
end

    
%%
PlotSpeed

%%
avgT50=1;
if avgT50 
    for i=1:numel(FilePath)
        if i==9
            idx = GF{i}>-1;
            T(i,1) = mean(t50{i}(idx));
            S(i,1) = std(t50{i}(idx))/sqrt(sum(idx));
        else
            idx = GF{i}>0.5;
            T(i,1) = mean(t50{i}(idx));
            S(i,1) = std(t50{i}(idx))/sqrt(sum(idx));
        end
        nSample(i,1) = sum(idx);
        lNames{i} = [LegNames{i} ' (n=' num2str(nSample(i,1)) ')'];

    end
    tau50 = T;
    xneg = S;
    xpos = S;
else
    xneg = tau50'-ci_t50(1,:);
    xpos = ci_t50(2,:)-tau50';
end

yneg = sVel';
ypos = sVel';

%%
figure
set(gcf,'position',[10 10 700 500])
hold on
mPot = [5 5 5 5 1 1 1 1 1];
shape = ['o','^','s','p','o','^','s','p','d'];
errorbar(tau50,mVel,yneg,ypos,xneg,xpos,'.','MarkerSize',10,'Color','k');
for i=1:numel(FilePath)
    sc(i) = scatter(tau50(i),mVel(i),125,'filled',shape(i));
    sc(i).CData = cl(mPot(i),:);
end
xlabel('t_5_0 (dyn/cm^2)')
ylabel('Cell Speed (um/hr)')
ylim([0 30])

x=1:600;
tbl = table(T,mVel);

% modelfun = @(b,x)b(1).*x.^b(2)+b(3);
% beta0  = [ 1.882e+04   -1.617 6.175]; 
% mdl = fitnlm(tbl,modelfun,beta0);
% y=modelfun(beta0,x);

mdl = fitlm(tbl);
y=mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2).*x;

plot(x,y,'k','LineStyle','--')

l=legend(sc(1:numel(lNames)),lNames,'Location','best');
set(gca,'fontsize',24,'FontName', 'Arial')
l.FontSize = 16;
l.FontName = 'Arial';
set(l,'color','none');
set(gca,'color','none')
box on



yneg = sDis';
ypos = sDis';

%%
figure
set(gcf,'position',[10 10 700 500])
hold on

errorbar(tau50,mDis,yneg,ypos,xneg,xpos,'.','MarkerSize',10,'Color','k');
for i=1:numel(lNames)
    sc(i) = scatter(tau50(i),mDis(i),125,'filled',shape(i));
    sc(i).CData = cl(mPot(i),:);
end

x=1:600;
tbl = table(T,mDis);

% modelfun = @(b,x)b(1).*x.^b(2)+b(3);
% beta0  = [ 438.9   -0.09948  -223.9];
% mdl = fitnlm(tbl,modelfun,beta0);
% y=modelfun(beta0,x);

mdl = fitlm(tbl);
y=mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2).*x;
plot(x,y,'k','LineStyle','--')

xlabel('t_5_0 (dyn/cm^2)')
ylabel('Displacement (um)')
ylim([0 100])
l=legend(sc(1:numel(lNames)),lNames,'Location','best');
set(gca,'fontsize',24,'FontName', 'Arial')
l.FontSize = 16;
l.FontName = 'Arial';
set(l,'color','none');
set(gca,'color','none')
box on




%%
box on
set(gca,'fontsize',18,'FontName', 'Arial')
set(gca,'color','none')


figure
set(gcf,'position',[250   100   800   600])
hold on
for i = 1:length(tau50)
    h=bar(i,tau50(i));
    set(h,'FaceColor',cl(mPot(i),:));
end
ylabel('t50 (dyn/cm2)')
errorbar(tau50,S,'LineStyle','none','Color','k')
xtickangle(45)
xticklabels(labels)
ylim([0 800])
set(gca,'fontsize',24,'FontName', 'Arial')
box on
