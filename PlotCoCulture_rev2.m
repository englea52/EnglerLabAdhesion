clear all
close all



FilePath = 'Z:\Ben Yeoman\Parallel Plate\Variable Width Chamber\Version 2\Co-Culture\'; 
date = {'Nov32021','Dec72021','Dec152021','Dec172021','Dec212021'};
pCancer = {'0_MDAMB231','25_MDAMB231','50_MDAMB231','75_MDAMB231','100_MDAMB231'};
                    
% load('Z:\Ben Yeoman\Matlab\ComsolShears')
% L=linspace(0,46.6528,61+1);
% t=t./max(t)*400; 
% t(l>max(L))=[];
% l(l>max(L))=[];

for F=1:numel(pCancer)
    nt=1;
    
    
    avgD = [];
    avgS = [];
    for d=3:5
        clear tPre tPost Shr
        nb=1;
        for i=1:4
             try
                FP = [FilePath date{d} '\' pCancer{F}];
                load([FP '\Slide_' num2str(i) '.mat'])
                nPre(:,nt) = Pre;
                

                if F==1 || F==5
                    tumor = 0;
                else
                    tumor = 1;
                end
                [f,gof] = fitAdhesion(shear,Post./Pre,tumor);
                %T(nt,F) = exp(f.m);
                R2(nt,F) = gof.adjrsquare;
                
                
                if tumor
                    P(nt,F) = f.p;
                else
                    P(nt,F) = 0;
                end
                if R2(nt,F)>0.3 && P(nt,F)<.99
                    
                    T(nt,F) = f.m*(-log(0.5))^(1/f.b); %#ok<*SAGROW>
                    S(nt,F) = f.b;
                    tPre(:,nb)=Pre;
                    tPost(:,nb)=Post;
                    Shr(:,nb)=shear;
                    nb=nb+1;
                else
                    P(nt,F) = NaN;
                    T(nt,F) = NaN;
                    S(nt,F) = NaN;
                end
                nt=nt+1;
            catch e
                T(nt,F) = NaN;
                P(nt,F) = NaN;
                disp(e.message)
                nt=nt+1;
            end
            
        end 
        try
            avgD = [avgD; mean(tPost./tPre,2)]; %#ok<*AGROW> 
            avgS = [avgS; shear];
        catch
            avgD = [avgD; NaN(size(shear))];
            avgS = [avgS; NaN(size(shear))];
        end
    end
    nCells{F} = sum(nPre);
%     try
    Fdet(:,F) = avgD;
    Fshr(:,F) = avgS;

    y = Fdet(:,F);
    y(y>1)=1;
    x = Fshr(:,F);
    

         
    if F==1 || F==5
        tumor = 0;
    else
        tumor = 1;
    end

    [f,gof] = fitAdhesion(x,y,tumor);


    
    pl(:,F) = f(linspace(1,max(Fshr(:,F)),183));
    %tau50(F,1) = exp(f.m);
    tau50(F,1) = f.m*(-log(0.5))^(1/f.b);
    sig(F,1) = f.b;
    %nCells(F,1) = sum(Fpre(:,F));
    try
        pop(F,1) = nanmean(P(:,F));
    catch
        if F==1
            pop(F,1) = 0;
        else
            pop(F,1) = 1;
        end
    end
    R(F,1) = gof.adjrsquare;
%     catch
%         pl(:,F) = zeros(size(shear));
%     end

end

%%
Ncells=[mean(nCells{1}) mean(nCells{2}) mean(nCells{3}) mean(nCells{4}) mean(nCells{5})]';
Pcnt=[0 .25 .5 .75 1]';
aP = 100*Pcnt*Ncells(end)./Ncells;

[S1,idx] = sort(Fshr);
for j=1:5
    S2(:,j) = mean(reshape(S1(:,j),3,61))';
    D1(:,j) = Fdet(idx(:,j),j);
    D2(:,j) = mean(reshape(D1(:,j),3,61))';
    Derr(:,j) = std(reshape(D1(:,j),3,61))';
end

figure(1)
hold on
cl = get(groot,'defaultAxesColorOrder');
set(gcf,'position',[10 100 800 600])
for i=2:4%:numel(pCancer)-1
y = Fdet(:,i);
y(y>1)=1;
% s = scatter(Fshr(1:61,i),Fdet(1:61,i),25,cl(i,:));
% s2 = scatter(Fshr(62:61*2,i),Fdet(62:61*2,i),25,cl(i,:),'filled');
% s3 = scatter(Fshr(61*2+1:61*3,i),Fdet(61*2+1:61*3,i),25,cl(i,:),'filled','d');

s1 = scatter(S2(:,i),D2(:,i),25,cl(i,:),'filled');
errorbar(S2(:,i),D2(:,i),Derr(:,i),'Color',cl(i,:),'LineStyle','none')
Y1 = pl(:,i);
p(i) = fancyPlot({linspace(1,max(Fshr(:,F)),183)},{pl(:,i)},{'xlabel','Shear Stress (dynes/cm^2)'},...
{'ylabel','Survival Fraction'},{'lineStyle','-'},{'ylim',[0,1.2]},...
{'spline'},{'color',cl(i,:)},{'lineWidth',2},{'fontSize',18});
set(gca, 'XScale', 'log')
end
pop=pop*100;
legName = {sprintf('%0.1f%%  %0.1f%%',aP(1),pop(1)),sprintf('%0.1f%%  %0.1f%%',aP(2),pop(2)),...
    sprintf('%0.1f%%  %0.1f%%',aP(3),pop(3)),sprintf('%0.1f%%  %0.1f%%',aP(4),pop(4)),sprintf('%0.1f%%  %0.1f%%',aP(5),pop(5))};
l=legend(p(2:4),{legName{2},legName{3},legName{4}},'Location','southwest');
%l=legend(p(2:4),{legName{1},legName{5}},'Location','southwest');
l.EdgeColor = [1,1,1];
l.FontSize = 18;
xlim([0 333])
set(l,'EdgeColor','none');
set(l,'color','none');

%%
figure(2)
hold on
cl = get(groot,'defaultAxesColorOrder');
set(gcf,'position',[10 100 800 600])
for i=[1,5]%:numel(pCancer)-1
y = Fdet(:,i);
y(y>1)=1;
% s = scatter(Fshr(1:61,i),Fdet(1:61,i),25,cl(i,:));
% s2 = scatter(Fshr(62:61*2,i),Fdet(62:61*2,i),25,cl(i,:),'filled');
% s3 = scatter(Fshr(61*2+1:61*3,i),Fdet(61*2+1:61*3,i),25,cl(i,:),'filled','d');

s1 = scatter(S2(:,i),D2(:,i),25,cl(i,:),'filled');
errorbar(S2(:,i),D2(:,i),Derr(:,i),'Color',cl(i,:),'LineStyle','none')
Y1 = pl(:,i);
p(i) = fancyPlot({linspace(1,max(Fshr(:,F)),183)},{pl(:,i)},{'xlabel','Shear Stress (dynes/cm^2)'},...
{'ylabel','Survival Fraction'},{'lineStyle','-'},{'ylim',[0,1.2]},...
{'spline'},{'color',cl(i,:)},{'lineWidth',2},{'fontSize',18});
set(gca, 'XScale', 'log')
end
legName = {sprintf('%0.1f%%  %0.1f%%',aP(1),pop(1)),sprintf('%0.1f%%  %0.1f%%',aP(2),pop(2)),...
    sprintf('%0.1f%%  %0.1f%%',aP(3),pop(3)),sprintf('%0.1f%%  %0.1f%%',aP(4),pop(4)),sprintf('%0.1f%%  %0.1f%%',aP(5),pop(5))};

l=legend(p([1,5]),{legName{1},legName{5}},'Location','southwest');
l.EdgeColor = [1,1,1];
l.FontSize = 18;
xlim([0 333])
set(l,'EdgeColor','none');
set(l,'color','none');

%%
figure(3)
load('SimulatedShear')
hold on
cl = get(groot,'defaultAxesColorOrder');
set(gcf,'position',[10 100 800 600])
for i=[1,5]
    y = Fdet(:,i);
    y(y>1)=1;
    
    s1 = scatter(S2(:,i),D2(:,i),25,cl(i,:),'filled');
    errorbar(S2(:,i),D2(:,i),Derr(:,i),'Color',cl(i,:),'LineStyle','none')
    Y1 = pl(:,i);
    p(i) = fancyPlot({linspace(1,max(Fshr(:,F)),183)},{pl(:,i)},{'xlabel','Shear Stress (dynes/cm2)'},...
    {'ylabel','Survival Fraction'},{'lineStyle','-'},{'ylim',[0,1.2]},...
    {'spline'},{'color',cl(i,:)},{'lineWidth',2},{'fontSize',18});
    set(gca, 'XScale', 'log')

end

for i=1:2
    s1 = scatter(Tau,adh(:,i),25,cl(i+1,:),'filled');
    [f,gof] = fitAdhesion(Tau,adh(:,i),0);
    p(i+1) = fancyPlot({Tau},{f(Tau)},{'xlabel','Shear Stress (dynes/cm2)'},...
        {'ylabel','Survival Fraction'},{'lineStyle','-'},{'ylim',[0,1.2]},...
        {'spline'},{'color',cl(i+1,:)},{'lineWidth',2},{'fontSize',18});
        set(gca, 'XScale', 'log')
end


l=legend(p([1,5,3,2]),{'MDAMB231','MCF10A','Simulated Cancer','Simulated Epithelial'},'Location','southwest');
l.EdgeColor = [1,1,1];
l.FontSize = 18;
xlim([0 333])
set(l,'EdgeColor','none');
set(l,'color','none');


%%
function [f,g] = fitAdhesion(x,y,tumor)

    if tumor
        eqn = 'p*exp(-(x/m)^b)+(1-p)*exp(-(x/810.4)^1.4154)';
        ft = fittype( eqn, 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [0.5 0 0];
        opts.StartPoint = [2 100 0.5];
        opts.Upper = [15 810.4 1];
    else
        eqn = 'exp(-(x/m)^b)';
        ft = fittype( eqn, 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.Lower = [0.1 0];
        opts.StartPoint = [1 1];
        opts.Upper = [15 Inf];
    end

    [xData, yData] = prepareCurveData( x, y);
    [f, g] = fit( xData, yData, ft, opts );
end
