clear all
load('sheardata.mat')


biop = 2;   %Biopsy, 1 for tumor, 2 for stroma, 3 for single fit tumor, 4 for single fit stroma
mThres = 2; %Number of mets for high metastatic potential


%Predictor variables
pred = [t50{biop},cFrac{biop}];

%Define classes for low or high metastatic potential
resp = mets{biop}>=mThres;       %Mets > threshold
for i=1:numel(resp)
    if resp(i)
        mpot(i,1) = {'HiMet'};
    else
        mpot(i,1) = {'LoMet'};
    end
end

%% Fit a logistic regression model
if size(pred,2)==2

    mdl = fitglm(pred,resp,'interactions','Distribution','binomial','Link','logit');
    b=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2) mdl.Coefficients.Estimate(3)...
        mdl.Coefficients.Estimate(4) ];
    xb=b(1)+b(2).*pred(:,1)+b(3).*pred(:,2)+b(4).*pred(:,1).*pred(:,2);

    %probability estimates from the logistic regression for low of high met potential
    scores = exp(xb)./(1 + exp(xb));  

%     mdl = fitglm(pred,mets{biop},'linear','Distribution','poisson','Link','log');
%     b=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2) mdl.Coefficients.Estimate(3)];...
%          %mdl.Coefficients.Estimate(4)];
%     xb=b(1)+b(2).*pred(:,1)+b(3).*pred(:,2)+b(4).*pred(:,1).*pred(:,2);
%     scores = exp(xb);  


else
    mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit');
    scores = mdl.Fitted.Probability;
end

[X,Y,T,AUC] = perfcurve(mpot,scores,'HiMet');




set(gcf,'position',[300 10 525 500])
plot(X,Y,'b','LineWidth',2)
hold on
scatter(X,Y,'b','filled')
% xlabel('False positive rate') 
% ylabel('True positive rate')
ylim([-.05 1.05])
xlim([-.05 1.05])
xticks([0 0.2 0.4 0.6 0.8 1])

%title('Adhesion + Cancer Fraction')
set(gca,'fontsize',24,'FontName', 'Arial')
set(gca,'color','none')
box on
grid off
%text(.5,.2,sprintf('AUC = %0.2f',AUC),'FontSize',24)


if size(pred,2)==2
%% Plot logistic regression probability estimates
sym = ['o','d','p'];
figure
set(gcf,'position',[300 10 800 750])
tau50 = linspace(0,350,20);
cF = linspace(0,1,20);
for i=1:20
    for j=1:20
        xb=b(1)+b(2)*tau50(i)+b(3)*cF(j)+b(4).*tau50(i).*cF(j);
        P(j,i) = exp(xb)./(1 + exp(xb)); 
    end
end

surf(tau50,cF,P)
hold on
c=colorbar;
c.Limits = [0,1];
cl = get(groot,'defaultAxesColorOrder');
% for i=1:numel(pred(:,1))
%    s=scatter3(pred(i,1), pred(i,2),100,250,'filled');
%    s.CData = cl(resp(i)+1,:);
%    s.Marker = sym(ad_idx{biop}(i));
% end


xlabel( 't50 (dyn/cm2)');
ylabel( 'Cancer Fraction');
view( 0, 90);
xlim([0 350])

set(gca,'fontsize',32,'FontName', 'Arial')
set(gca,'color','none')
box on
grid off

%%
figure
set(gcf,'position',[300 10 1000 400])
sym = ['o','d','p'];
subplot(1,2,1)
hold on
for i=1:numel(resp)
    s=scatter(pred(i,1),mets{biop}(i),125,sym(ad_idx{biop}(i)),'filled');
    s.CData = cl(resp(i)+1,:);
end
x=1:350;
tbl = table(pred(:,1),mets{biop});
mdl1 = fitlm(tbl);
y=mdl1.Coefficients.Estimate(1)+mdl1.Coefficients.Estimate(2).*x;
plot(x,y,'k','LineStyle','--')
ylim([0 90])
xlim([0 350])

title('Stroma Biopsy')
ylabel( 'GFP+ Lung Nodules');
xlabel('t50 (dyn/cm2)')
box on
set(gca,'fontsize',20,'FontName', 'Arial')
set(gca,'color','none')

subplot(1,2,2)
hold on
for i=1:numel(resp)
    s=scatter(pred(i,2),mets{biop}(i),125,sym(ad_idx{biop}(i)),'filled');
    s.CData = cl(resp(i)+1,:);
end
x=0:0.01:1;
tbl = table(pred(:,2),mets{biop});
mdl2 = fitlm(tbl);
y=mdl2.Coefficients.Estimate(1)+mdl2.Coefficients.Estimate(2).*x;
plot(x,y,'k','LineStyle','--')
ylim([0 90])

ylabel( 'GFP+ Lung Nodules');
xlabel('Cancer Fraction')
box on
set(gca,'fontsize',20,'FontName', 'Arial')
set(gca,'color','none')
xlim([0 1])
xticks([0 0.25 0.5 0.75 1.0])
set(gca,'color','none')
end