close all
clearvars

folder = 'Z:\Ben Yeoman\Parallel Plate\Variable Width Chamber\Version 2\CCD Experiment\T47D9\FN2_CD80\Plate2';
slide = 4;
maxTau = 400*5/6;
fitType = 2;

imPre = readnd2([folder '\' 'Pre.nd2']);
imPost = readnd2([folder '\' 'Post.nd2']);

Dpre = imPre(:,:,slide,1);
Dpost = imPost(:,:,slide,1);
Bpre = imPre(:,:,slide,2);
Bpost = imPost(:,:,slide,2);

b1 = Bpre;
b2 = Bpost;
[c1,~] = ProcessImage(Dpre);
[d1,~] = ProcessImage(Dpost);

res = size(Dpre,2)/3;
ylen = length(Dpre);


%Get Shear Values from Comsol
nL = floor(ylen/res);
Len = nL*res*1.6/1000;
load('Z:\Ben Yeoman\Matlab\ComsolShears')
t=t./max(t)*maxTau;
L=linspace(0,Len,nL+1);
t(l>max(L))=[];
l(l>max(L))=[];

try
    load([folder '\Slide_' num2str(slide) '.mat']);
    dataExists = 1;
    
catch
    dataExists = 0;

    c2 = imbinarize(c1);
    d2 = imbinarize(d1);
    r1 = imdilate(c2,strel('disk',20,0));
    c3 = bwpropfilt(c2,'Area',[30 800]);
    d3 = bwpropfilt(d2,'Area',[30 800]);
    d3(~r1) = 0;

    %remove cells on gasket
    pos_triangle = [1  24900 1 29158 325 29158];
    G = insertShape(zeros(size(b1)),'FilledPolygon',{pos_triangle});
    G = G(:,:,1);
    G(G>0)=1;
    d3(G==1) = 0;
    
    %Get post shear positions
    
    cen = vertcat(regionprops(d3,'Centroid').Centroid);
    x2=cen(:,1);
    y2=cen(:,2);

    r2 = insertShape(zeros(size(c3)),'FilledCircle',[x2 y2 40*ones(size(x2))]);
    r2 = r2(:,:,1);
    r2(r2>0)=1;
    c3(r2==1) = 0;
    c3(G==1) = 0;
    
    %Get pre shear positions
    cen = vertcat(regionprops(c3,'Centroid').Centroid);
    x1=cen(:,1);
    y1=cen(:,2);
    
    
    %Get cell counts
    Ypre=1.6*[y1;y2]/1000;
    Ypost=1.6*y2/1000;
    Pre = sum(Ypre<=L(2:end)&Ypre>L(1:end-1))';
    Post = sum(Ypost<=L(2:end)&Ypost>L(1:end-1))';  

    loc1 = 0;
    loc2 = 0;
end
shear=t(ceil(linspace(1,numel(l),nL)));

%%
f2 = figure(2);
h2 = imshow(brighten(rescale(b1),0.2));
hSP2 = imscrollpanel(f2,h2);
hold on
s2 = gobjects;
set(gcf,'position',[1 41 1920 963])
api2 = iptgetapi(hSP2);
api2.setVisibleLocation(0.5,loc1)

f1 = figure(1);
h1 = imshow(brighten(rescale(b2),0.2));
hSP1 = imscrollpanel(f1,h1);
hold on
s1 = gobjects;
set(gcf,'position',[1 41 1920 963])
api1 = iptgetapi(hSP1);
api1.setVisibleLocation(0.5,loc1)

mvin = 0;
while 1
    fig = gcf;
    currentFig = fig.Number; 
    
    %Get positions for Post image
    switch currentFig
        case 1
            r = api1.getVisibleImageRect();
            loc1=r(2);
        
            delete(s1(:))
            s1(1)=scatter(x2,y2,10,'filled','c');
            if mvin>numel(x2); mvin=0; end
            if mvin > 0 
                s1(2) = scatter(x2(mvin),y2(mvin),15,'m','filled');
            end
            
            [x2,y2,out,mvin] = UserInterface(x2,y2,mvin);
            Ypre=1.6*[y1;y2]/1000;
            Ypost=1.6*y2/1000;
            Pre = sum(Ypre<=L(2:end)&Ypre>L(1:end-1))';
            Post = sum(Ypost<=L(2:end)&Ypost>L(1:end-1))';
            save([folder '\Slide_' num2str(slide) '.mat'],'x1','x2','y1','y2','Pre','Post','shear','loc1','loc2')
            switch out
                case 'enter'
                    plotAdhesion(shear,Post./Pre,fitType);
                    figure(1)

                case 'up'
                    loc1=loc1-300;
                    api1.setVisibleLocation(0.5,loc1)
                    api2.setVisibleLocation(0.5,loc1)
        
                case 'down'
                    loc1=loc1+300;
                    api1.setVisibleLocation(0.5,loc1)
                    api2.setVisibleLocation(0.5,loc1)

                case 'right'
                    %Get rid of double counts
                    f2 = figure(2);
                    pPre = [x1,y1];
                    pPost = [x2,y2];
                    D = pdist2(pPost,pPre);
                    idx = min(D)<20;
                    x1(idx) = [];
                    y1(idx) = [];
                    mvin = 0;

                case 'space'
                    break
            end

        %Get positions for Pre image
        case 2

            r = api2.getVisibleImageRect();
            loc1=r(2);

            delete(s2(:))
            x = [x2;x1];
            y = [y2;y1];

            s2(2)=scatter(x2,y2,10,'filled','c');
            s2(1)=scatter(x1,y1,10,'filled','r');
            if mvin>numel(x); mvin=0; end
            if mvin > 0
                s2(3) = scatter(x(mvin),y(mvin),15,'m','filled');
            end

            X=x;
            Y=y;
            [x,y,out,mvout] = UserInterface(x,y,mvin);
            if numel(x)>numel([x2;x1])
                x1(end+1)=x(end);
                y1(end+1)=y(end);
            elseif numel(x)<numel([x2;x1])
                if mvin>numel(x2)
                    x1(mvin-numel(x2))=[];
                    y1(mvin-numel(y2))=[];
                else
                    x2(mvin) = [];
                    y2(mvin) = [];
                end
            else
                if mvin > 0
                    if mvin>numel(x2)
                        x1(mvin-numel(x2))=x(mvin);
                        y1(mvin-numel(y2))=y(mvin);
                    else
                        x2(mvin) = x(mvin);
                        y2(mvin) = y(mvin);
                    end
                end
            end
            mvin = mvout;

            Ypre=1.6*[y1;y2]/1000;
            Ypost=1.6*y2/1000;
            Pre = sum(Ypre<=L(2:end)&Ypre>L(1:end-1))';
            Post = sum(Ypost<=L(2:end)&Ypost>L(1:end-1))';
            save([folder '\Slide_' num2str(slide) '.mat'],'x1','x2','y1','y2','Pre','Post','shear','loc1','loc2')
            switch out
                case 'enter'
                    plotAdhesion(shear,Post./Pre,fitType);
                    figure(2)

                case 'up'
                    loc1=loc1-300;
                    api1.setVisibleLocation(0.5,loc1)
                    api2.setVisibleLocation(0.5,loc1)

                case 'down'
                    loc1=loc1+300;
                    api1.setVisibleLocation(0.5,loc1)
                    api2.setVisibleLocation(0.5,loc1)

                case 'left'
                    f1 = figure(1);
                    mvin = 0;

                case 'space'
                    break
                    
            end
    end
end
[f,g]=plotAdhesion(shear,Post./Pre,fitType);
if fitType==1
    disp([exp(f.m) f.p])
else
    disp(exp(f.m))
end




%%
function [x,y,out,mvout] = UserInterface(xin,yin,mvin)

out = '';
mvout = 0;


%Get Operating System
str = computer;
if strcmp(str,'PCWIN64')
    del = 127;
elseif strcmp(str,'MACI64')
    del = 8;
end

x = xin;
y = yin;
n = numel(x);

clear midx dmin
dmin = zeros;


%Setup Pointer
mouse = NaN(32);
mouse(16:17,[1:14,19:32]) = 2;
mouse([1:14,19:32],16:17) = 2;
set(gcf,'Pointer','custom')
set(gcf,'PointerShapeCData',mouse)
set(gcf,'PointerShapeHotSpot',[16 16])


%Gets Input
waitforbuttonpress;
crd = get(gca, 'CurrentPoint');
ytemp = crd(1,2);
xtemp = crd(1,1);
click = get(gcf,'SelectionType');
key = double(get(gcf,'CurrentCharacter'));

if strcmpi(click,'normal') && isempty(key)
   btn = 1;
elseif strcmpi(click,'open') && isempty(key)
   btn = 2;
elseif strcmpi(click,'extend') && isempty(key)
   btn = 2;
elseif strcmpi(click,'alt') && isempty(key)
   btn = 3;
else
   btn = key;
end
set(gcf,'CurrentCharacter',char(0))


if isempty(btn)
    btn = -1;
    if isempty(xtemp)
        btn = -2;
        xtemp = 0;
        ytemp = 0;
    end
end


if xtemp <= 478*3 && xtemp >= 0 && ytemp <= 29543 && ytemp >= 0
    
    switch btn
        
        %Select cell (Left click)
        case 1
            out = 'flip';
            dmin=pdist2([xtemp,ytemp],[x y]);
            mvout = find(dmin == min(dmin),1);
            if mvout == mvin
                mvout = 0;
            end
            
        %Count new cell (Shift + Left click)
        case 2
            out = 'add';
            n=n+1;
            x(n,1) = round(xtemp);
            y(n,1) = round(ytemp);
            
            %Repositions existing cell (Right click)
        case 3
            out = 'move';
            mvout = mvin;
            if mvin > 0
                x(mvin) = round(xtemp);
                y(mvin) = round(ytemp);
            end
            
            %Open closest cell to editing (e)
        case 101
            out = 'edit';
            dmin=pdist2([xtemp,ytemp],[x y]);
            mvout = find(dmin == min(dmin),1);
            if mvout == mvin
                mvout = 0;
            end
            
            %Delete
        case del
            if mvin > 0
                x(mvin) = [];
                y(mvin) = [];
            end
            dmin=pdist2([xtemp,ytemp],[x y]);
            mvout = find(dmin == min(dmin),1);
            if mvout == mvin
                mvout = 0;
            end
            out = 'delete';
            
            %Ctrl v
        case 22
            out = 'check';
            
            %o
        case 111
            out = 'overlay';
            
            %Ctrl z
        case 26
            out = 'undo';

            %Up arrow
        case 30
            out = 'up';

            %Down arrow
        case 31
            out = 'down';
            
        %Ctrl k
        case 11
            out = 'all';
            dmin=pdist2([xtemp,ytemp],[x y]);
            mvout = find(dmin == min(dmin),1);
            if mvout == mvin
                mvout = 0;
            end
            
            %Left arrow
        case 28
            out = 'left';
            mvout = mvin;
            
            %Right arrow
        case 29
            out = 'right';
            mvout = mvin;

            %p
        case 112
            out = 'plot';
            
            %Shift n
        case 78
            out = 'next';
            
            % Zoom in (Shift +)
        case 43
            out = 'zoomin';
            zoom on
            pause()
            zoom off
            
            
            % Zoom out (Shift -)
        case 95
            out = 'zoomout';
            zoom on
            zoom out
            
            
            %Space Bar
        case 32
            out = 'space';

            
            %Enter
        case 13
            out = 'enter';
    end
end
end

%% 
function [I3,I2] = ProcessImage(I)
    n=1;
    img = zeros(461,478);
    for j=2:63
        for i=0:2
            img(:,:,n) = I(461*j+1:461*j+461,478*i+1:478*i+478);
            n=n+1;
        end
    end

    M = uint16(repmat(mean(img,3),65,3));
    M(size(I,1)+1:end,:)=[];
    I2 = I-M;
    I3 = imadjust(I2,[0 0.2],[],1.5);
end


%% 
function [f,gof] = plotAdhesion(shear,y1,fitType)

if fitType == 1
    muEpi = log(866.24);     %ln(tau50) for epithelial cells
    sigEpi = 0.44656;   %Scale parameter for epithelial cell distribution
    a1 = num2str(muEpi);
    b1 = num2str(sigEpi);
    eqn = ['(1-p)*(1-(x.^(1/' b1 '))./(exp(' a1 './' b1 ')+x.^(1./' b1 ')))+p*(1-(x.^(1/b))./(exp(m./b)+x.^(1./b)))'];
    ft = fittype( eqn, 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0 0 0];
    opts.StartPoint = [0.6 log(50) 0.1];
    opts.Upper = [1 muEpi 1];

elseif fitType == 2
    eqn = '(1-(x.^(1/b))./(exp(m./b)+x.^(1./b)))';
    ft = fittype( eqn, 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [0 0];
    opts.StartPoint = [0.2 1];
    opts.Upper = [Inf Inf ];
end

    figure(3)
    clf
    hold on
    set(gcf,'position',[200 200 600 400])
    cl = get(groot,'defaultAxesColorOrder');
    y1(y1>1) = 1;
    scatter(shear,y1,10,cl(1,:),'filled')
    [xData, yData] = prepareCurveData( shear(~isnan(y1)), y1(~isnan(y1)));


    [f, gof] = fit( xData, yData, ft, opts );
    p1 = f(shear);

    fancyPlot({shear},{p1},{'xlabel','Shear Stess (dynes/cm^2)'},...
        {'ylabel','Survival Fraction'},{'lineStyle','-'},{'ylim',[0,1.2]},...
        {'spline'},{'lineWidth',2},{'fontSize',18},{'color',cl(1,:)});
    set(gca, 'XScale', 'log')
    
end