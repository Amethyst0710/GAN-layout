%%%% last copied: 01/29 
%%%% this file works on windows for process sth.
%% not from json but draw myself for test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% start & global setting
clc;
clear all;
close all;
%% config
global test_show_im
test_show_im=10;     
% 0 not showed, 1 show all filtered, 
% 2 show all opc process; 10 only show final result & original pic  
global xpcnt;global ypcnt;
% xpcnt=0.33;ypcnt=0.33;
xpcnt=0.5;ypcnt=0.5;
global skip;global opc_width;
skip=10;     %采样点间隔
opc_width=10;
global minEPE_rate;
minEPE_rate=0.1;%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for test ways
% if idx==2
%     input('k');
%     imshow(img_filtered,[]);
% end
% idx           
% if idx==2;disp(w);disp(EPE_w);end
% X = ['w=',num2str(type*w),',EPE_w=',num2str(EPE_w),',EPE_min_w=',num2str(EPE_min_w),',EPE_cmp',num2str(EPE_cmp)];
% disp(X)
%% original pic draw
w=200;h=200;
img_target=zeros(w,h);
% in draw sys -->x  |y
% draw rec -| 1-
w1=20;h1=50;
x1=70;y1=70;
x2=x1+w1;y2=y1+h1;
img_target=draw_rec(img_target,x1,y1,x2,y2,1);  %%%%%%%
% draw rec -| 2|
w2=50;h2=20;
x3=x2;y4=y2;
x4=x3+w2;y3=y4-h2;
x2=70;y2=60; %L横线的右上角
img_target=draw_rec(img_target,x3,y3,x4,y4,1);

img_target=draw_rec(img_target,105,90,120,110,1);  %%%%%%%

img_target=cloned_part_img(img_target,130,40,img_target,x1,y1,40,40);
img_target=cloned_part_img(img_target,20,140,img_target,x1,y1,40,60);
img_target=cloned_part_img(img_target,1,1,img_target,x1,y1,20,70);
img_target=cloned_part_img(img_target,160,150,img_target,x1,y1,40,40);
% % OPC ear
% r=0;
% img_target=ear4(img_target,x1,y1,x3,y3,r,1);

% img_target=draw_rec(img_target,1,1,20,20,1);
% % 靠近边缘环形测试
% xxx=5;yyy=95;www=10;
% img_target=draw_rec(img_target,xxx,xxx,xxx+www,yyy,1);         %%%%%%%%%%%%%%%
% img_target=draw_rec(img_target,xxx,xxx,yyy,xxx+www,1);
% img_target=draw_rec(img_target,yyy-www,xxx,yyy,yyy,1);        
% img_target=draw_rec(img_target,xxx,yyy-www,yyy,yyy,1);

% % L形凹进去
% r=3;
% img_target=draw_rec(img_target,x3-r,y2-r,x3+r,y2+r,0);

img_target_wb=1-img_target;        
% img是黑0白1，现取反，为看起来方便画图区域应该是黑
% ？取消这步之后没有边框问题了,
% 所以实际处理时应该白1表示画了的掩膜图形,即处理过程中显示的是黑底
if test_show_im==10
    figure();
    imshow(img_target_wb);
    show_edge(img_target,'r',true);
end

%% 低通滤波
img_process=filtering(img_target);
% figure,
% imshow(img_process,[]);
% cut center pic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img_target_cut=cut_center_img(img_target);
img_process_cut=cut_center_img(img_process);
if test_show_im==1
    figure()
    subplot(1,2,1),imshow(img_target_cut);
    show_edge(img_target_cut,'r',true);
    subplot(1,2,2),imshow(img_process_cut,[]);
    show_edge(img_process_cut,'b',true);
end
%% OPC
% OPC(img_target_cut,img_process_cut);
%%%%%%%%%%%%%%%%%%%55test%%%%%%
% clc;
% close all;
global test_edge_img
test_edge_img=img_target_cut;
OPC(img_target);
% %%
% close all
% %%
% % img_process only 0/1
% bw_p=1-logical(img_process);      % 黑白显色
% imwrite(bw_p,save_p_path,'png');

%% show result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('name','result','color','w'),
% show_edge(img_process,'b',true);
% % 原图形轮廓
% show_edge(img_target,'r',true);
% set(gca,'YDir','reverse');        %将x轴方向设置为反向(从上到下递增)。
% axis equal;
% axis off;

%%
%---------------function---------------%
% 从gds2ascii的结果文件中读取数据到cell中
function data_cell=read_data_from(fpath)
    % fpath='./json/2.json';
    fid=fopen(fpath,'r');
    
    %  如果读取不成功
    if fid==-1
        fprintf("ERROR: can't read from %s\n",fpath);
        data_cell={};
        return
    end
    
    %  如果读取成功
    xy_find=0;
    data_cell={};
    data_cnt=1;
    xy_array=[];
    xy_array_cnt=1;

    while feof(fid) ~= 1                %用于判断文件指针p在其所指的文件中的位置，如果到文件末，函数返回1，否则返回0  
        line = fgetl(fid);            %获取文档第一行    

        %    "XY",      [   100,    200,    300,    400     ]
        % 开始读取XY
        if contains(line,'"XY"')
            xy_find=1;
            line = fgetl(fid);   % 跳过该行
        end

        if xy_find==1 
            % 正式开始读取
            if contains(line,'[')
                line = fgetl(fid);   % 跳过该行
            end

            % 结束这次到XY读取
            if contains(line,']')
                xy_find=0;
                data_cell{data_cnt}=xy_array;
                data_cnt=data_cnt+1;
                xy_array=[];
                xy_array_cnt=1;
            else
                % 读入到cell
                x=str2double(line(isstrprop(line,'digit')));
                line = fgetl(fid);            %获取Y
                y=str2double(line(isstrprop(line,'digit')));

                xy_array(xy_array_cnt,:)=[x,y];
                xy_array_cnt=xy_array_cnt+1;
            end     % end of real start getting XY
        end     % end of start getting XY
    end     % end of while

    fclose(fid);
end         % end of function :  read_data_from

%轮廓提取
function p=show_edge(img,color,is_scatter)
    BW = imbinarize(img);
    % imshow(BW);
    [B,L] = bwboundaries(BW);
    hold on
    global skip
    for k = 1:length(B)
        boundary = B{k};     %包含最外的矩形框
        p=plot(boundary(:,2), boundary(:,1), color, 'LineWidth', 1);    %plot坐标(0,0)左下角
        if is_scatter
            scatter(boundary(1:skip:length(boundary),2), boundary(1:skip:length(boundary),1)) %采样点
        end
    end
end

function x=check_in_range(x,min,max)
    if x<min ;x=min; end
    if x>max ;x=max; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function img=draw_rec(img,x1,y1,x3,y3,z)
% (i,j)矩阵第i行第j列，等价于MATLAB坐标位置(j,MAMY-i)，设置左上角为(0,0)则(j,i)
% 输入坐标四种都可       
% x,y in plot is diff from row,coloum (x=c, y=r) 
% 调用时候若输入以右为x,下为y的坐标，img=draw_rec(img,y1,x1,y3,x3,z)
% 但实际画的时候会少画一个点以此避免重叠的问题，
% function img=draw_rec(img,y1,x1,y3,x3,z)
    [a,b] = size(img);

    x1=check_in_range(x1,1,a);
    x3=check_in_range(x3,1,a);
    y1=check_in_range(y1,1,b);
    y3=check_in_range(y3,1,b);
    
    
    if(x1>x3 && y1>y3)  % 右下左上
        t=x1;x1=x3;x3=t;
        t=y1;y1=y3;y3=t;
    elseif(x1<x3 && y1>y3)  % 右上左下
        t=y1;y1=y3;y3=t;
    elseif(x1>x3 && y1<y3)  % 左下右上
        t=x1;x1=x3;x3=t;
    end
    

    % 左上右下
    for i=x1:1:x3
        for j=y1:1:y3
            img(i,j)=z;
        end
    end
end

function img=ear4(img,x1,y1,x3,y3,r,z)
    img=draw_rec(img,x1-r,y1-r,x1+r,y1+r,z);
    img=draw_rec(img,x3-r,y1-r,x3+r,y1+r,z);
    img=draw_rec(img,x1-r,y3-r,x1+r,y3+r,z);
    img=draw_rec(img,x3-r,y3-r,x3+r,y3+r,z);
end

% clone part from img_source(xs,ys) to img_dest(xd,yd) 
% here x,y means row,col
% 现在不能超出矩阵范围，以后可优化只粘贴重叠部分 ###TODO-n
function img=cloned_part_img(img_dest,xd,yd,img_source,xs,ys,w,h)
    img=img_dest;
    img(xd:xd+w-1,yd:yd+h-1)=img_source(xs:xs+w-1,ys:ys+h-1);
end

function img=cut_center_img(img_source)
    global xpcnt;global ypcnt;
    [a,b]=size(img_source);
    w=floor(a*xpcnt);
    h=floor(b*ypcnt);
    x1=floor((a-w)/2);
    y1=floor((b-h)/2);
    img=img_source(x1:x1+w,y1:y1+h);
end

function img_source=restore_center_img(img_source,img_filled)
    global xpcnt;global ypcnt;
    [a,b]=size(img_source);
    w=floor(a*xpcnt);
    h=floor(b*ypcnt);
    x1=floor((a-w)/2);
    y1=floor((b-h)/2);
    img_source(x1:x1+w,y1:y1+h)=img_filled;
end

function img_process=filtering(img_target)
    % 低通滤波
    % 调用lybo.m中的"好像比较靠谱 但是这个没有转灰度，"
    global test_show_im
    if test_show_im==1
        figure,
        subplot(2,2,1),imshow(img_target);
        title('原图像');
    end

    % 将低频移动到图像的中心，这个也很重要
    s=fftshift(fft2(img_target));

    if test_show_im==1
        subplot(2,2,3),imshow(log(abs(s)),[]);% imagesc(abs(s));
        title('图像傅里叶变换取对数所得频谱');
    end

    % 求解变换后的图像的中心，我们后续的处理是依据图像上的点距离中心点的距离来进行处理
    [a,b] = size(img_target);
    a0 = round(a/2);
    b0 = round(b/2);
    d = min(a0,b0)/12;      %12  此处决定距离中心多远的频率不要  
    % d = 5;
    d = d^2;
    low_filter=zeros(a,b);
    for i=1:a
        for j=1:b
            distance = (i-a0)^2+(j-b0)^2;
            if distance<d
                low_filter(i,j) = s(i,j);
            end
        end
    end

    if test_show_im==1
        subplot(2,2,4),imshow(log(abs(low_filter)),[]);% imagesc(abs(low_filter));
        title('低通滤波频谱');
    end

    img_process = uint8(real(ifft2(ifftshift(low_filter))));

    if test_show_im==1
        subplot(2,2,2),imshow(img_process,[]);
        title('低通滤波后的图像');
    end
end

% here x,y means row,col
function flag=judge_corner(img,x,y)
    % judge by sum on the pix aroud it. 
    % standard L should be 4/9, may 4/6, 4/4.(2/6 2/4 ?
    [a,b] = size(img);
    xl=check_in_range(x-1,1,a);
    yl=check_in_range(y-1,1,b);
    xr=check_in_range(x+1,1,a);
    yr=check_in_range(y+1,1,b);
    all=(xr-xl+1)*(yr-yl+1);
    s=sum(sum(img(xl:xr,yl:yr)));
%     if (sum([s,all]~=[4,9])==0 || ...
%         sum([s,all]~=[4,4])==0 || ...
%         sum([s,all]~=[4,6])==0 || ...
%         sum([s,all]~=[2,6])==0 || ...
%         sum([s,all]~=[2,4])==0 )


    situations=[[4,9];[4,4];[4,6];[2,6];[2,4]];
    if ismember([s,all],situations,'rows')
%     if (isequal([s,all],[4,9])|| ...
%         isequal([s,all],[4,4])||...
%         isequal([s,all],[4,6])||...
%         isequal([s,all],[2,6])||...
%         isequal([s,all],[2,4]))
        flag=true;
    else
        flag=false;
    end
end

function [img,flag]=cal_opc_w(img,img_source_i,bs,type,k,w)
% do the opc according to the current idx
%     1. arr & idx -->p1,p2
%     2. figure out to this line what is in/out
%     3. draw_rect to change img
%     4. save situation to cell bs
    flag=true;
    % img is the base of cal_opc & return
    idx=bs{k,4};
    % if finished
    if idx==0
        return
    end
    
%     type
    if type==0; 
        return; 
    end
    
    % idx in (1~length)
    x1=bs{k,1}(idx);
    y1=bs{k,2}(idx);
    if idx==1
        idx2=length(bs{k,3});   % if last one
    else
        idx2=idx-1; % idx的顺时针方向的第一个（倒退一个
    end
    x2=bs{k,1}(idx2);
    y2=bs{k,2}(idx2);

    % here x,y means row,col
%     % if 2 point already no
%     if img(x1,y1)==off && img(x2,y2)==off
%         return
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%
    [img,flag]=cal_opc_w_process(img,img_source_i,x1,x2,y1,y2,type,w);
    
    global test_show_im
    global test_edge_img
    if test_show_im==2
        figure()
        imshow(img,[]);
%         show_edge(img,'r',true);
        show_edge(test_edge_img,'r',true);
    end
end

function [img,flag]=cal_opc_w_process(img,img_source_i,x1,x2,y1,y2,type,width)
    [a,b] = size(img);
    on=1;off=0;
    flag=true;
    %%%%%%%%%%TODO, 加一条，边缘不画
    if(x1==x2)
        % ——
        if(x1==a||x1==1);flag=false;return;end     % 边缘不画
        xx=x1+width-1; %%%%%%%%%%%%%%%%%%%%
        xx=check_in_range(xx,1,a);
        if(img_source_i(xx,floor((y1+y2)/2))==on)
            % down is on
            direc=1;
        else
            % up is on
            direc=0;
        end
        if(y1<y2);fit=1;else;fit=-1;end  % 调整重叠问题
        switch (type+direc*10)
            case 1  % type=1, direc=0 draw down on
                img=draw_rec(img,x1+width,y1+fit,x2,y2,on);
            case -1 % type=-1, direc=0 draw up off
                img=draw_rec(img,x1-width,y1+fit,x2,y2,off);
            case 11 % type=1, direc=1 draw up on
                img=draw_rec(img,x1-width,y1+fit,x2,y2,on);
            case 9  % type=-1, direc=1 draw down off
                img=draw_rec(img,x1+width,y1+fit,x2,y2,off);
        end
    elseif(y1==y2)
        % |
        if(y1==b||y1==1);flag=false;return;end     % 边缘不画
        yy=y1+width-1; %%%%%%%%%%%%%%%%%%%%
        yy=check_in_range(yy,1,b);
        if(img_source_i(floor((x1+x2)/2),yy)==on)
            % right is on
            direc=1;
        else
            % left is on
            direc=0;
        end
        if(x1<x2);fit=1;else;fit=-1;end  % 调整重叠问题
%         type+direc*10
        switch (type+direc*10)
            case 1  % type=1, direc=0 draw right on
                img=draw_rec(img,x1+fit,y1+width,x2,y2,on);
            case -1 % type=-1, direc=0 draw left off
                img=draw_rec(img,x1+fit,y1-width,x2,y2,off);
            case 11 % type=1, direc=1 draw left on
                img=draw_rec(img,x1+fit,y1-width,x2,y2,on);
            case 9  % type=-1, direc=1 draw right off
                img=draw_rec(img,x1+fit,y1+width,x2,y2,off);
        end
    elseif (type~=0)    % 理论0上不会传进来
        % \ / 
        if type==1
            f=on;fn=off;
        elseif type==-1
            f=off;fn=on;
        end
%         img=draw_rec(img,x1,y1,x2,y2,on);
        % 以顶点为中心画方块？%%%%%%%%%%%%%%%%%%%%
        % 存在问题，判断点刚好是顶点？
        d=floor(width);

        corner=[[1,b];[a,b];[a,1];[1,1]]; % 边缘不画
    
        % one off, liek L left bottom is on
        % both two on, liek L right top is on
        if(img_source_i(x2,y1)==fn || judge_corner(img_source_i,x1,y2))
            if ismember([x1,y2],corner,'rows');flag=false;return;end % 边缘不画
            % center point is x1,y2
            img=draw_rec(img,x1-d,y2-d,x1+d,y2+d,f);
        elseif(img_source_i(x1,y2)==fn || judge_corner(img_source_i,x2,y1))
            if ismember([x2,y1],corner,'rows');flag=false;return;end % 边缘不画
            % center point is x2,y1
            img=draw_rec(img,x2-d,y1-d,x2+d,y1+d,f);
        else
            disp("both on, & both not corner, cant't judge.");%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp({x1,y1;x2,y2});disp(width);disp(type);
            flag=false;
%                 judge_corner(img,x2,y1)
%             judge_corner(img,x1,y2)
            %%%%%%%%%%%%%% TODO 
        end

    end
end

function [kp,idx,vmin]=find_nearest_value(same,data,cell_list,find_type)
    vmin=1000;idx=-1;kp=-1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % idx not impossible except refer -1 wrong
    if find_type=='x'
        x=1;y=2;
    elseif find_type=='y'
        x=2;y=1;
    end  
    [a,b]=size(cell_list);
    for k = 1:a
        xlist=cell_list{k,x};
        ylist=cell_list{k,y};
        
%         same
        xx=xlist(ylist==same);% all the possible data
        if isempty(xx);continue;end
        [v,i] = min(abs(xx-data));
        
%         if find_type=='x'
%             xx=xlist(ylist==same);% all the possible data
%             [v,i] = min(abs(xx-data))
%         elseif find_type=='y'
%             yy=ylist(xlist==same);% all the possible data
%             [v,i] = min(abs(yy-data))
%         end                
        if v<vmin;vmin=v;idx=i;kp=k;end
    end
end

function EPE=cal_EPE(img_source,img_process)
% all EPE
% 1. get checkpoint in source & process boundaries
% 2. foreach cpoint find point have same x/y in boundaries
% 3. calculate EPE

    global skip;
    BW = imbinarize(img_source);
    [B,L] = bwboundaries(BW);
    BWp = imbinarize(img_process);
    [Bp,Lp] = bwboundaries(BWp);
    bs=cell(length(B),2);   % no.|| x  y || k*2
    bsp=cell(length(Bp),2);   % no.|| x  y || k*2
    %cell2mat(p(1,1)) to get data /// or bs{k,3}
    % now see that k=1 ---- only has 1 boundary
    for k = 1:length(B)
       boundary = B{k};     %包含最外的矩形框
       % here x,y means row,col
       xs=boundary(1:skip:length(boundary),1);
       ys=boundary(1:skip:length(boundary),2);
       bs(k,1)={xs};
       bs(k,2)={ys};
    end
    for k = 1:length(Bp)
       boundary = Bp{k};     %包含最外的矩形框
       % here x,y means row,col
       xs=boundary(:,1);
       ys=boundary(:,2);
       bsp(k,1)={xs};
       bsp(k,2)={ys};
    end
    
    EPE=0;
    not_find_EPE=10; %%%%%%%%%%%%%%%%%%%
    [a,b]=size(bs);
    global opc_width
    for k = 1:a
        % bs(k,1)={xs};bs(k,2)={ys};
        xs=bs{k,1};
        ys=bs{k,2};
        for i = 1:length(bs{k,1})
            if i==length(bs{k,1});j=1;else;j=i+1;end
            % x-same ----,find nearest y; y...
            % xy-diff, 
            if xs(i)==xs(j)
%                 [kp,idx,vmin]=find_nearest_value(xs(i),ys(i),bsp,'y');
                [kp,idx,vmin]=find_nearest_value(ys(i),xs(i),bsp,'x');
            elseif ys(i)==ys(j)
%                 [kp,idx,vmin]=find_nearest_value(ys(i),xs(i),bsp,'x');
                [kp,idx,vmin]=find_nearest_value(xs(i),ys(i),bsp,'y');
            else
                % 4 situation, 2 class
                if(xs(i)>xs(j) && ys(i)>ys(j)) || (xs(i)<xs(j) && ys(i)<ys(j))
                    % L or -|
                    [kp,idx,vmin]=find_nearest_value(ys(i),xs(i),bsp,'x');
                else
                    % |- or _|
                    [kp,idx,vmin]=find_nearest_value(xs(i),ys(i),bsp,'y');
                end
            end
            
            
            if idx==-1 || vmin>=opc_width   % not find
                vmin=not_find_EPE;%%%%%%%%%%;NOT FIND   TODO
            end
            
            EPE=EPE+vmin;
        end 
    end    

end

function [img_process,bs,flag,EPE]=opc_process(bs,img_source,img_process_base,EPE_min,knum)
    % 保证knum>=1时：
    [img_process,bs,flag,EPE]=opc_process_k(bs,img_source,img_process_base,EPE_min,1);
    for k = 2:knum
        [img_process,bs,flag,EPE]=opc_process_k(bs,img_source,img_process,EPE,k); %EPE_min
        if flag; return; end
    end
end

function [img_process,bs,flag,EPE_min]=opc_process_k(bs,img_source,img_process_base,EPE_min,k)
    flag=false;   %
    global minEPE_rate
    
    img_process=img_process_base;
    %%%%%%%%%%%%%% above seemed error
    idx=bs{k,4};     % next 
    if( idx==0 )
        return
    end
    
    global opc_width
    len=length(bs{k,3});  
    EPE_cmp=minEPE_rate*len*opc_width;
        
%     flag=false;
    choose_type=0;  % -1
    choose_w=1000;
    
    for type = [0, 1,-1]
        w=0;
        if(type==0)
            EPE=EPE_min;
        else
            img_source_i=cut_center_img(img_source);
            % calculate which width will minimum EPE, from 2 - opc_width
            EPE_min_w=EPE_min;  % opc_width:-1:2    2:opc_width
            %%%% 这里 w顺序不一样结果不一样是因为相同时取得是先得到最小值的
            for w = 2:opc_width
                [img_process,opc_flag]=cal_opc_w(img_process_base,img_source_i,bs,type,k,w);
%                 imshow(img_process,[]);
%                 qqq=input("input ");

                img_process_o=restore_center_img(img_source,img_process);       
                img_filtered_o=filtering(img_process_o);        
                img_filtered=cut_center_img(img_filtered_o);
                
%                 if idx==2
%                     input('k');
%                     imshow(img_filtered,[]);
%                 end

                EPE_w=cal_EPE(img_source_i,img_filtered);
%                 idx           
%                 if idx==2;disp(w);disp(EPE_w);end
%                 X = ['w=',num2str(type*w),',EPE_w=',num2str(EPE_w),',EPE_min_w=',num2str(EPE_min_w),',EPE_cmp',num2str(EPE_cmp)];
%                 disp(X)
                if( EPE_w<EPE_cmp  ) 
                    disp('EPE_w<EPE_cmp')
                    choose_w=w;EPE_min_w=EPE_w;flag=true;
                    break;
                elseif(EPE_w < EPE_min_w)
%                     disp('EPE_w < EPE_min_w')
                    choose_w=w;EPE_min_w=EPE_w;choose_type=type;
%                     if idx==52;disp(type);end
                end  
            end
            EPE=EPE_min_w;
            if flag; break;end
        end
        
%         type
%         EPE
        if( EPE<EPE_cmp  ) 
            EPE_min=EPE;flag=true;
%             EPE_min
            break;
        elseif(EPE < EPE_min)
%             X = ['EPE=',num2str(EPE),'<EPE_min=',num2str(EPE_min)];
%             disp(X)
            EPE_min=EPE;
%             EPE_min
%             disp("EPE < EPE_min");disp(type*w);disp(EPE_min);
        end  
    end    

    
    
    if flag
        bs(k,4)={idx-1};  
        arr=bs{k,3};
        if opc_flag==true
            arr(idx)=choose_type*choose_w;
        else
            arr(idx)=0;
        end
        bs(k,3)={arr};
    else
        % flag==false
        [img_process,opc_flag]=cal_opc_w(img_process_base,img_source_i,bs,choose_type,k,choose_w);
        bs(k,4)={idx-1};  
        arr=bs{k,3};
        if opc_flag==true
            arr(idx)=choose_type*choose_w;
        else
            arr(idx)=0;
        end
        bs(k,3)={arr};    
%         arr
        [img_process,bs,flag,EPE_min]=opc_process_k(bs,img_source,img_process,EPE_min,k);        
    end
    


end

function OPC(img_source)
    img_process_o=img_source;
    img_source_i=cut_center_img(img_source);
%     img_process=img_source_i;   % --> img_opc not filtered
    % _o后缀是原大小图片、未切割
    % img_source_i 原图片切割后结果
    % img_source_o 切割前，滤波前，一直在变化
    
    
    % to get initial work, sample & stack
    global skip;
    BW = imbinarize(img_source_i);
    [B,L] = bwboundaries(BW);  
    bs=cell(length(B),4);   % no.|| x  y  stack_to_record_sample index || k*4
    %cell2mat(p(1,1)) to get data /// or bs{k,3}
    for k = 1:length(B)
        boundary = B{k};     % 顺时针顺序
        % here x,y means row,col
        xs=boundary(1:skip:length(boundary),1);
        ys=boundary(1:skip:length(boundary),2);
        
%         a=[xs(1),ys(1);xs(end),ys(end)];
%         d=pdist(a,'euclidean')
%         if d<skip/2
%             ;
%         end
                
        bs(k,1)={xs};
        bs(k,2)={ys};
%         bs(k,1)={boundary(1:skip:length(boundary),:)} %%%%(y,x)
        bs(k,3)={zeros(1,length(xs))};
        bs(k,4)={length(xs)};   % 因为从最后开始，所以是逆时针顺序
       
    end

    img_filtered_o=filtering(img_process_o);  % filtering need no cut image    
    img_filtered=cut_center_img(img_filtered_o);
    EPE=cal_EPE(img_source_i,img_filtered);
    
    EPE_min=EPE;
    % cal_min_EPE and its img_process 
    [img_process,bs,flag,EPE]=opc_process(bs,img_source,img_source_i,EPE_min,length(B));  % get: img_process EPE
    
    if flag
        disp("find EPE="+int2str(EPE)+"< EPE_origin="+int2str(EPE_min));% find < EPE_min
    else
        disp("this is the minimum: "+int2str(EPE)+" & EPE_origin="+int2str(EPE_min));% this is the minest
    end
    
    % end
%     img_process_o=restore_center_img(img_source,img_process);       
%     img_filtered_o=filtering(img_process_o);        
%     img_filtered=cut_center_img(img_filtered_o);
%     
%     EPE=cal_EPE(img_source_i,img_filtered);
    
    global test_show_im
    if test_show_im==10
        % show result
        subp_num_r=2;
        subp_num_c=2;
        figure();
        % img source cut line & its filtered result
        subplot(subp_num_r,subp_num_c,1),imshow(img_source_i);
        show_edge(img_source_i,'r',true);   
%         img_filtered_so=filtering(img_source);  % filtering need no cut image    
%         img_filtered_s=cut_center_img(img_filtered_so);
%         show_edge(img_filtered_s,'g',false); 
        title("img source & its checkpoint");

        % img source & its processed result
        subplot(subp_num_r,subp_num_c,2),imshow(img_process);
        show_edge(img_process,'g',false);
        show_edge(img_source_i,'r',false);
        title("img source & its processed result");

    %     img_target_wb=1-img_filtered;% img是黑0白1，现取反，为看起来方便画图区域应该是黑
    %     subplot(1,3,3),imshow(img_target_wb,[]);
    %     show_edge(img_source_i,'r',false);

        % img source & its processed filtered result
        img_process_o=restore_center_img(img_source,img_process);       
        img_filtered_o=filtering(img_process_o);        
        img_filtered=cut_center_img(img_filtered_o);
        subplot(subp_num_r,subp_num_c,3),imshow(img_filtered,[]); % processed img filter
        p1=show_edge(img_filtered,'g',false);       % original img filtered  
        
        img_filtered_so=filtering(img_source);  % filtering need no cut image    
        img_filtered_s=cut_center_img(img_filtered_so);
        p2=show_edge(img_filtered_s,'b',false);         
        
        p3=show_edge(img_source_i,'r',false);  % original img
        legend([p1,p2,p3],{'process filtered','original filtered','original img'},'Location','best')
        title({['img source & its filtered result'],['& its processed filtered result']});   
        
        subplot(subp_num_r,subp_num_c,4);
        result_add = imadd(2*img_filtered_s,img_filtered);
%         imshow(result_add);
        imagesc(result_add);
        axis off;
        axis equal;
        title('2 filtered result');
        [a,b]=size(result_add);
        text(0,b+15,{'bluer-processed','yellower-original'})
    end

    
%     EPE;
    
end