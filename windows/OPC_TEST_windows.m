%%%% last copied: 01/29 
%%%% this file works on windows for process sth.
%% not from json but draw myself for test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% start & global setting
clc;
clear all;
close all;
%% config
global limit;
limit=5;
global test_show_im
test_show_im=10;     
% 0 not showed, 1 show all filtered, 
% 2 show all opc process; 10 only show final result & original pic  
global xpcnt;global ypcnt;
% xpcnt=0.33;ypcnt=0.33;
xpcnt=0.5;ypcnt=0.5;

% global skip;skip=5;     %��������
global opc_width;
opc_width=5;
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
x2=70;y2=60; %L���ߵ����Ͻ�
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
% % ������Ե���β���
% xxx=5;yyy=95;www=10;
% img_target=draw_rec(img_target,xxx,xxx,xxx+www,yyy,1);         %%%%%%%%%%%%%%%
% img_target=draw_rec(img_target,xxx,xxx,yyy,xxx+www,1);
% img_target=draw_rec(img_target,yyy-www,xxx,yyy,yyy,1);        
% img_target=draw_rec(img_target,xxx,yyy-www,yyy,yyy,1);

% % L�ΰ���ȥ
% r=3;
% img_target=draw_rec(img_target,x3-r,y2-r,x3+r,y2+r,0);

img_target_wb=1-img_target;        
% img�Ǻ�0��1����ȡ����Ϊ���������㻭ͼ����Ӧ���Ǻ�
% ��ȡ���ⲽ֮��û�б߿�������,
% ����ʵ�ʴ���ʱӦ�ð�1��ʾ���˵���Ĥͼ��,�������������ʾ���Ǻڵ�
if test_show_im==10
    figure();
    imshow(img_target_wb);
    show_edge(img_target,'r',true);
end

%% ��ͨ�˲�
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
img_process=OPC(img_target);
% %%
% close all
% %%
% % img_process only 0/1
% bw_p=1-logical(img_process);      % �ڰ���ɫ
% imwrite(bw_p,save_p_path,'png');

%% show result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('name','result','color','w'),
% show_edge(img_process,'b',true);
% % ԭͼ������
% show_edge(img_target,'r',true);
% set(gca,'YDir','reverse');        %��x�᷽������Ϊ����(���ϵ��µ���)��
% axis equal;
% axis off;

%%
%---------------function---------------%
% ��gds2ascii�Ľ���ļ��ж�ȡ���ݵ�cell��
function data_cell=read_data_from(fpath)
    % fpath='./json/2.json';
    fid=fopen(fpath,'r');
    
    %  �����ȡ���ɹ�
    if fid==-1
        fprintf("ERROR: can't read from %s\n",fpath);
        data_cell={};
        return
    end
    
    %  �����ȡ�ɹ�
    xy_find=0;
    data_cell={};
    data_cnt=1;
    xy_array=[];
    xy_array_cnt=1;

    while feof(fid) ~= 1                %�����ж��ļ�ָ��p������ָ���ļ��е�λ�ã�������ļ�ĩ����������1�����򷵻�0  
        line = fgetl(fid);            %��ȡ�ĵ���һ��    

        %    "XY",      [   100,    200,    300,    400     ]
        % ��ʼ��ȡXY
        if contains(line,'"XY"')
            xy_find=1;
            line = fgetl(fid);   % ��������
        end

        if xy_find==1 
            % ��ʽ��ʼ��ȡ
            if contains(line,'[')
                line = fgetl(fid);   % ��������
            end

            % ������ε�XY��ȡ
            if contains(line,']')
                xy_find=0;
                data_cell{data_cnt}=xy_array;
                data_cnt=data_cnt+1;
                xy_array=[];
                xy_array_cnt=1;
            else
                % ���뵽cell
                x=str2double(line(isstrprop(line,'digit')));
                line = fgetl(fid);            %��ȡY
                y=str2double(line(isstrprop(line,'digit')));

                xy_array(xy_array_cnt,:)=[x,y];
                xy_array_cnt=xy_array_cnt+1;
            end     % end of real start getting XY
        end     % end of start getting XY
    end     % end of while

    fclose(fid);
end         % end of function :  read_data_from

%������ȡ
function p=show_edge(img,color,is_scatter)
    BW = imbinarize(img);
    % imshow(BW);
    [B,L] = bwboundaries(BW);
    hold on
    global skip
    for k = 1:length(B)
        boundary = B{k};     %��������ľ��ο�
        p=plot(boundary(:,2), boundary(:,1), color, 'LineWidth', 1);    %plot����(0,0)���½�
        if is_scatter
            scatter(boundary(1:skip:length(boundary),2), boundary(1:skip:length(boundary),1)) %������
        end
    end
end

function show_seg_chk(seg_pnt,chk_pnt)
    for i=1:length(seg_pnt(:,1))
        scatter(seg_pnt{i,2}, seg_pnt{i,1},'x');
    end
    for i=1:length(chk_pnt(:,1))
        scatter(chk_pnt{i,2}, chk_pnt{i,1},'p');
    end
end

function x=check_in_range(x,min,max)
    if x<min ;x=min; end
    if x>max ;x=max; end
end

function b=check_in_range_arr(a,min,max)
    % b=a={1,2,3,4}
    b=cell(1,length(a));
    for i=1:length(a)
        b(i)={check_in_range(a{i},min,max)};
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function img=draw_rec(img,x1,y1,x3,y3,z)
% (i,j)�����i�е�j�У��ȼ���MATLAB����λ��(j,MAMY-i)���������Ͻ�Ϊ(0,0)��(j,i)
% �����������ֶ���       
% x,y in plot is diff from row,coloum (x=c, y=r) 
% ����ʱ������������Ϊx,��Ϊy�����꣬img=draw_rec(img,y1,x1,y3,x3,z)
% ��ʵ�ʻ���ʱ����ٻ�һ�����Դ˱����ص������⣬
% function img=draw_rec(img,y1,x1,y3,x3,z)
    [a,b] = size(img);

    x1=check_in_range(x1,1,a);
    x3=check_in_range(x3,1,a);
    y1=check_in_range(y1,1,b);
    y3=check_in_range(y3,1,b);
    
    
    if(x1>x3 && y1>y3)  % ��������
        t=x1;x1=x3;x3=t;
        t=y1;y1=y3;y3=t;
    elseif(x1<x3 && y1>y3)  % ��������
        t=y1;y1=y3;y3=t;
    elseif(x1>x3 && y1<y3)  % ��������
        t=x1;x1=x3;x3=t;
    end
    

    % ��������
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
% ���ڲ��ܳ�������Χ���Ժ���Ż�ֻճ���ص����� ###TODO-n
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
    % ��ͨ�˲�
    % ����lybo.m�е�"����ȽϿ��� �������û��ת�Ҷȣ�"
    global test_show_im
    if test_show_im==1
        figure,
        subplot(2,2,1),imshow(img_target);
        title('ԭͼ��');
    end

    % ����Ƶ�ƶ���ͼ������ģ����Ҳ����Ҫ
    s=fftshift(fft2(img_target));

    if test_show_im==1
        subplot(2,2,3),imshow(log(abs(s)),[]);% imagesc(abs(s));
        title('ͼ����Ҷ�任ȡ��������Ƶ��');
    end

    % ���任���ͼ������ģ����Ǻ����Ĵ���������ͼ���ϵĵ�������ĵ�ľ��������д���
    [a,b] = size(img_target);
    a0 = round(a/2);
    b0 = round(b/2);
    d = min(a0,b0)/12;      %12  �˴������������Ķ�Զ��Ƶ�ʲ�Ҫ  
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
        title('��ͨ�˲�Ƶ��');
    end

    img_process = uint8(real(ifft2(ifftshift(low_filter))));

    if test_show_im==1
        subplot(2,2,2),imshow(img_process,[]);
        title('��ͨ�˲����ͼ��');
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
        idx2=idx-1; % idx��˳ʱ�뷽��ĵ�һ��������һ��
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
    global limit;
    %%%%%%%%%%TODO, ��һ������Ե����
    if(x1==x2)
        % ����
        if(x1==a||x1==1);flag=false;return;end     % ��Ե����
        if(abs(y1-y2)<=limit);flag=false;return;end    % <limit����
        xx=x1+width-1; %%%%%%%%%%%%%%%%%%%%
        xx=check_in_range(xx,1,a);
        if(img_source_i(xx,floor((y1+y2)/2))==on)
            % down is on
            direc=1;
        else
            % up is on
            direc=0;
        end
        if(y1<y2);fit=1;else;fit=-1;end  % �����ص�����
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
        if(y1==b||y1==1);flag=false;return;end     % ��Ե����
        if(abs(x1-x2)<=limit);flag=false;return;end    % <limit����
        yy=y1+width-1; %%%%%%%%%%%%%%%%%%%%
        yy=check_in_range(yy,1,b);
        if(img_source_i(floor((x1+x2)/2),yy)==on)
            % right is on
            direc=1;
        else
            % left is on
            direc=0;
        end
        if(x1<x2);fit=1;else;fit=-1;end  % �����ص�����
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
    elseif (type~=0)    % ����0�ϲ��ᴫ����
        % \ / 
        if type==1
            f=on;fn=off;
        elseif type==-1
            f=off;fn=on;
        end
%         img=draw_rec(img,x1,y1,x2,y2,on);
        % �Զ���Ϊ���Ļ����飿%%%%%%%%%%%%%%%%%%%%
        % �������⣬�жϵ�պ��Ƕ��㣿
        d=floor(width);

        corner=[[1,b];[a,b];[a,1];[1,1]]; % ��Ե����

        tmp=check_in_range_arr({x1-1,x1+1,x2-1,x2+1},1,a);
        [x1n,x1p,x2n,x2p]=deal(tmp{:});
        tmp=check_in_range_arr({y1-1,y1+1,y2-1,y2+1},1,b);
        [y1n,y1p,y2n,y2p]=deal(tmp{:});
        
        % one off, liek L left bottom is on
        % both two on, liek L right top is on
        % 1 ���������淽����-1���������ж�---
%         if(img_source_i(x2,y1)==fn && type==1 || ...
%                 judge_corner(img_source_i,x1,y2)||...
%                 (type==-1 && img_source_i(x2n,y1n)==fn && img_source_i(x2p,y1p)==fn))
%             if ismember([x1,y2],corner,'rows');flag=false;return;end % ��Ե����
%                     
%             % center point is x1,y2
%             img=draw_rec(img,x1-d,y2-d,x1+d,y2+d,f);
%             
%         elseif(img_source_i(x1,y2)==fn && type==1 || ...
%                 judge_corner(img_source_i,x2,y1)||...
%                 (img_source_i(x1n,y2n)==fn && img_source_i(x1p,y2p)==fn && type==-1))
%             if ismember([x2,y1],corner,'rows');flag=false;return;end % ��Ե����
%             
%             % center point is x2,y1
%             img=draw_rec(img,x2-d,y1-d,x2+d,y1+d,f);
        %%%%%%%%%%%%%NEW%%%%%%%%%%%%
        if(img_source_i(x2,y1)==fn && type==1 || ...
                judge_corner(img_source_i,x1,y2)||...
                (type==-1 && ~judge_corner(img_source_i,x2,y1) && img_source_i(x2n,y1n)==fn && img_source_i(x2p,y1p)==fn))
            if ismember([x1,y2],corner,'rows');flag=false;return;end % ��Ե����
%             if ~judge_corner(img_source_i,x2,y1);flag=false;return;end % ��Ե����
                    
            % center point is x1,y2
            img=draw_rec(img,x1-d,y2-d,x1+d,y2+d,f);
            
        elseif(img_source_i(x1,y2)==fn && type==1 || ...
                judge_corner(img_source_i,x2,y1)||...
                (img_source_i(x1n,y2n)==fn && img_source_i(x1p,y2p)==fn && type==-1 && ~judge_corner(img_source_i,x1,y2)))
            if ismember([x2,y1],corner,'rows');flag=false;return;end % ��Ե����
%             if ~judge_corner(img_source_i,x1,y2);flag=false;return;end % ��Ե����
            
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

%     global skip;
    BW = imbinarize(img_source);
    [B,L] = bwboundaries(BW);
    BWp = imbinarize(img_process);
    [Bp,Lp] = bwboundaries(BWp);
    
    [~,check_pnt]=process_coordinate(B);
%     bs=cell(length(B),2);   % no.|| x  y || k*2
    bsp=cell(length(Bp),2);   % no.|| x  y || k*2
%     bs(:,1:2)=check_pnt;  % x y
    bs=check_pnt;
%     %cell2mat(p(1,1)) to get data /// or bs{k,3}
%     % now see that k=1 ---- only has 1 boundary
%     for k = 1:length(B)
%         boundary = B{k};     %��������ľ��ο�
%         % here x,y means row,col
%         xs=boundary(1:skip:length(boundary),1);
%         ys=boundary(1:skip:length(boundary),2);
%         bs(k,1)={xs};
%         bs(k,2)={ys};
%     end
%     %%%%%%%%%%%%%%%%%%% here for center in line test
%     bsc=cell(length(B),2);   % no.|| x  y || k*2
%     for k = 1:length(B)
%         xs=bs{k,1};
%         ys=bs{k,2};
%         xx=zeros(1,length(xs));
%         yy=zeros(1,length(xs));
%         for idx = 1:length(xs)
%             x1=xs(idx);
%             y1=ys(idx);
%             if idx==1
%                 idx2=length(xs);   
%             else
%                 idx2=idx-1;
%             end
%             x2=bs{k,1}(idx2);
%             y2=bs{k,2}(idx2);
%             if x1==x2
%                 x=x1;y=(y1+y2)/2;
%             elseif y1==y2
%                 x=(x1+x2)/2;y=y1;
%             else
%                 %%%%%%%%���������ô��???????
%                 x=x1;y=y1;
% %                 if (abs(x1-x2)>abs(y1-y2))
% %                     
% %                 else
% %                     
% %                 end         
%             end
%             xx(idx)=x;
%             yy(idx)=y;
%         end
%         bsc{k,1}=xx;
%         bsc{k,2}=yy;
%     end
%     bs=bsc;
%     %%%%%%%%%%%%%%%%%%%%%%
    for k = 1:length(Bp)
       boundary = Bp{k};     %��������ľ��ο�
       % here x,y means row,col
       xs=boundary(:,1);
       ys=boundary(:,2);
       bsp(k,1)={xs};
       bsp(k,2)={ys};
    end
    
    EPE=0;
    not_find_EPE=100; %%%%%%%%%%%%%%%%%%%
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
            
            
            if idx==-1 || vmin>=opc_width*1.3   % not find
                vmin=not_find_EPE;%%%%%%%%%%;NOT FIND   TODO
            end
            
            EPE=EPE+vmin;
        end 
    end    

end

function [img_process,bs,flag,EPE]=opc_process(bs,img_source,img_process_base,EPE_min,knum)
    % ��֤knum>=1ʱ��
    [img_process,bs,flag,EPE]=opc_process_k(bs,img_source,img_process_base,EPE_min,1);
    for k = 2:knum
        [img_process,bs,flag,EPE]=opc_process_k(bs,img_source,img_process,EPE,k); %EPE_min
        if flag; return; end
    end
end

function [img_process,bs,flag,EPE_min]=opc_process_k(bs,img_source,img_process_base,EPE_min,k)
    flag=false;   %
    opc_flag = true;
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
            %%%% ���� w˳��һ�������һ������Ϊ��ͬʱȡ�����ȵõ���Сֵ��
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

function [seg_pnt,chk_pnt]=process_coordinate(B)
    seg_pnt=cell(length(B),2); 
    chk_pnt=cell(length(B),2); 
    % B - every point
    % limit for seg_rule
    for k = 1:length(B)
        boundary = B{k};     % ˳ʱ��˳��
        % here x,y means row,col
        xs=boundary(:,1);
        ys=boundary(:,2);
        corner_k_x=zeros(length(xs),1);
        corner_k_y=zeros(length(xs),1);
        corner_k_x(1)=xs(1);
        corner_k_y(1)=ys(1);
%         check_pnt_k=[];
        cnt_k=1;
        is_row=true;
        
        % find the cornor
        for i = 2:length(xs)
            j=i-1;
            if xs(i)==xs(j)
                if ~is_row
                    % a cornor  |_ j i 
                    cnt_k=cnt_k+1;
                    corner_k_x(cnt_k)=xs(j);
                    corner_k_y(cnt_k)=ys(j);
                    is_row=true;
                end
            elseif ys(i)==ys(j)
                if is_row
                    % a cornor  _| j i 
                    cnt_k=cnt_k+1;
                    corner_k_x(cnt_k)=xs(j);
                    corner_k_y(cnt_k)=ys(j);
                    is_row=false;
                end
            else    % not consider it
                % / \ boundaries caused error
                if (abs(xs(i)-xs(j))==1 && abs(ys(i)-ys(j))==1)
                    if is_row
                        cnt_k=cnt_k+1;
                        corner_k_x(cnt_k)=xs(j);
                        corner_k_y(cnt_k)=ys(i);
                        is_row=false;
                        i=i+1;
                    else
                        cnt_k=cnt_k+1;
                        corner_k_x(cnt_k)=xs(i);
                        corner_k_y(cnt_k)=ys(j);
                        is_row=true;
                        i=i+1;
                    end
                end
                
                
            end            
        end
        corner_k_x=corner_k_x(1:cnt_k);
        corner_k_y=corner_k_y(1:cnt_k);
        
        % process coordinate
        [seg_pnt_k,chk_pnt_k]=process_coordinate_k(corner_k_x,corner_k_y);
        seg_pnt(k,1)={seg_pnt_k(:,1)};
        seg_pnt(k,2)={seg_pnt_k(:,2)};
        chk_pnt(k,1)={chk_pnt_k(:,1)};
        chk_pnt(k,2)={chk_pnt_k(:,2)};         
       
    end
end

function [seg,chk]=process_coordinate_k(xs,ys)
    cnt=length(xs);
    seg=zeros(cnt,2);
    chk=zeros(cnt,2);
    global limit;
    
    % if <3*limit ,ֱ�Ӽ�¼����
    % if >3*limit, seg seg����¼����-limit
    
    seg(1,:)=[xs(1),ys(1)];
    cnt_seg=1;
    last_vtx=true;  % last vtx recorded
    hammer_chk=false; %
    cnt_chk=0;
    first_seg=false; % if >3*limit & seg
    for i = 2:cnt+1
        j=i-1;
        %%%%%
%         if ys(j)==41 && xs(j)==56
%             disp("weew");
%         end
        %%%%
        if i==cnt+1;i=1;end
        if xs(i)==xs(j)
            x=xs;y=ys;
            
            % j i
            if abs(y(i)-y(j))>3*limit
                if i==2;first_seg=true;end
                % seg
                cnt_seg=cnt_seg+1;
                if y(j)>y(i);pn=-1;else;pn=1;end
                seg(cnt_seg,:)=[x(i),y(j)+pn*limit];% ����j������
                cnt_seg=cnt_seg+1;
                seg(cnt_seg,:)=[x(i),y(i)-pn*limit];% ����i������
                last_vtx=false;

                % chk
                cnt_chk=cnt_chk+1;
                chk(cnt_chk,:)=[x(i), floor((y(i)+y(j))/2)]; 
                if ~hammer_chk
                    % record 2 hammer_chk
                    cnt_chk=cnt_chk+1;
                    chk(cnt_chk,:)=[x(i), floor(y(j)+pn*limit/2)];% ����j������
                    hammer_chk=true;
                end
                if~(first_seg && j==cnt)
                    cnt_chk=cnt_chk+1;
                    chk(cnt_chk,:)=[x(i), floor(y(i)-pn*limit/2)];% ����i������  
                else
                    % remove first vertex
                    seg(1,:)=[];
                    cnt_seg=cnt_seg-1;
                end
                 
            else
                % seg
                if ~last_vtx
                     % remove last
                     seg(cnt_seg,:)=[x(j),y(j)];
                     last_vtx=true;
                else
                    % chk
                    cnt_chk=cnt_chk+1;
                 end
                 
                 if i~=1
                     cnt_seg=cnt_seg+1;
                     seg(cnt_seg,:)=[x(i),y(i)];
                 end

                 % chk
                 chk(cnt_chk,:)=[x(i), floor((y(i)+y(j))/2)];  
                 hammer_chk=false;
            end
            
        elseif ys(i)==ys(j)
            % same as xs(i)==xs(j)
            y=xs;x=ys;  % here diff
            
            % j i
            if abs(y(i)-y(j))>3*limit   %%%%%%%%%%%%%%%%%%%%%%TODO
                if i==2;first_seg=true;end
                % seg
                cnt_seg=cnt_seg+1;
                if y(j)>y(i);pn=-1;else;pn=1;end
                seg(cnt_seg,:)=[y(j)+pn*limit,x(i)];% ����j������
                cnt_seg=cnt_seg+1;
                seg(cnt_seg,:)=[y(i)-pn*limit,x(i)];% ����i������
                last_vtx=false;

                % chk
                cnt_chk=cnt_chk+1;
                chk(cnt_chk,:)=[floor((y(i)+y(j))/2),x(i) ]; 
                if ~hammer_chk
                    % record 2 hammer_chk
                    cnt_chk=cnt_chk+1;
                    chk(cnt_chk,:)=[ floor(y(j)+pn*limit/2),x(i)];% ����j������
                    hammer_chk=true;
                end
                if~(first_seg && j==cnt)
                    cnt_chk=cnt_chk+1;
                    chk(cnt_chk,:)=[floor(y(i)-pn*limit/2),x(i) ];% ����i������  
                else
                    % remove first vertex
                    seg(1,:)=[];
                    cnt_seg=cnt_seg-1;
                end
                 
            else
                % seg
                if ~last_vtx
                     % remove last
                     seg(cnt_seg,:)=[y(j),x(j)];
                     last_vtx=true;
                else
                    % chk
                    cnt_chk=cnt_chk+1;
                 end
                 
                 if i~=1
                     cnt_seg=cnt_seg+1;
                     seg(cnt_seg,:)=[y(i),x(i)];
                 end

                 % chk
                 chk(cnt_chk,:)=[floor((y(i)+y(j))/2),x(i)];  
                 hammer_chk=false;
            end
            
        else
         % pass
        end
        
    end
    seg=seg(1:cnt_seg,:);
    chk=chk(1:cnt_chk,:);
end

function img_process=OPC(img_source)
    img_process_o=img_source;
    img_source_i=cut_center_img(img_source);
%     img_process=img_source_i;   % --> img_opc not filtered
    % _o��׺��ԭ��СͼƬ��δ�и�
    % img_source_i ԭͼƬ�и����
    % img_source_o �и�ǰ���˲�ǰ��һֱ�ڱ仯
    
    
    % to get initial work, sample & stack
%     global skip;
    BW = imbinarize(img_source_i);
    [B,L] = bwboundaries(BW);  
    %%%%%%%%%%%%%%%%%%%%%5
    [seg_pnt,chk_pnt]=process_coordinate(B);
%     figure();
%     imshow(img_source_i);hold on;
%     scatter(seg_pnt{2,2}, seg_pnt{2,1},'x');
%     scatter(check_pnt{2,2}, check_pnt{2,1},'p');

    bs=cell(length(B),4);   % no.|| x  y  stack_to_record_sample index || k*4
    bs(:,1:2)=seg_pnt;  % x y
    for k = 1:length(B)
        xs=seg_pnt{k,1};
        bs(k,3)={zeros(1,length(xs))};
        bs(k,4)={length(xs)};   
    end
    %%%%%%%%%%%%%%%%%%%%%
%     bs=cell(length(B),4);   % no.|| x  y  stack_to_record_sample index || k*4
%     %cell2mat(p(1,1)) to get data /// or bs{k,3}
%     for k = 1:length(B)
%         boundary = B{k};     % ˳ʱ��˳��
%         % here x,y means row,col
%         xs=boundary(1:skip:length(boundary),1);
%         ys=boundary(1:skip:length(boundary),2);
%         
% %         a=[xs(1),ys(1);xs(end),ys(end)];
% %         d=pdist(a,'euclidean')
% %         if d<skip/2
% %             ;
% %         end
%                 
%         bs(k,1)={xs};
%         bs(k,2)={ys};
% %         bs(k,1)={boundary(1:skip:length(boundary),:)} %%%%(y,x)
%         bs(k,3)={zeros(1,length(xs))};
%         bs(k,4)={length(xs)};   % ��Ϊ�����ʼ����������ʱ��˳��
%        
%     end

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
        subplot(subp_num_r,subp_num_c,1),imshow(img_source_i);hold on;
        show_seg_chk(seg_pnt,chk_pnt);
%         show_edge(img_source_i,'r',true);   
% %         img_filtered_so=filtering(img_source);  % filtering need no cut image    
% %         img_filtered_s=cut_center_img(img_filtered_so);
% %         show_edge(img_filtered_s,'g',false); 
        title("img source & its checkpoint");

        % img source & its processed result
        subplot(subp_num_r,subp_num_c,2),imshow(img_process);
        show_edge(img_process,'g',false);
        show_edge(img_source_i,'r',false);
        title("img source & its processed result");

    %     img_target_wb=1-img_filtered;% img�Ǻ�0��1����ȡ����Ϊ���������㻭ͼ����Ӧ���Ǻ�
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