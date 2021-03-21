%%%% last copied: 01/29 
%%%% this file works on windows for process sth.
%% not from json but draw myself for test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
%% start & global setting
clc;
clear all;
%% config
global test_show_im
test_show_im=0;
global xpcnt;global ypcnt;
xpcnt=0.33;ypcnt=0.33;
global skip;global opc_width;
skip=30;     %��������
opc_width=5;
global minEPE;
minEPE=0.5;%%%%%%%%%%%%%%%%%%%%%%%%%%%?
%% original pic draw
w=200;h=200;
img_target=zeros(w,h);
% in draw sys -->x  |y
% draw rec L-1
w1=20;h1=50;
x1=70;y1=70;
x2=x1+w1;y2=y1+h1;
img_target=draw_rec(img_target,x1,y1,x2,y2,1);  %%%%%%%
% draw rec L-2
w2=50;h2=20;
x3=x2;y4=y2;
x4=x3+w2;y3=y4-h2;
x2=70;y2=60; %L���ߵ����Ͻ�
img_target=draw_rec(img_target,x3,y3,x4,y4,1);

img_target=cloned_part_img(img_target,140,20,img_target,x1,y1,40,40);
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
imshow(img_target_wb);
show_edge(img_target,'r');
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
    show_edge(img_target_cut,'r');
    subplot(1,2,2),imshow(img_process_cut,[]);
    show_edge(img_process_cut,'b');
end
%% OPC
% OPC(img_target_cut,img_process_cut);
%%%%%%%%%%%%%%%%%%%55test%%%%%%
% clc;
% close all;
global ttt
ttt=0;
global test_edge_img
test_edge_img=img_target_cut;
OPC(img_target);
% %%
% close all
% %%
% % img_process only 0/1
% bw_p=1-logical(img_process);      % �ڰ���ɫ
% imwrite(bw_p,save_p_path,'png');

%% show result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure('name','result','color','w'),
% show_edge(img_process,'b');
% % ԭͼ������
% show_edge(img_target,'r');
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
function show_edge(img,color)
    BW = imbinarize(img);
    % imshow(BW);
    [B,L] = bwboundaries(BW);
    hold on
    global skip
    for k = 1:length(B)
       boundary = B{k};     %��������ľ��ο�
       plot(boundary(:,2), boundary(:,1), color, 'LineWidth', 1)    %plot����(0,0)���½�
%        skip=10;
       scatter(boundary(1:skip:length(boundary),2), boundary(1:skip:length(boundary),1))    %������
    end
end

function x=check_in_range(x,min,max)
    if x<min ;x=min; end
    if x>max ;x=max; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function img=draw_rec(img,x1,y1,x3,y3,z)
% (i,j)�����i�е�j�У��ȼ���MATLAB����λ��(j,MAMY-i)���������Ͻ�Ϊ(0,0)��(j,i)
% ���ڱ�����������������or�������£��Ժ���Ը��Ż�, ��� ###TODO
% x,y in plot is diff from row,coloum (x=c, y=r) 
% ����ʱ������������Ϊx,��Ϊy�����꣬img=draw_rec(img,y1,x1,y3,x3,z)
% function img=draw_rec(img,y1,x1,y3,x3,z)
    [a,b] = size(img);

    x1=check_in_range(x1,1,a);
    x3=check_in_range(x3,1,a);
    y1=check_in_range(y1,1,b);
    y3=check_in_range(y3,1,b);
    
    if(x1>x3 && y1>y3)
        t=x1;x1=x3;x3=t;
        t=y1;y1=y3;y3=t;
    end

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
% ���ڲ��ܳ�������Χ���Ժ���Ż�ֻճ���ص����� ###TODO
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
    if (sum([s,all]~=[4,9])==0 || ...
        sum([s,all]~=[4,4])==0 || ...
        sum([s,all]~=[4,6])==0 || ...
        sum([s,all]~=[2,6])==0 || ...
        sum([s,all]~=[2,4])==0 )
        flag=true;
    else
        flag=false;
    end
    
end

function [img,bs]=cal_opc(img_source,bs,type)
%     1. arr & idx -->p1,p2
%     2. figure out to this line what is in/out
%     3. draw_rect to change img
%     4. save situation to cell bs

    img=img_source; % do nothing
%     type
    % now see that k=1
    k=1;
    idx=bs{k,4};
    % if finished
    if idx==0
        return
    end
    arr=bs{k,3};
    
    % idx in (1~length)
    x1=bs{k,1}(idx);
    y1=bs{k,2}(idx);
    if idx==1
        idx2=length(bs{k,3});   % if last one
    else
        idx2=idx-1;
    end
    x2=bs{k,1}(idx2);
    y2=bs{k,2}(idx2);
        
    on=1;off=0;
    global opc_width
    [a,b] = size(img);
%     2. figure out to this line what is in/out
%     3. draw_rect to change img
    % here x,y means row,col
    if(x1==x2)
        % ����
        xx=x1+opc_width-1; %%%%%%%%%%%%%%%%%%%%
        xx=check_in_range(xx,1,a);
        if(img_source(xx,floor((y1+y2)/2))==on)
            % down is on
            direc=1;
        else
            % up is on
            direc=0;
        end
        switch (type+direc*10)
            case 1  % type=1, direc=0 draw down on
                img=draw_rec(img,x1+opc_width,y1,x2,y2,on);
            case -1 % type=-1, direc=0 draw up off
                img=draw_rec(img,x1-opc_width,y1,x2,y2,off);
            case 11 % type=1, direc=1 draw up on
                img=draw_rec(img,x1-opc_width,y1,x2,y2,on);
            case 9  % type=-1, direc=1 draw down off
                img=draw_rec(img,x1+opc_width,y1,x2,y2,off);
        end
    elseif(y1==y2)
        % |
        yy=y1+opc_width-1; %%%%%%%%%%%%%%%%%%%%
        yy=check_in_range(yy,1,b);
        if(img_source(floor((x1+x2)/2),yy)==on)
            % right is on
            direc=1;
        else
            % left is on
            direc=0;
        end
        switch (type+direc*10)
            case 1  % type=1, direc=0 draw right on
                img=draw_rec(img,x1+opc_width,y1,x2,y2,on);
            case -1 % type=-1, direc=0 draw left off
                img=draw_rec(img,x1-opc_width,y1,x2,y2,off);
            case 11 % type=1, direc=1 draw left on
                img=draw_rec(img,x1-opc_width,y1,x2,y2,on);
            case 9  % type=-1, direc=1 draw right off
                img=draw_rec(img,x1+opc_width,y1,x2,y2,off);
        end
    elseif (type~=0)
        % \ / 
        if type==1
%             img=draw_rec(img,x1,y1,x2,y2,on);
            % �Զ���Ϊ���Ļ����飿%%%%%%%%%%%%%%%%%%%%
            % �������⣬�жϵ�պ��Ƕ��㣿
            d=floor(opc_width);
            
            % one off, liek L left bottom is on
            % both two on, liek L right top is on
            if(img_source(x2,y1)==off | judge_corner(img,x1,y2))
                % center point is x1,y2
                img=draw_rec(img,x1-d,y2-d,x1+d,y2+d,on);
            elseif(img_source(x1,y2)==off | judge_corner(img,x2,y1))
                % center point is x2,y1
                img=draw_rec(img,x2-d,y1-d,x2+d,y1+d,on);
            else
                disp("both on, & both not corner, cant't judge.")%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
            end
            
        elseif type==-1
            img=draw_rec(img,x1,y1,x2,y2,off);
        end
        
    end

    bs(k,4)={idx-1};  
    arr(idx)=type;
    bs(k,3)={arr};
    

    global test_show_im
    global test_edge_img
    if test_show_im==10
        figure()
        imshow(img,[]);
%         show_edge(img,'r');
        show_edge(test_edge_img,'r');
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
%     EPE=0.3; % <minEPE
%     EPE=0.6; % >minEPE
%     EPE=rand(1)

%     global ttt
%     if ttt<35
%         ttt=ttt+1;
%         EPE=1;
%     else
%         EPE=0.2;
%     end
%     return


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
       boundary = B{k};     %��������ľ��ο�
       % here x,y means row,col
       xs=boundary(1:skip:length(boundary),1);
       ys=boundary(1:skip:length(boundary),2);
       bs(k,1)={xs};
       bs(k,2)={ys};
    end
    for k = 1:length(Bp)
       boundary = Bp{k};     %��������ľ��ο�
       % here x,y means row,col
       xs=boundary(:,1);
       ys=boundary(:,2);
       bsp(k,1)={xs};
       bsp(k,2)={ys};
    end
    
    EPE=0;
    not_find_EPE=10; %%%%%%%%%%%%%%%%%%%
    [a,b]=size(bs);
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
            
            if idx==-1
                vmin=not_find_EPE;%%%%%%%%%%;NOT FIND
            end
            
            EPE=EPE+vmin
        end 
    end    

end

function [img_process,bs,flag]=opc_process(bs,img_source,img_process_base)
    flag=false;   %
    global minEPE
    % now see that k=1
    k=1;
    
    img_process=img_process_base;
    
    % ��ֹ�������ﵽminEPE or ȫ�����������
    function flag=types_opc_process(type)
%         % now see that k=1
%         k=1;
        idx=bs{k,4};     % next 
%         bs{k,3}         % now
        % sum(sum(bs{k,3}==0))==0 
        flag=false;
        if( idx==1 )
            %   Ŀǰimg_process_base����������base��ʵ�ʻ�������img_source
            [img_process,bs]=cal_opc(img_process_base,bs,type);  
            
            img_process_o=restore_center_img(img_source,img_process);       
            img_filtered_o=filtering(img_process_o);        
            img_filtered=cut_center_img(img_filtered_o);
            img_source_i=cut_center_img(img_source);

            EPE=cal_EPE(img_source_i,img_filtered); 
            if( EPE<minEPE )  
                flag=true;
            end
        else
            [img_process,bs]=cal_opc(img_process_base,bs,type);            
            [img_process,bs,flag]=opc_process(bs,img_source,img_process);
        end
    end
    
    for type = [0,1,-1]
        if types_opc_process(type)
            flag=true;
            return
        end
        bs(k,4)={bs{k,4}+1};
    end
    
    % if all failed TODO
    
    % cal_opc: choose a line, draw new img, 0
    % cal_opc: choose a line, draw new img, 1
    % cal_opc: choose a line, draw new img, -1
end

function OPC(img_source)
    img_process_o=img_source;
    img_source_i=cut_center_img(img_source);
%     img_process=img_source_i;   % --> img_opc not filtered
    % _o��׺��ԭ��СͼƬ��δ�и�
    % img_source_i ԭͼƬ�и����
    % img_source_o �и�ǰ���˲�ǰ��һֱ�ڱ仯
    
    
    % to get initial work, sample & stack
    global skip;
    BW = imbinarize(img_source_i);
    [B,L] = bwboundaries(BW);
    bs=cell(length(B),4);   % no.|| x  y  stack_to_record_sample index || k*4
    %cell2mat(p(1,1)) to get data /// or bs{k,3}
    % now see that k=1 ---- only has 1 boundary
    for k = 1:length(B)
        boundary = B{k};     %��������ľ��ο�
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
        bs(k,4)={length(xs)};
       
    end

%     img_filtered_o=filtering(img_process_o);  % filtering need no cut image    
%     img_filtered=cut_center_img(img_filtered_o);
%     
%     EPE=cal_EPE(img_source_i,img_filtered);
    
    %%%% cal_min_EPE and its img_process TODO
    [img_process,bs,flag]=opc_process(bs,img_source,img_source_i);  % get: img_process EPE
    
    % end
    img_process_o=restore_center_img(img_source,img_process);       
    img_filtered_o=filtering(img_process_o);        
    img_filtered=cut_center_img(img_filtered_o);
    
    EPE=cal_EPE(img_source_i,img_filtered);
    
    % show result
    figure();
    subplot(1,3,1),imshow(img_source_i);
    show_edge(img_source_i,'r');    
    subplot(1,3,2),imshow(img_process);
    show_edge(img_process,'r');
    subplot(1,3,3),imshow(img_filtered,[]);
    show_edge(img_filtered,'r');
    EPE
    
end