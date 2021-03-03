%%%% last copied: 01/29 
%%%% this file works on windows for process sth.
%% not from json but draw myself for test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start & global setting
clc;
clear all;
global test_show_im
test_show_im=1;
global xpcnt;global ypcnt;
xpcnt=0.33;ypcnt=0.33;
global skip;global opc_width;
skip=10;     %��������
opc_width=5;
global minEPE;
minEPE=0.5;%%%%%%%%%%%%%%%%%%%%%%%%%%%?
%% original pic draw
w=200;h=200;
img_target=zeros(w,h);
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
%% cut center pic
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
OPC(img_target_cut,img_process_cut);

%%
% img_process only 0/1
bw_p=1-logical(img_process);      % �ڰ���ɫ
imwrite(bw_p,save_p_path,'png');

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
    for k = 1:length(B)
       boundary = B{k};     %��������ľ��ο�
       plot(boundary(:,2), boundary(:,1), color, 'LineWidth', 1)    %plot����(0,0)���½�
       skip=10;
       scatter(boundary(1:skip:length(boundary),2), boundary(1:skip:length(boundary),1))    %������
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function img=draw_rec(img,x1,y1,x3,y3,z)
% (i,j)�����i�е�j�У��ȼ���MATLAB����λ��(j,MAMY-i)���������Ͻ�Ϊ(0,0)��(j,i)
% ���ڱ����������������£��Ժ�����Ż� ###TODO
function img=draw_rec(img,y1,x1,y3,x3,z)
    % [a,b] = size(img);
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

function [img,bs]=cal_opc(img_source,bs)

end

function EPE=cal_EPE(img_source,img_process)

end

function [img_process,bs]=opc_process(img_process,bs,img_source,img_source_i,EPE)
    global minEPE
    % ��ֹ�������ﵽminEPE or ȫ�����������
    function flag=types_opc_process(type)
        [img_process,bs]=cal_opc(img_process,bs,type);

        img_process_o=restore_center_img(img_source,img_process);       
        img_filtered_o=filtering(img_process_o);        
        img_filtered=cut_center_img(img_filtered_o);

        EPE=cal_EPE(img_source_i,img_filtered);
        
        if( EPE>minEPE || sum(sum(cell2mat(bs(:,3))==0))==0 )  
            flag=true
            return
        else
            [img_process,bs]=opc_process(img_process,bs,img_source_o,img_source_i,EPE);
        end
        flag=false
    end
    
    for type = [0,1,-1]
        if types_opc_process(type)
            return
        end
    end
    % if all failed TODO
    
    % cal_opc: choose a line, draw new img, 0
    % cal_opc: choose a line, draw new img, 1
    % cal_opc: choose a line, draw new img, -1
end

function OPC(img_source)
    img_process_o=img_source;
    img_source_i=cut_center_img(img_source);
    img_process=img_source_i;   % --> img_opc not filtered
    % _o��׺��ԭ��СͼƬ��δ�и�
    % img_source_i ԭͼƬ�и����
    % img_source_o �и�ǰ���˲�ǰ��һֱ�ڱ仯
    
    
    % to get initial work, sample & stack
    global skip;
    BW = imbinarize(img_process);
    [B,L] = bwboundaries(BW);
    bs=cell(length(B),4);   % no.|| x  y  stack_to_record_sample index || k*4
    %cell2mat(p(1,1)) to get data
    for k = 1:length(B)
       boundary = B{k};     %��������ľ��ο�
       xs=boundary(1:skip:length(boundary),2);
       ys=boundary(1:skip:length(boundary),1);
       bs(k,1)={xs};
       bs(k,2)={ys};
       bs(k,3)={zeros(1,len(xs))};
       bs(k,4)={len(xs)};
    end
     
    img_filtered_o=filtering(img_process_o);  % filtering need no cut image    
    img_filtered=cut_center_img(img_filtered_o);
    
    EPE=cal_EPE(img_source_i,img_filtered);
    
    %%%% cal_min_EPE and its img_process TODO
    opc_process(img_process,bs,img_source,img_source_i,EPE);  % get: img_process EPE
    
%     global minEPE
%     % ��ֹ�������ﵽminEPE or ȫ�����������
%     while( EPE>minEPE || sum(sum(cell2mat(bs(:,3))==0))==0 )  
%         [img_process,bs]=cal_opc(img_process,bs);
%         
%         img_process_o=restore_center_img(img_source,img_process);       
%         img_filtered_o=filtering(img_process_o);        
%         img_filtered=cut_center_img(img_filtered_o);
% 
%         EPE=cal_EPE(img_source_i,img_filtered);
%     end
    
    % end
    img_process_o=restore_center_img(img_source,img_process);       
    img_filtered_o=filtering(img_process_o);        
    img_filtered=cut_center_img(img_filtered_o);
    
    % show result
    figure();
    subplot(1,2,1),imshow(img_source_i);
    show_edge(img_source_i,'r');    
    subplot(1,3,2),imshow(img_process);
    show_edge(img_process,'r');
    subplot(1,3,3),imshow(img_filtered,[]);
    show_edge(img_filtered,'r');
    EPE
    
end