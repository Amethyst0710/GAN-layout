% input: 
% ascii_dir='ascii/demo3_162nm';
% fname='demo3_162nm';

picdir='./original_plot/';
resizedir='./original_data/';
p_datadir='./process_data/';
    
ascii_data_dir=strcat('./',ascii_dir,'/');
subdir  = dir( ascii_data_dir );
% 1, 2是list包含了. 和 .. 路径。
for i = 3 : length( subdir )
    subdirpath = fullfile( ascii_data_dir, subdir( i ).name);
    file = subdir( i ) ;              % 子文件
    filename=file.name;
    
    [startIndex, endIndex] = regexp(filename, '([0-9]+-){2}[0-9]+') ;
    split_arr=filename(startIndex:endIndex);    % 3300-5225-2000
    out = regexp(split_arr, '-', 'split');
    xx=str2double(cell2mat(out(1)));
    yy=str2double(cell2mat(out(2)));
    radius=str2double(cell2mat(out(3)));
    
  
    % eg. filename='clips_ringo_layer1-0_2700-6890-200.json';
    fpath=strcat(ascii_data_dir,filename);
    save_o_path=strcat(picdir,fname,'/o_',filename(1:end-5),'.jpg');
    save_r_path=strcat(resizedir,fname,'/r_',filename(1:end-5),'.jpg');
    save_p_path=strcat(p_datadir,fname,'/p_',filename(1:end-5),'.jpg');

    % imshow figure
    % disp(test_show_im) ;

    i
    run OPC_TEST.m

end