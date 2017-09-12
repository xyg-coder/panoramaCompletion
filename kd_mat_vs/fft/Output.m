function [  ] = Output( matrix,path_name )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
file=fopen(path_name,'wt');
[row,col]=size(matrix);
for i=1:row
    for j=1:col
        fprintf(file,'%d ',matrix(i,j));
    end;
    fprintf(file,'\n');
end;

end

