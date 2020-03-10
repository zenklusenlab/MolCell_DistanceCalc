function [out,vals]=colocalize3(rna1, rna2, rna3, pixelsize,radiuss)
%%Return the number of colocalize rna from rna1 and rna2 and rna3
global res
global coloc_ind
global min_coloc_dif
radius=radiuss/pixelsize;
coloc_ind={};
n= numel(rna1(:,1));
res=zeros(n,9);
min_coloc_dif=zeros(n,9);
doColoc(rna1, rna2, rna3,radius);
% waitfor(hbut, 'UserData')
out=res;
out(:,9) = out(:,2)+out(:,5);
out(find(out(:,9)<2),:) = [];

out(find(out(:,3)>radius),:) = [];
out(find(out(:,6)>radius),:) = [];
out(find(out(:,8)>radius),:) = [];

u=unique(out(find(~isnan(out(:,4))>0),4));
n=histc(out(:,4),u);
d = u(n > 1);
vals = find(ismember(out(:,4),d));

out(vals,:) = [];

u=unique(out(find(~isnan(out(:,7))>0),7));
n=histc(out(:,7),u);
d = u(n > 1);
vals = find(ismember(out(:,7),d));

out(vals,:) = [];

out(:,9) = [];

each_coloc=coloc_ind;
min_coloc=min_coloc_dif;
min_coloc(vals,:) = [];


    function updateVal(~, ~)
        pix=get(hsl, 'Value');
        set(h,'MarkerSize', pix)
        set(htext, 'String', [num2str(pix*100*pixelsize/100), 'nm']);
    end

end

function doColoc(rna1, rna2, rna3,radius)
global res
global coloc_ind
global min_coloc_dif
% slider=findobj(0,'Tag','gslider');
% radius=get(slider, 'Value');
% disp(['Radius : ', num2str(radius*pixelsize),'nm']);
for i=1:length(res(:,1))
    dst1= double((rna2(:,1)-rna1(i,1)).^2 + (rna2(:,2)-rna1(i,2)).^2);
    coloc1= dst1<=(radius.^2);
%     summ1 = sum(coloc1);
    
    dst2= double((rna3(:,1)-rna1(i,1)).^2 + (rna3(:,2)-rna1(i,2)).^2);
    coloc2= dst2<=(radius.^2);
%     summ2 = sum(coloc2);
    
    [min_val1,min_index1] = min(dst1);
    [min_val2,min_index2] = min(dst2);
    
    dst3 = double((rna3(min_index2,1)-rna2(min_index1,1)).^2 + (rna3(min_index2,2)-rna2(min_index1,2)).^2);
    
    res(i,1:4)=double([i, logical(nnz(coloc1)), sqrt(min_val1), min_index1]);
    res(i,5:8)=double([logical(nnz(coloc2)), sqrt(min_val2), min_index2 sqrt(dst3)]);
    coloc_ind{i}=find(coloc1);
%     min_coloc_dif(i) = min_index;
%      rna2(min_index,1)-rna1(i,1)
%      rna2(min_index,2)-rna1(i,2)
    min_coloc_dif(i,1:4)=double([i, logical(nnz(coloc1)),rna2(min_index1,1)-rna1(i,1),rna2(min_index1,2)-rna1(i,2)]);
    min_coloc_dif(i,5:8)=double([logical(nnz(coloc2)),rna3(min_index2,1)-rna1(i,1),rna3(min_index2,2)-rna1(i,2) sqrt(dst3)]);
end
% close gcf
end