function [out,min_coloc,vals]=colocalize2_3d(mrna1, mrna2, pixelsize,zpixelsize,radiuss)
%%Return the number of colocalize arn from rna1 and rna2
global res
global coloc_ind
global min_coloc_dif
rna1 = mrna1;
rna2 = mrna2;
radius=radiuss;
coloc_ind={};
n= numel(rna1(:,1));
res=zeros(n,4);
min_coloc_dif=zeros(n,4);
doColoc(rna1, rna2, pixelsize,zpixelsize, radius);
% waitfor(hbut, 'UserData')
out=res;

u=unique(out(find(~isnan(out(:,4))>0),4));
n=histc(out(:,4),u);
d = u(n > 1);
vals = find(ismember(out(:,4),d));

out(vals,:) = [];

out(find(out(:,2)<1),:) = [];

each_coloc=coloc_ind;
min_coloc=min_coloc_dif;
min_coloc(vals,:) = [];


    function updateVal(~, ~)
        pix=get(hsl, 'Value');
        set(h,'MarkerSize', pix)
        set(htext, 'String', [num2str(pix*100*pixelsize/100), 'nm']);
    end

end

function doColoc(rna1, rna2, pixelsize,zpixelsize,radius)
global res
global coloc_ind
global min_coloc_dif
% slider=findobj(0,'Tag','gslider');
% radius=get(slider, 'Value');
% disp(['Radius : ', num2str(radius*pixelsize),'nm']);
for i=1:length(res(:,1))
    dst= double(((rna2(:,1)-rna1(i,1))*pixelsize).^2 + ((rna2(:,2)-rna1(i,2))*pixelsize).^2+ ((rna2(:,3)-rna1(i,3))*zpixelsize).^2);
    coloc= dst<=(radius.^2);
    summ = sum(coloc);
    [min_val,min_index] = min(dst);
    res(i,1:4)=double([i, sum(coloc), sqrt(min_val), min_index]);
    coloc_ind{i}=find(coloc);
%     min_coloc_dif(i) = min_index;
%      rna2(min_index,1)-rna1(i,1)
%      rna2(min_index,2)-rna1(i,2)
    min_coloc_dif(i,1:5)=double([i, sum(coloc),rna2(min_index,1)-rna1(i,1),rna2(min_index,2)-rna1(i,2),rna2(min_index,3)-rna1(i,3)]);
end
% close gcf
end