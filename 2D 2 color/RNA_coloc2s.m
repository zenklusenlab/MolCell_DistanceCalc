function [cytval,nucval] = RNA_coloc2s(mask1, mask2, mrna5file, mrna3file, pixelsize, radius, dist, alchk,nucval,cytval)
%ARN_coloc, took 4 to 6 , Use:
%**ARN_coloc(imgfile, mrna3file, midmrnafile, mrna5file)
%**ARN_coloc(imgfile, mrna3file, midmrnafile, mrna5file,k)
%**ARN_coloc(imgfile, mrna3file, midmrnafile, mrna5file,k, threshold)
%  - imgfile = a string indicating the path to the nucleus mask file (tif or tiff)
%  - mrna3file= a string indicating the path to the sense erna loc file
%  - midmrnafile= a string indicating the path to the anti-sense erna loc file
%  - mrna5file= a string indicating the path to the mrna loc file
%  - k, this is the intensity coefficient to use for transcription site
%  detection, not mandatory, default=1.5
%  - threshold, this is not mandatory,it's used to discard any spot with intensity
%below threshold, default value =-1
%
%3 files will be saved:
%  - 'final_label.png' : an image labelling nucleus and spot:
%        +red spot = mrna
%        +blue spot = sense erna
%        +green spot = anti-sense erna
%  - 'trans_coloc_analysis.txt' : output with only transcription site mrna
%  - 'spot_coloc_analysis.txt' : output with all spot analysis
%
% EXAMPLE:
% ARN_coloc('mask.tif', 'mrna3.loc', 'midmrna.loc', 'mRNA.loc' )


narginchk(4,10);
threshold=-1;
if ~exist('radius', 'var') || isempty(radius)
    radius=300;
end

if ~exist('pixelsize', 'var') || isempty(pixelsize)
    pixelsize=39.682539;
end

% if sum(mask2=='\')>=1
%     mask2=[];
% end

img1=imread(mask1); % read the image tif file - can be the mask 
img1=bwlabel(im2bw(mat2gray(img1),0),4); %labeling different regions

if ~isempty(mask2)
img2=imread(mask2);
img2=bwlabel(im2bw(mat2gray(img2),0),4); %labeling different regions
else
    img2=[];
end

%TRouver la nouvelle distribution en intensité des noyaux
inten_dist1=int64(unique(sort(img1(img1>0)))); %almost useless
inten_dist2=int64(unique(sort(img2(img2>0))));
indexing=floor(1:3);

%Lecture des fichier de coordonnées des arn
xxx=[];

if(~isempty(mrna3file))
    mrna3=load(mrna3file);
    mrna3=mrna3((mrna3(:,3)>=threshold),indexing);
    
    for i = 1:length(mrna3)
        D = zeros(size(mrna3,1),1);
        D(i+1:end) = double((mrna3(i+1:end,1)-mrna3(i,1)).^2 + (mrna3(i+1:end,2)-mrna3(i,2)).^2);
        D = sqrt(D);
        D = D*pixelsize;
        out =find(D(i+1:end)<dist);
        if ~isempty(out)
            xxx = [xxx;i;out+i];
       end
    end
    xxx = unique(xxx);
    mrna3(xxx,:) = [];
    mrna3cyt=nucleus(mrna3, img1, inten_dist1);
    mrna3cyt=mrna3cyt(mrna3cyt(:,4)>=1,1:3);
    if ~isempty(img2)
        mrna3nuc=nucleus(mrna3, img2, inten_dist2);
        mrna3nuc=mrna3nuc(mrna3nuc(:,4)>=1,1:3);
    else
        mrna3nuc=[];
    end
end

xxx=[];

if(~isempty(mrna5file))
    mrna5=load(mrna5file);
    mrna5=mrna5((mrna5(:,3)>=threshold),indexing);
   
       for i = 1:length(mrna5)
        D = zeros(size(mrna5,1),1);
      	D(i+1:end) = double((mrna5(i+1:end,1)-mrna5(i,1)).^2 + (mrna5(i+1:end,2)-mrna5(i,2)).^2);
        D = sqrt(D);
        D = D*pixelsize;
        out =find(D(i+1:end)<dist);
        if ~isempty(out)
            xxx = [xxx;i;out+i];
        end
       end
    xxx = unique(xxx);
    mrna5(xxx,:) = [];
    if ~isempty(img2)
        mrna5nuc=nucleus(mrna5, img2, inten_dist2);  
        mrna5nuc=mrna5nuc(mrna5nuc(:,4)>=1,1:3);
    else
        mrna5nuc =[];
    end
    
    mrna5cyt=nucleus(mrna5, img1, inten_dist1);
    mrna5cyt=mrna5cyt(mrna5cyt(:,4)>=1,1:3);
end

if(alchk==1)
    if (~isempty(mrna5file) &&  ~isempty(mrna3file))
         if ~isempty(img2)
             vals=[];
             [mrna5_coloc_mrna3_nuc,mrna5_mrna3_coloc_val_nuc,vals]= colocalize2(mrna5nuc, mrna3nuc,pixelsize,radius);
             mrna5nuc(vals,:) = [];
         else
             [mrna5_coloc_mrna3_nuc] = [];
             [mrna5_mrna3_coloc_val_nuc] = [];
         end
        [mrna5_coloc_mrna3_cyt,mrna5_mrna3_coloc_val_cyt,vals]= colocalize2(mrna5cyt, mrna3cyt,pixelsize,radius);
        mrna5cyt(vals,:) = [];
        twospotInput(mrna5cyt, mrna3cyt, mrna5_coloc_mrna3_cyt, mrna5nuc, mrna3nuc, mrna5_coloc_mrna3_nuc, pixelsize,img2);
    end
    if ~isempty(img2)
        mrna5_mrna3_coloc_val_nuc = mrna5_mrna3_coloc_val_nuc(mrna5_mrna3_coloc_val_nuc(:,2)==1,:);
    else
        mrna5_mrna3_coloc_val_nuc =[];
    end
    mrna5_mrna3_coloc_val_cyt = mrna5_mrna3_coloc_val_cyt(mrna5_mrna3_coloc_val_cyt(:,2)==1,:);

    if ~isempty(img2)
        mrna5_mrna3_coloc_val_nuc_fin(:,3:4) =  mrna5_mrna3_coloc_val_nuc(:,3:4).*pixelsize;
    else
        mrna5_mrna3_coloc_val_nuc_fin=[];
    end
    
    mrna5_mrna3_coloc_val_cyt_fin(:,3:4) =  mrna5_mrna3_coloc_val_cyt(:,3:4).*pixelsize;

    if ~isempty(img2)
        pixel_shift = [mrna5_mrna3_coloc_val_nuc_fin(:,3:4);mrna5_mrna3_coloc_val_cyt_fin(:,3:4)];
    else
         pixel_shift = [mrna5_mrna3_coloc_val_cyt_fin(:,3:4)];
    end
    
%     header={'#','yshift','xshift'};
    pixshiftval= zeros(size(mrna5_mrna3_coloc_val_nuc_fin,1)+size(mrna5_mrna3_coloc_val_cyt_fin,1),4);
    pixshiftval(:,1) = 1:1:size(mrna5_mrna3_coloc_val_nuc_fin,1)+size(mrna5_mrna3_coloc_val_cyt_fin,1);
    pixshiftval(:,2:3) =  pixel_shift;
    
    csvwrite('Pixel shift.csv',pixshiftval);

    

elseif(alchk==0)
    if (~isempty(mrna5file) &&  ~isempty(mrna3file))
        if ~isempty(img2) 
            vals=[];
            [mrna5_coloc_mrna3_nuc,mincol,vals]= colocalize2(mrna5nuc, mrna3nuc,pixelsize,radius);
%             mrna5nuc(vals,:) = [];
        else
            [mrna5_coloc_mrna3_nuc] = [];
        end
%         
        vals=[];
        [mrna5_coloc_mrna3_cyt,mincol,vals]= colocalize2(mrna5cyt, mrna3cyt,pixelsize,radius);
%         mrna5cyt(vals,:) = [];
        twospotInput(mrna5cyt, mrna3cyt, mrna5_coloc_mrna3_cyt, mrna5nuc, mrna3nuc, mrna5_coloc_mrna3_nuc, pixelsize,img2);
    end

end

mrna5_coloc_mrna3_nuc = mrna5_coloc_mrna3_nuc(find(mrna5_coloc_mrna3_nuc(:,2)>0),:);
mrna5_coloc_mrna3_cyt = mrna5_coloc_mrna3_cyt(find(mrna5_coloc_mrna3_cyt(:,2)>0),:);


nucval = [nucval;mrna5nuc(mrna5_coloc_mrna3_nuc(:,1),1:2) mrna3nuc(mrna5_coloc_mrna3_nuc(:,4),1:2) mrna5_coloc_mrna3_nuc(:,3).*pixelsize];
cytval = [cytval;mrna5cyt(mrna5_coloc_mrna3_cyt(:,1),1:2) mrna3cyt(mrna5_coloc_mrna3_cyt(:,4),1:2) mrna5_coloc_mrna3_cyt(:,3).*pixelsize];

disp('Done')

end

%% Cas de colocalisation a 2 (erna, spot1)
function twospotInput(spot1cyt, spot2cyt, spot1_coloc_spot2_cyt,spot1nuc, spot2nuc, spot1_coloc_spot2_nuc, pixelsize,img2)

%data: 'Nuc', 'mrna3-spot1', #spot1	#mrna3
%trans_data: 'Nuc', 'mrna3-spot1', trans_number, nascent_mrna, #mrna3, mRNAnascent_coloc_+message,	#mRNAnascent_noColoc_erna
% mrnaDATAcyt=double([(1:size(spot1cyt,1))' spot1cyt spot1_coloc_spot2_cyt(:,end-1:end)]);
% spot1_w_spot2_cyt=mrnaDATAcyt(spot1_coloc_spot2_cyt(:,2)==1,:);
% spot1_wo_spot2_cyt= mrnaDATAcyt(spot1_coloc_spot2_cyt(:,1)==0,1:end-2);

% if ~isempty(img2)
%     mrnaDATAnuc=double([(1:size(spot1nuc,1))' spot1nuc spot1_coloc_spot2_nuc(:,end-1:end)]);
%     spot1_w_spot2_nuc=mrnaDATAnuc(spot1_coloc_spot2_nuc(:,2)==1,:);
% else
%     mrnaDATAnuc=[];
%     spot1_w_spot2_nuc = [];
% end
% spot1_wo_spot2_nuc= mrnaDATAnuc(spot1_coloc_spot2_nuc(:,1)==0,1:end-2);

%writing trans with spot2
% 
% header={'#','5 intensity','3 intensity','Cytoplasmic_Distance' };
mrna5_3Data1= zeros(size(spot1_coloc_spot2_cyt,1),4);
a=size(spot1_coloc_spot2_cyt,1);
if(a>0)
    mrna5_3Data1(:,1) = 1:1:size(spot1_coloc_spot2_cyt,1);
    mrna5_3Data1(:,2) = spot1cyt(spot1_coloc_spot2_cyt(:,1),3);
    mrna5_3Data1(:,3) = spot2cyt(spot1_coloc_spot2_cyt(:,4),3);
    mrna5_3Data1(:,4) = spot1_coloc_spot2_cyt(:,3).*pixelsize;
end

if ~isempty(img2)
csvwrite('Cytoplasmic Distances.csv',mrna5_3Data1);
else
csvwrite('Distances.csv',mrna5_3Data1);    
end

% header={'#','5 intensity','3 intensity','Nuclear_Distance' };

if ~isempty(img2)
    mrna5_3Data2= zeros(size(spot1_coloc_spot2_nuc,1),4);
    a=size(spot1_coloc_spot2_nuc,1);
    if(a>0)
        mrna5_3Data2(:,1) = 1:1:size(spot1_coloc_spot2_nuc,1);
        mrna5_3Data2(:,2) = spot1nuc(spot1_coloc_spot2_nuc(:,1),3);
        mrna5_3Data2(:,3) = spot2nuc(spot1_coloc_spot2_nuc(:,4),3);
        mrna5_3Data2(:,4) = spot1_coloc_spot2_nuc(:,3).*pixelsize;
    end
csvwrite('Nuclear Distances.csv',mrna5_3Data2);

else
    mrna5_3Data2 = [];
end

end

function coor=nucleus(coor, label_img, nuc_int)
coor(:,end+1)= zeros(size(coor(:,end)));
for i=1:length(coor(:,1))
    i_nuc=label_img(round(coor(i,2)),round(coor(i,1)));
    proche_nuc=zeros(2,2);
    if(coor(i,2)>1 && coor(i,1)>1)
        proche_nuc= label_img(round(coor(i,2))-1:round(coor(i,2))+1,round(coor(i,1))-1:round(coor(i,1))+1);
    end
    if(i_nuc>0)
        coor(i,end)= find(nuc_int==i_nuc);
    elseif sum(proche_nuc(:))>0
        coor(i, end)= find(nuc_int==proche_nuc(find(proche_nuc~=0,1)));
    end
end
end