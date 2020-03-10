function RNA_coloc3(mask1, mask2, mrnafile , pixelsize, radius, dist, alchk,reference)
%%%%%%%%%% Latest version 3 color RNA quantification %%%%%%%%%%%%%%%%%
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
% 
% mrnas = nchoosek(mrnafile,2);
% filenames = ["Cy5_Cy3";"Cy5_Dy488";"Cy3_Dy488"];
mrna5file = mrnafile(1);
mrnamidfile = mrnafile(2);
mrna3file = mrnafile(3);
% files = filenames(qu);

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
if(~isempty(mrnamidfile))
    mrnamid=load(mrnamidfile);
    mrnamid=mrnamid((mrnamid(:,3)>=threshold),indexing);
    
    for i = 1:length(mrnamid)
        D = zeros(size(mrnamid,1),1);
        D(i+1:end) = double((mrnamid(i+1:end,1)-mrnamid(i,1)).^2 + (mrnamid(i+1:end,2)-mrnamid(i,2)).^2);
        D = sqrt(D);
        D = D*pixelsize;
        out =find(D(i+1:end)<dist);
        if ~isempty(out)
            xxx = [xxx;i;out+i];
       end
    end
    xxx = unique(xxx);
    mrnamid(xxx,:) = [];
    mrnamidcyt=nucleus(mrnamid, img1, inten_dist1);
    mrnamidcyt=mrnamidcyt(mrnamidcyt(:,4)>=1,1:3);
    if ~isempty(img2)
        mrnamidnuc=nucleus(mrnamid, img2, inten_dist2);
        mrnamidnuc=mrnamidnuc(mrnamidnuc(:,4)>=1,1:3);
    else
        mrnamidnuc=[];
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
% 
% clustercyt = [mrna5cyt(:,1:2) ones(length(mrna5cyt),1)*100000;mrnamidcyt(:,1:2) ones(length(mrnamidcyt),1)*100000+100000; mrna3cyt(:,1:2) ones(length(mrna3cyt),1)*100000+200000];
% aa = kmeans(clustercyt,max([length(mrna5cyt) length(mrna3cyt) length(mrnamidcyt)]));
% aaa = aa(length(mrna5cyt)+1:length(mrna5cyt)+length(mrnamidcyt));
% length(unique(aaa))

if(reference == "mid")
    rna1cyt = mrnamidcyt;
    rna2cyt = mrna5cyt;
    rna3cyt = mrna3cyt;
    
    rna1nuc = mrnamidnuc;
    rna2nuc = mrna5nuc;
    rna3nuc = mrna3nuc;
    
elseif(reference == "5p")
    rna1cyt = mrna5cyt;
    rna2cyt = mrnamidcyt;
    rna3cyt = mrna3cyt;
    
    rna1nuc = mrna5nuc;
    rna2nuc = mrnamidnuc;
    rna3nuc = mrna3nuc;  
elseif(reference == "3p")
    rna1cyt = mrna3cyt;
    rna2cyt = mrna5cyt;
    rna3cyt = mrnamidcyt;
    
    rna1nuc = mrna3nuc;
    rna2nuc = mrna5nuc;
    rna3nuc = mrnamidnuc;
else
    cprintf('err','Set reference');
    return 
end

if(alchk==1)
    if (~isempty(mrna5file) &&  ~isempty(mrna3file))
         if ~isempty(img2)
             vals=[];
             [mrna5_coloc_mrna3_nuc,mrna5_mrna3_coloc_val_nuc,vals]= colocalize2_3colornew(mrna5nuc, mrna3nuc,pixelsize,radius);
             mrna5nuc(vals,:) = [];
         else
             [mrna5_coloc_mrna3_nuc] = [];
             [mrna5_mrna3_coloc_val_nuc] = [];
         end
        [mrna5_coloc_mrna3_cyt,mrna5_mrna3_coloc_val_cyt,vals]= colocalize2_3colornew(mrna5cyt, mrna3cyt,pixelsize,radius);
        mrna5cyt(vals,:) = [];
        twospotInput(mrna5cyt, mrna3cyt, mrna5_coloc_mrna3_cyt, mrna5nuc, mrna3nuc, mrna5_coloc_mrna3_nuc, pixelsize,img2,files);
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
    if (~isempty(mrna5file) &&  ~isempty(mrna3file) && ~isempty(mrna3file))
        if ~isempty(img2) 
            vals=[];
            [coloc_mrna5mid3_nuc,vals]= colocalize3(rna1nuc, rna2nuc,rna3nuc,pixelsize,radius);
%             mrna5nuc(vals,:) = [];
        else
            [coloc_mrna5mid3_nuc] = [];
        end
%         
        vals=[];
        [coloc_mrna5mid3_cyt,vals]= colocalize3(rna1cyt, rna2cyt,rna3cyt, pixelsize,radius);
       
        %           mrna5cyt(vals,:) = [];
        rna1cyts = rna1cyt(coloc_mrna5mid3_cyt(:,1),:);
        rna2cyts = rna2cyt(coloc_mrna5mid3_cyt(:,4),:);
        rna3cyts = rna3cyt(coloc_mrna5mid3_cyt(:,7),:);
        
        rna1nucs = rna1nuc(coloc_mrna5mid3_nuc(:,1),:);
        rna2nucs = rna2nuc(coloc_mrna5mid3_nuc(:,4),:);
        rna3nucs = rna3nuc(coloc_mrna5mid3_nuc(:,7),:);
        
        twospotInput(rna1cyts, rna2cyts, rna3cyts, coloc_mrna5mid3_cyt, rna1nucs, rna2nucs,rna3nucs, coloc_mrna5mid3_nuc, pixelsize,img2,reference);
    end

end

disp('Done')

end

%% Cas de colocalisation a 2 (erna, spot1)
function twospotInput(spot1cyt, spot2cyt, spot3cyt, spot1_spot2_spot3_cyt,spot1nuc, spot2nuc, spot3nuc, spot1_spot2_spot3_nuc, pixelsize,img2,reference)

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
mrna5_3Data1= zeros(size(spot1cyt,1),7);
mrna_spots_cyt= zeros(size(spot1cyt,1),7);
a=size(spot1cyt,1);
if(a>0)
    mrna5_3Data1(:,1) = 1:1:size(spot1cyt,1);
    mrna5_3Data1(:,2) = spot1cyt(:,3);
    mrna5_3Data1(:,3) = spot2cyt(:,3);
    mrna5_3Data1(:,4) = spot3cyt(:,3);
    mrna5_3Data1(:,5) = spot1_spot2_spot3_cyt(:,3).*pixelsize;
    mrna5_3Data1(:,6) = spot1_spot2_spot3_cyt(:,6).*pixelsize;
    mrna5_3Data1(:,7) = spot1_spot2_spot3_cyt(:,8).*pixelsize;

    mrna_spots_cyt(:,1) = 1:1:size(spot1cyt,1);
    mrna_spots_cyt(:,2:3) = spot1cyt(:,1:2);
    mrna_spots_cyt(:,4:5) = spot2cyt(:,1:2);
    mrna_spots_cyt(:,6:7) = spot3cyt(:,1:2);

    spotss11 = find(mrna5_3Data1(:,7)>300);
    mrna5_3Data1(spotss11,:) = [];
    mrna_spots_cyt(spotss11,:) = [];
    mrna531 = num2cell(mrna5_3Data1);
    mrnaspotscyt = num2cell(mrna_spots_cyt);

    C = cell(length(mrna531), 1);
    if(reference == "mid")
        C(:) = {'mid'};
    elseif(reference == "5p") 
        C(:) = {'5p'};
    elseif(reference == "3p") 
        C(:) = {'3p'};    
    end
    
    mrna531(:,8) = C;
    mrnaspotscyt(:,8) = C;
    
end



if ~isempty(img2)
% filesname = sprintf(,files);
% csvwrite('Cytoplasmic Distances.csv',mrna5_3Data1);
cell2csv('Cytoplasmic Distances.csv',mrna531,','); 
cell2csv('Cytoplasmic Spots.csv',mrnaspotscyt);
else
csvwrite('Distances.csv',mrna5_3Data1);    
end

% header={'#','5 intensity','3 intensity','Nuclear_Distance' };

if ~isempty(img2)
mrna5_3Data2= zeros(size(spot1nuc,1),6);
mrna_spots_nuc= zeros(size(spot1nuc,1),7);
a=size(spot1nuc,1);
if(a>0)
    mrna5_3Data2(:,1) = 1:1:size(spot1nuc,1);
    mrna5_3Data2(:,2) = spot1nuc(:,3);
    mrna5_3Data2(:,3) = spot2nuc(:,3);
    mrna5_3Data2(:,4) = spot3nuc(:,3);
    mrna5_3Data2(:,5) = spot1_spot2_spot3_nuc(:,3).*pixelsize;
    mrna5_3Data2(:,6) = spot1_spot2_spot3_nuc(:,6).*pixelsize;
    mrna5_3Data2(:,7) = spot1_spot2_spot3_nuc(:,8).*pixelsize;
    
    mrna_spots_nuc(:,1) = 1:1:size(spot1nuc,1);
    mrna_spots_nuc(:,2:3) = spot1nuc(:,1:2);
    mrna_spots_nuc(:,4:5) = spot2nuc(:,1:2);
    mrna_spots_nuc(:,6:7) = spot3nuc(:,1:2);
    
    spotss22 = find(mrna5_3Data2(:,7)>300);
    mrna5_3Data2(spotss22,:) = [];
    mrna_spots_nuc(spotss22,:) = [];
    
%     mrna5_3Data2(find(mrna5_3Data1(:,7)>300),:) = [];
    mrna532 = num2cell(mrna5_3Data2);
    mrnaspotsnuc = num2cell(mrna_spots_nuc);


    C = cell(length(mrna532), 1);
    if(reference == "mid")
        C(:) = {'mid'};
    elseif(reference == "5p") 
        C(:) = {'5p'};
    elseif(reference == "3p") 
        C(:) = {'3p'};    
    end
    
    mrna532(:,8) = C;
    mrnaspotsnuc(:,8) = C;
    
    
    
% filesname2 = sprintf('Nuclear Distances.csv',files);
% csvwrite('Nuclear Distances.csv',mrna5_3Data2);
cell2csv('Nuclear Distances.csv',mrna532,',');
cell2csv('Nuclear Spots.csv',mrnaspotsnuc);

else
    mrna5_3Data2 = [];
end

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