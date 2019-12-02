image_in = imread('Image A.png');
image_in = image_in(1:256,1:256,:);
h1=[1/2 0 0 0 1/2 0 0 0; 
    1/2 0 0 0 -1/2 0 0 0;
    0 1/2 0 0 0 1/2 0 0;
    0 1/2 0 0 0 -1/2 0 0;
    0 0 1/2 0 0 0 1/2 0;
    0 0 1/2 0 0 0 -1/2 0;
    0 0 0 1/2 0 0 0 1/2;
    0 0 0 1/2 0 0 0 -1/2

    ];


h2=[1/2 0 1/2 0 0 0 0 0;
    1/2 0 -1/2 0 0 0 0 0;
    0 1/2 0 1/2 0 0 0 0;
    0 1/2 0 -1/2 0 0 0 0;
    0 0 0 0 1 0 0 0;
    0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 1
    ];

h3=[1/2 1/2 0 0 0 0 0 0;
    1/2 -1/2 0 0 0 0 0 0;
    0 0 1 0 0 0 0 0;
    0 0 0 1 0 0 0 0;
    0 0 0 0 1 0 0 0;
    0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 1
    
    ];

h1 = normc(h1);
h2 = normc(h2);
h3 = normc(h3);


h = h1 * h2 * h3;

haar = @(x) (h)'*x*(h);
image_in_R = image_in(:,:,1);
image_in_G = image_in(:,:,2);
image_in_B = image_in(:,:,3);

[rows,columns]=size(image_R);

image_in_cell_R=mat2cell( double(image_in_R) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_haar_cell_R=cellfun( haar, image_in_cell_R , 'UniformOutput',false);
image_haar_R= cell2mat(image_haar_cell_R);

image_in_cell_G=mat2cell( double(image_in_G) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_haar_cell_G=cellfun( haar, image_in_cell_G , 'UniformOutput',false);
image_haar_G= cell2mat(image_haar_cell_G);

image_in_cell_B=mat2cell( double(image_in_B) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_haar_cell_B=cellfun( haar, image_in_cell_B , 'UniformOutput',false);
image_haar_B= cell2mat(image_haar_cell_B);

w_R = reshape(image_haar_R,[],1);
w_R = compress(w_R,5);
[x_R,z_R] = encode(w_R);
w_R = reshape(w_R,256,256);

w_G = reshape(image_haar_G,[],1);
w_G = compress(w_G,5);
[x_G,z_G] = encode(w_G);
w_G = reshape(w_G,256,256);

w_B = reshape(image_haar_B,[],1);
w_B = compress(w_B,5);
[x_B,z_B] = encode(w_B);
w_B = reshape(w_B,256,256);

fileID = fopen('img.bin','w','a');
fwrite(fileID,z_R)
fwrite(fileID,z_G)
fwrite(fileID,z_B)

% fid = fopen('img.bin','r');
% format = 'string';
% b = fread(fid,Inf,format); % this one works
% fclose(fid);
% 
% b_R = b(1:65536);
% b_G = b(65537:65536*2);
% b_B = b(131073:65536*3);
% x_R = decompress(b_R);
% x_G = decompress(b_G);
% x_B = decompress(b_B);

%w = reshape(h,512,512);
invhaar = @(x) (h')'*x*(h');

image_haar_cell_R=mat2cell( double(w_R) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_out_cell_R=cellfun( invhaar, image_haar_cell_R , 'UniformOutput',false);
image_out_R= uint8(cell2mat(image_out_cell_R));

image_haar_cell_G=mat2cell( double(w_G) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_out_cell_G=cellfun( invhaar, image_haar_cell_G , 'UniformOutput',false);
image_out_G= uint8(cell2mat(image_out_cell_G));

image_haar_cell_B=mat2cell( double(w_B) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_out_cell_B=cellfun( invhaar, image_haar_cell_B , 'UniformOutput',false);
image_out_B= uint8(cell2mat(image_out_cell_B));

image_out = image_in;
image_out(:,:,1) = image_out_R;
image_out(:,:,2) = image_out_G;
image_out(:,:,3) = image_out_B;

[psnr_out,snr_out] = psnr((image_in),uint8(image_out))
figure,imshow(uint8(image_out));
imwrite(image_out,'compressed.png')

function[x] = compress(wavelet,k)
l = uint64(k*length(wavelet)/200);
v1 = maxk(wavelet,l);
v2 = mink(wavelet,l);
for i=1:length(wavelet)
    
        if ismember(wavelet(i),v1)==0 && ismember(wavelet(i),v2)==0
            wavelet(i) = 0;
        end
end
x = wavelet;
end

function [x,z] = encode(w)
i=1;
j=1;
l = length(w);
while i<=l
    if w(i)~=0
        y(j) = w(i)+0.0001;
        j=j+1;
        i=i+1;
    end
   c=0;
   if i<=l && w(i)==0
   while i<=l && w(i)==0
       c=c+1;
       i=i+1;
   end
   y(j)=uint8(c);
   j=j+1;
   end
   
end
u = num2str(y(1),8);
k=length(y)
for i=2 : k
   u = strcat(' ',u,num2str(y(i),8));
end
x=y;
z=u;
end

function [x] = decompress(w)
l = length(w);
i=1;
j=1;
while i<=l
    if ceil(w(i)) ~= w(i)
        y(j)= w(i);
        j=j+1;
    else
        for k=1:w(i)
            y(j) = 0;
            j=j+1;
        end
    end
end
x=y;
end


function y=f_d2b(n)
   
    
  
    
     
    strn=strtrim(num2str(n));
     
    %------------------------------------------------
    if isempty(find(strn=='.'))   
        y=d2b(n);
        if length(y)<12
         y = strcat(y,'.00000000000');
        end
        return;        
    else
    %------------------------------------------------
        k = find(strn =='.');
    end
    
    
    %Retrieving INTEGER and FRACTIONAL PARTS as strings
    i_part=strn(1:k-1);
    f_part=strn(k:end);
    
    
    %Converting the strings back to numbers
    ni_part=str2num(i_part);
    nf_part=str2num(f_part);
    
    ni_part=d2b(ni_part);
      
    strtemp='';
    
    
    temp=nf_part;
    %-------------------------------------------------
    t='1';s='0';
    
    while nf_part>= 0
        nf_part=nf_part*2;
        if (nf_part==1) || (nf_part==temp)
            strtemp=strcat(strtemp,t);
            break;
        elseif nf_part>1
            strtemp=strcat(strtemp,t);
            nf_part=nf_part-1;
        else
            strtemp=strcat(strtemp,s);
        end
    end
    
    
     
     if(i_part=='0')
        y=strcat('0.',strtemp);
         
     else
        y=strcat(ni_part,'.',strtemp);
     end
     if length(y)<12
         y = strcat(y,'00000000000');
     end
     
     %------------------------------------------------
end
    
function y=d2b(n)
strtemp='';
strn=strtrim(num2str(n));
if n<0
    fprintf(' %f is not a valid number\n',n)
 
end
  while n~=0
    strtemp=strcat(num2str(mod(n,2)),strtemp);
    n=floor(n/2);
  end
  y=strtemp;
end

function y=f_b2d(n)
    
   
    n=strtrim(n);
    numadd=0;
    t=length(n);
    %----------------------------------------------------------------------
    if isempty(find(n=='.'))
        y=bin2dec(n);     
    else
    %----------------------------------------------------------------------
        j=find(n=='.');
        i_part=n(1:j-1);
        f_part=n(j+1:end);
        
        % Look for the indices of the fraction part which are ones and
        % a simple computation on it to convert to decimal
        d_fpart=sum(0.5.^find(f_part=='1'));
        
        y=f_b2d(i_part)+d_fpart; % Concatenate the results
    %----------------------------------------------------------------------
    end
end

