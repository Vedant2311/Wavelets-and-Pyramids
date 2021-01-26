imageA = imread('nikhil.png');
imageB = imread('goel.png');
mask = imread('mask1.png');

imageA_R = imageA(:,:,1);
imageA_R = double(imageA_R);
imageB_R = imageB(:,:,1);
imageB_R = double(imageB_R);

x = size(imageA_R)

imageA_G = imageA(1:257,1:257,2);
imageA_G = double(imageA_G);
imageB_G = imageB(1:257,1:257,2);
imageB_G = double(imageB_G);

imageA_B = imageA(1:257,1:257,3);
imageA_B = double(imageA_B);
imageB_B = imageB(1:257,1:257,3);
imageB_B = double(imageB_B);

mask = mask(1:257,1:257);
mask = double(mask);


[rows,columns] = size(imageA_R);

%image_in = image + 30 * randn(rows);
image_inA_R = imageA_R;
image_inB_R = imageB_R;

image_inA_G = imageA_G;
image_inB_G = imageB_G;

image_inA_B = imageA_B;
image_inB_B = imageB_B;

mask_in = mask; 
%figure,imshow(uint8(image_in));

a = 0.4;
b = 1/4;
c = 1/4 - a/2;
w_x = [c,b,a,b,c];
w_y = w_x';

w = kron(w_x,w_y);

l = floor(log(rows)/log(2));
size(image_inA_R)
image_gaussianA_R = apply_gaussian(image_inA_R,l,w);
image_laplacianA_R = apply_laplacian(image_gaussianA_R,a,b,c);

image_gaussianA_G = apply_gaussian(image_inA_G,l,w);
image_laplacianA_G = apply_laplacian(image_gaussianA_G,a,b,c);

image_gaussianA_B = apply_gaussian(image_inA_B,l,w);
image_laplacianA_B = apply_laplacian(image_gaussianA_B,a,b,c);

image_gaussianB_R = apply_gaussian(image_inB_R,l,w);
image_laplacianB_R = apply_laplacian(image_gaussianB_R,a,b,c);

image_gaussianB_G = apply_gaussian(image_inB_G,l,w);
image_laplacianB_G = apply_laplacian(image_gaussianB_G,a,b,c);

image_gaussianB_B = apply_gaussian(image_inB_B,l,w);
image_laplacianB_B = apply_laplacian(image_gaussianB_B,a,b,c);

mask_gaussian = apply_gaussian(mask_in,l,w);
mask_laplacian = apply_laplacian(mask_gaussian,a,b,c);

LS_R = {};
for i=1:9
    for j=1:length(image_laplacianA_R{i})
        for k=1:length(image_laplacianA_R{i})
           % mask_gaussian{i}(j,k)
            LS_R{i}(j,k) = ((mask_gaussian{i}(j,k)/255)*image_laplacianB_R{i}(j,k)) + ((1-(mask_gaussian{i}(j,k))/255)*image_laplacianA_R{i}(j,k));
        end
    end
end

LS_G = {};
for i=1:9
    for j=1:length(image_laplacianA_G{i})
        for k=1:length(image_laplacianA_G{i})
           % mask_gaussian{i}(j,k)
            LS_G{i}(j,k) = ((mask_gaussian{i}(j,k)/255)*image_laplacianB_G{i}(j,k)) + ((1-(mask_gaussian{i}(j,k))/255)*image_laplacianA_G{i}(j,k));
        end
    end
end

LS_B = {};
for i=1:9
    for j=1:length(image_laplacianA_B{i})
        for k=1:length(image_laplacianA_B{i})
           % mask_gaussian{i}(j,k)
            LS_B{i}(j,k) = ((mask_gaussian{i}(j,k)/255)*image_laplacianB_B{i}(j,k)) + ((1-(mask_gaussian{i}(j,k))/255)*image_laplacianA_B{i}(j,k));
        end
    end
end

image_out_R = reconstruct(LS_R,a,b,c);
image_out_G = reconstruct(LS_G,a,b,c);
image_out_B = reconstruct(LS_B,a,b,c);
image_out = imageA;
image_out(:,:,1) = image_out_R;
image_out(:,:,2) = image_out_G;
image_out(:,:,3) = image_out_B;
figure,imshow(uint8(image_out));
imwrite(image_out,'spline.png')

%imwrite(uint8(image_out),'recons.png')
%[psnr_out,snr_out] = psnr(uint8(image),uint8(image_out));
%psnr_out

%print_cell(image_laplacian);

function E = reconstruct(cellM,a,b,c)

    [~,s] = size(cellM);
    
    for i = s-1 : -1 : 1
        A = cellM(1,i);
        A = [A{:}];
        
        B = cellM(1,i+1);
        B = [B{:}];
        
        F = A + expand(B,A,a,b,c);
       % F = hard(F,20);
        cellM(1,i) = {F};
        
    end
    
    E = cellM(1,1);
    
    E = [E{:}];
end

function A = apply_gaussian(image,l,w)

     A = cell(1,1+l);
     C = image;
     
     A(1,1) = {C};
     
     for c=1:l
        [rows,columns] = size(C);
    
         B = zeros((rows-1)/2+1 ,(columns-1)/2 +1 );
         
         for i=1:((rows-1)/2 +1 )
            
             for j=1:((columns-1)/2 +1 )
                 sum = 0;
                 for m=-2:2
                     for n=-2:2
                         
                         if ((2*i-1+m >=1) && (2*i-1+m <=rows) && (2*j-1+n >=1) && (2*j-1+n <=columns))
                             
                            sum = sum +  w(m+3,n+3) * C(2*i-1+m,2*j-1+n);
                         end
                     end
                 end
                 
                 B(i,j) = sum;
                 
             end
         end
         
         A(1,c+1) = {B};
         C = B;
         
        
     end
    
end

function print_cell(cell)

    [~,s] = size(cell);
    
    for i=1:s
        
        A = cell(1,i);
        A = [A{:}];
        figure,imshow(uint8(A));
    end

end

function B = apply_laplacian(cellM,a,b,c)
     
    [~,s] = size(cellM);
    B = cell(1,s);
    
    K = cellM(1,s);
    K = [K{:}];
    
    
    B(1,s) = {K};
    for i=s-1:-1:1
        
        A = cellM(1,i);
        A = [A{:}];
        
        C = cellM(1,i+1);
        C = [C{:}];
        
        D = expand(C,A,a,b,c);
        
        E = A -D;
        
        
        % For hard, the best level for thresholding was achieved  at level 2.
        % The best value at 120: PSNR = 23.7838
        
        if i<2
            E = soft(E,120);
        end
      
        
        
        % For soft, the best level: 2
        % Best value at 60: PSNR = 24.0497
        %{
        if i<2
            E = soft(E,60);
        end
        %}
        
        B(1,i) = {E};
        
        
    end 

end

function B = expand(C,A,a,b,c)

    k_x = [c b a b c];
    kernel = kron(k_x,k_x') * 4;

    ker00 = kernel(1:2:5,1:2:5); %2*i,2*i
    ker01 = kernel(1:2:5,2:2:5); %2*i,2*i+1
    
    ker10 = kernel(2:2:5,1:2:5); %2*i+1,2*i
    ker11 = kernel(2:2:5,2:2:5); %2*i+1, 2*i+1

    [rows,columns] = size(A);
    B = zeros(rows,columns);

    C_h = padarray(C,[0 1],'replicate'); 
	C_v = padarray(C,[1 0],'replicate'); 
	
	img00 = imfilter(C,ker00);
    
    
	img01 = conv2(C_v,ker01,'valid') ;
    
	img10 = conv2(C_h,ker10,'valid');
	img11 = conv2(C,ker11,'valid');
	
	B(1:2:rows,1:2:columns) = img00;
	B(2:2:rows,1:2:columns) = img10;
	B(1:2:rows,2:2:columns) = img01;
	B(2:2:rows,2:2:columns) = img11;
    
    %figure,imshow(uint8(B));
end


function A = hard(w,t)

[r,c] = size(w);

A = zeros(r,c);

for i=1:r
    for j=1:c
        
        if (abs(w(i,j)) > t)
            A(i,j) = w(i,j);
        else
            A(i,j) = 0;
        end
        
    end
end

end

function A = soft(w,t)

[r,c] = size(w);

A = zeros(r,c);

for i=1:r
    for j=1:c
        
        if (w(i,j)) > t
            A(i,j) = w(i,j) -t;
        elseif w(i,j) < -t
            A(i,j) = w(i,j) +t;
        else
            A(i,j) = 0;
        end
        
    end
end

end
