image = imread('buttress.jpg');

image = image(1:257,1:257);
image = double(image);

[rows,columns] = size(image);

%image_in = image + 30 * randn(rows);
image_in = image;
%figure,imshow(uint8(image_in));

a = 0.4;
b = 1/4;
c = 1/4 - a/2;
w_x = [c,b,a,b,c];
w_y = w_x';

w = kron(w_x,w_y);

l = floor(log(rows)/log(2));

image_gaussian = apply_gaussian(image_in,l,w);
image_laplacian = apply_laplacian(image_gaussian,a,b,c);

image_out = reconstruct(image_laplacian,a,b,c);
figure,imshow(uint8(image_out));
imwrite(uint8(image_out),'recons.png')
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
