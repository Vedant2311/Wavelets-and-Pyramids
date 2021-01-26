image = imread('buttress.jpg');

orig = image;
[rows,~] = size(orig);
    
flag = 0;

G = randn(rows);
image = double(image) + 30*G;
flag = 1;


lambda_0 = image;
k = log(rows)/log(2);

%figure,imshow(uint8(image));

output = cell(1,k);
gammas = cell(1,k);

% Pyramid formation
for p=1:k
    [rows,columns] = size(lambda_0);
    lambda_1 = lambda_0(2:2:rows,2:2:columns);

    temp = lambda_0(1:2:rows,1:2:columns);
    gamma_1 = zeros(rows/2,columns/2);

    for i=1:rows/2

        for j=1:columns/2
            gamma_1(i,j) = temp(i,j) - 1/2 * (lambda_1(i,j) + lambda_1((min(i+1,rows/2)),(min(j+1,columns/2))));
        end
    end

    for i=1:rows/2

        for j=1:columns/2
            lambda_1(i,j) = lambda_1(i,j) + 1/4*(gamma_1(max(1,i-1),max(1,j-1)) + gamma_1(i,j));
        end
    end

    lambda_0 = lambda_1;
    output(1,p) = {lambda_1};
    gammas(1,p) = {gamma_1};
    
    %figure,imshow(lambda_1);
end

% Reconstruction

if flag==0 
    l = 1;
else 
    l = 2;
end

A = output(1,l);
lambda = [A{:}];

B = gammas(1,l);
gamma = [B{:}];
size(lambda)
lambda = reshape(lambda,[],1);
lambda = compress(lambda,20);
lambda = reshape(lambda,128,128);


for p= l : -1 : 1
    
    %{
    HARD thresholding: 62 -> 19.4499 (THE BEST VALUE)
                       90 ->  19.4213
    %}
    
    %{
    
    SOFT thresholding: 40 -> 20.2382 (THE BEST VALUE)
                       90 -> 19.3612
    
    %}
    
    %gamma = soft(gamma,28);

    [rows,columns] = size(lambda);
    rows = rows * 2;
    columns = columns * 2;
    
    for i=1:rows/2

        for j=1:columns/2
            lambda(i,j) = lambda(i,j) - 1/4*(gamma(max(1,i-1),max(1,j-1)) + gamma(i,j));
        end
    end

    temp = zeros(rows/2,columns/2);

    for i=1:rows/2

        for j=1:columns/2
            temp(i,j) = gamma(i,j) + 1/2*(lambda(i,j) + lambda(min(i+1,rows/2),min(j+1,columns/2)));
        end
    end

    recons = zeros(rows,columns);

    for i=1:rows/2
        for j=1:columns/2
            recons(2*i-1,2*j-1) = temp(i,j);
            recons(2*i,2*j) = lambda(i,j);
            recons(2*i,2*j-1) = (temp(i,j) + lambda(i,j))/2;
            recons(2*i-1,2*j) = (temp(i,j) + lambda(i,j))/2;        
        end
    end
    
    lambda = recons;

    B = gammas(1,max(p-1,1));
    gamma = [B{:}];

end
imwrite(uint8(recons),'ext.png')
figure,imshow(uint8(recons));

[psnr_out,snr_out] = psnr(uint8(recons),uint8(orig))

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

function[x] = compress(wavelet,k)
l = uint64(k*(512*512)/200);
v1 = maxk(wavelet,l);
v2 = mink(wavelet,l);
for i=1:length(wavelet)
    
        if ismember(wavelet(i),v1)==0 && ismember(wavelet(i),v2)==0
            wavelet(i) = 0;
        end
end
x = wavelet;
end
