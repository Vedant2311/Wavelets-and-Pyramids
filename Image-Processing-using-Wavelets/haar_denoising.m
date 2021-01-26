image = imread('buttress.jpg');

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
% image_R = image(:,:,1);
% image_G = image(:,:,2);
% image_B = image(:,:,3);

[rows,columns]=size(image);

G = randn(rows);
image_in = double(image) + 30*G;

%figure, imshow(uint8(image_in));
imwrite(uint8(image_in),'noisy.png')
image_in_cell=mat2cell( double(image_in) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_haar_cell=cellfun( haar, image_in_cell , 'UniformOutput',false);
image_haar= cell2mat(image_haar_cell);

%{
Threshold for soft: PSNR for 30 -> 23.2046
                         for 60 -> 22.4426
                         for 32 -> 23.3037 (BUT THE BEST IMAGE)
                         for 300 -> 14.8038 (But very pixelated!)

Threshold for hard: PSNR for 100 -> 22.7650 
                         for 60 ->  21.5848
                         for 67 -> 22.0990 (BUT THE BEST IMAGE)
                         for 300 -> 19.7624 (But very pixelated)

%}                          
image_haar = hard(image_haar,60);

invhaar = @(x) (h')'*x*(h');

image_haar_cell=mat2cell( double(image_haar) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_out_cell=cellfun( invhaar, image_haar_cell , 'UniformOutput',false);
image_out= uint8(cell2mat(image_out_cell));

%figure, imshow(uint8(image_out));
imwrite(uint8(image_out),'wave67.png')
[psnr_out,snr_out] = psnr(uint8(image),uint8(image_out))


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
