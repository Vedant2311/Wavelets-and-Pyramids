image = imread('Image A.png');
image = image(1:256,1:256,:);
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
image_R = image(:,:,1);
image_G = image(:,:,2);
image_B = image(:,:,3);

[rows,columns]=size(image_R);

G = randn(rows);
image_in_R = double(image_R) + 30*G;
image_in_G = double(image_G) + 30*G;
image_in_B = double(image_B) + 30*G;

image_in = image;

image_in(:,:,1) = image_in_R;
image_in(:,:,2) = image_in_G;
image_in(:,:,3) = image_in_B;


%figure, imshow(uint8(image_in));
imwrite(uint8(image_in),'noisy.png')

image_in_cell_R=mat2cell( double(image_in_R) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_haar_cell_R=cellfun( haar, image_in_cell_R , 'UniformOutput',false);
image_haar_R= cell2mat(image_haar_cell_R);

image_in_cell_G=mat2cell( double(image_in_G) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_haar_cell_G=cellfun( haar, image_in_cell_G , 'UniformOutput',false);
image_haar_G= cell2mat(image_haar_cell_G);

image_in_cell_B=mat2cell( double(image_in_B) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_haar_cell_B=cellfun( haar, image_in_cell_B , 'UniformOutput',false);
image_haar_B= cell2mat(image_haar_cell_B);

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
image_haar_R = hard(image_haar_R,80);
image_haar_G = hard(image_haar_G,80);
image_haar_B = hard(image_haar_B,80);

invhaar = @(x) (h')'*x*(h');

image_haar_cell_R=mat2cell( double(image_haar_R) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_out_cell_R=cellfun( invhaar, image_haar_cell_R , 'UniformOutput',false);
image_out_R= uint8(cell2mat(image_out_cell_R));

image_haar_cell_G=mat2cell( double(image_haar_G) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_out_cell_G=cellfun( invhaar, image_haar_cell_G , 'UniformOutput',false);
image_out_G= uint8(cell2mat(image_out_cell_G));

image_haar_cell_B=mat2cell( double(image_haar_B) , 8*ones(1,rows/8), 8*ones(1,columns/8) );
image_out_cell_B=cellfun( invhaar, image_haar_cell_B , 'UniformOutput',false);
image_out_B= uint8(cell2mat(image_out_cell_B));

%figure, imshow(uint8(image_out));
image_out = image;
image_out(:,:,1) = image_out_R;
image_out(:,:,2) = image_out_G;
image_out(:,:,3) = image_out_B;
figure,imshow(uint8(image_out));

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
