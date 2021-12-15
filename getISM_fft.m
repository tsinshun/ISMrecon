%% Efficient ISM image reconstruction
%  Author: shun qin
%  Contact: shun.qin@outlook.com
%
% input:
% imstack: raw image of real sample, F: index for each frame (from 1 to N),
% (X,Y) is the coordinate of each scanning position
% scale: image scale factor, not necessary but still available to use
% Npixel: the window size (Npixel*Npixel) used to cut each light-spot 
% Bg: background 

% output:
% Iism: reconstructed ISM image with measured image
% Ipat: reconstructed ISM image with the same localization data, but
% use a standard Gaussian function as a light-spot to do the ISM reconstruction
%% 
function [Iism,Ipat] = getISM_fft(imstack, F, X, Y,scale, Npixel,Bg)
% Npixel: Width of PSF

[M,N,Lz] = size(imstack);
% imstack(imstack-Bg<0)=0;


Nframe = max(F(:));
Iism = 0;
% Isup = zeros(M*scale*2,N*scale*2,Lz);
Ipat = 0;
parfor i = 1:Nframe  
    
    if i==37
       1; 
    end

   idx = find(F==i);
%    F1 = F(idx);
   X1 = X(idx);
   Y1 = Y(idx);
   
   frame = imstack(:,:,i);
   frame = frame - Bg; frame(frame<0) = 0; % remove background
   if scale>1
      I = imresize(frame,scale,'Method','lanczos2'); % lanczos2 bicubic
   else
       I = frame;
   end
%    I = frame-frame + 100;
%    I = I/max(I(:));
   [Iism0,Iism1] = copyPSF1(I,Npixel*scale,[X1 Y1]*scale,scale);
%    Isup(:,:,i) = Iism0;
%    Ipat(:,:,i) = Iism1;
   Iism = Iism + Iism0;
   Ipat = Ipat + Iism1;
   fprintf('Frame number = %d \n',i);
   imagesc(Iism0);title(['nFrame: ' num2str(i)]);
   axis image;drawnow
   
end
if scale >1
    Iism = imresize(Iism,2*[M N]);
end
% Iism = sum(Isup,3);
% Ipat = sum(Ipat,3);

% imshow(mat2gray(Iism))
% 
% Iavg = mean(imstack1,3);
% Iavg = imresize(Iavg,2/scale,'Method','lanczos2');


function [Iism, Iism1]= copyPSF1(I,Npixel,Position,Scale)

[M, N] = size(I);
Np = floor(Npixel/2);
% Inew = zeros(M+Np,N+Np);
Iism = zeros(2*M,2*N);
Iism1 = zeros(2*M,2*N);

Inew = padarray(I, [Np Np]);
Iism = padarray(Iism, [Np Np]);  %Iism = Iism + nan;
Iism1 = padarray(Iism1, [Np Np]);

X = Position(:,1);
Y = Position(:,2);
 
L = length(X);

% PSF = gaussian2d(2*Np,3); PSF = PSF./max(max(PSF));%imagesc((PSF));

for i=1:L
    n = round(X(i)+0.5); m = round(Y(i)+0.5);
    x0 = n + Np;  % plus 0.5 to get real pixel index,i.e. (m,n)
    y0 = m + Np; 
    x1 = 2*n + Np;  
    y1 = 2*m + Np;
    V1 = 0.5+[X(i)-n  Y(i)-m];
    
    if x0>0 && y0>0 && x1>0 && y1>0
        if (x1-Np<1 || x1+Np>2*N+2*Np|| y1-Np<1|| y1+Np>2*M+2*Np) || (x0-Np<1 || x0+Np>2*N+2*Np|| y0-Np<1|| y0+Np>2*M+2*Np)
%             disp('window for PSF copy out of idx...');
            continue;
        end
        
%         if y0-Np>0 && y0-Np<=M && x0-Np>0 && x0-Np<=N && G(y0-Np,x0-Np)<=1            
%         end    
         %Iism(y1-Np:y1+Np,x1-Np:x1+Np) =  PSF;
%         W = fspecial('gaussian',2*Np+1,1.09);
%         W = ones(2*Np+1,2*Np+1);

%         [x,y] = meshgrid(x0-Np:x0+Np,y0-Np:y0+Np);
%         sigma = 1.09;
%         r3 = (x-(X(i)+0.5 + Np)).^2 + (y-(Y(i)+0.5 + Np)).^2;
%         W = exp(-r3/(2*sigma^2));

        W = Inew(y0-Np:y0+Np,x0-Np:x0+Np);
        [x,y] = meshgrid(x0-Np:x0+Np,y0-Np:y0+Np);
        r2 = (x-x0).^2 + (y-y0).^2;
        Id = r2>Np^2;
        W(Id)=0;               
                
        W1 = W; 
        W1= ShiftPSF(W,V1); % subpixel shifting by FFT
        Id = Id | (W1<0);     
        W1(Id)=0;     
        Iism(y1-Np:y1+Np,x1-Np:x1+Np) = W1;
        
%         V2 = [x1 y1] - (2*[X(i)+0.5 Y(i)+0.5]+Np);
%         Iism1(y1-Np:y1+Np,x1-Np:x1+Np) =  ShiftPSF(getPSF1(2*Np+1, 1.09),-V2);
        
        [x,y] = meshgrid(x1-Np:x1+Np,y1-Np:y1+Np);
        sigma = 1.09*Scale;
        r3 = (x-(2*(X(i)+0.5) + Np)).^2 + (y-(2*(Y(i)+0.5) + Np)).^2;
        W = exp(-r3/(2*sigma^2));       
        Iism1(y1-Np:y1+Np,x1-Np:x1+Np) = W;
    end

end
Iism = Iism(Np+1:end-Np,Np+1:end-Np);
Iism1 = Iism1(Np+1:end-Np,Np+1:end-Np);


function Ish= ShiftPSF(W,V)

% W1 = padarray(W0,[10 10]);
% W = circshift(W1,[5 10]); 


%%
% V=-[10 5];
[~,~,kx,ky] = getCoordinate(size(W),[1 1]);
Fw = fftshift(fft2(W));
F1 = Fw.*exp(-1i*kx*V(1)-1i*ky*V(2));
Ish = real(ifft2(ifftshift(F1)));
% imagesc(Ish)




function [x1,y1,kx,ky] = getCoordinate(Siz,Delta)
% [M,N] = size(I);
M = Siz(1); N = Siz(2);
dx = Delta(1); dy = Delta(2);
[x0,y0]=meshgrid(1:N,1:M);
x1 = x0-round((N+1)/2);
y1 = y0-round((M+1)/2);
x = x1*dx; % pixel size dx=16um
y = y1*dy;



Lx=max(x(1,:))-min(x(1,:));
Ly=max(y(:,1))-min(y(:,1));
% fx=x1/Lx; fy=y1/Ly;
kx = 2*pi*x1/Lx; ky = 2*pi*y1/Ly; 

