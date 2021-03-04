Path = 'F:\GWDG\ISM\data\IHCs_p14_calretinin488_Vglut568_Otof647\red_1';
Fname = 'red_1_MMStack_Pos0.ome.tif';

Img = read_tiff(strcat(Path,'\',Fname),1);

 %%
fileName = 'red_1_MMStack_Pos0.ome.txt'; % localization dat by RapidSTORM

Dat = importdata([Path '\' fileName]);
P = Dat.data;
Fn1 = P(:,3)+1;     % for localization data by RapidSTORM
X1 = P(:,1)/100;    % don't forget divided by the pixel-size set in the localization software
Y1 = P(:,2)/100;

%%

[Iism,Ipat] = getISM_OC(Img, Fn1, X1, Y1, 1, 8, 100);
figure;
subplot(121); imagesc(mean(Img,3));axis image
subplot(122);imagesc(Iism); axis image