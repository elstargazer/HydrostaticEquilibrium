ccc
%% Plotting settings
fntsize = 12;
fntsize_sm = 8;
im_size=[0 0 13 9];


%% Reading data
listname = 'ShapeList';
in = fopen(listname);

i=1;
filename=fgetl(in);
while (filename~=-1)  
    [pathstr,name,ext] = fileparts(filename);
    filecell{i}=name;
    lmcosi_shape{i}=ReadSHNimmo(filename);  
    i=i+1;
    filename=fgetl(in);    
end

fclose(in);

R = [561.4 252.1 198.2 763.5 531.0];

for i=1:numel(lmcosi_shape);
    lmcosi_shape{i}(1,3:4)=lmcosi_shape{i}(1,3:4)*R(i)*1000;
    lmcosi_shape{i}(2:end,3:4)=lmcosi_shape{i}(2:end,3:4)*1000;
end


%% Computing and plotting power spectra
figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

set(gca,'XScale','log');
set(gca,'YScale','log');


ylabel('Frequency [cycles/km]','FontSize',fntsize);
xlabel('Topography power [km^{3}]','FontSize',fntsize);

ccj = jet(numel(filecell));

for i=1:numel(lmcosi_shape);    
    [k_icy,sdl_icy]=PowerSpectrum(lmcosi_shape{i}); 
    plot(k_icy,sdl_icy,'-o','Color',ccj(i,:));    
end




