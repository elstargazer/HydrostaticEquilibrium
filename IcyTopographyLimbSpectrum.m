ccc
%% Plotting settings
fntsize = 12;
fntsize_sm = 8;
im_size=[0 0 13 9];

fig_folder='~/Dawn/Papers/CeresPaper1/';



%% Reading data
listname = 'List';
listmaxname = 'MaxList';
listminname = 'MinList';

in = fopen(listname);
in_min = fopen(listminname);
in_max = fopen(listmaxname);

i=1;

filename = fgetl(in);
filename_min = fgetl(in_min);
filename_max = fgetl(in_max);

PlanetName={'Dione','Enceladus','Europa','Iapetus','Mimas','Rhea','Tethys'};

while (filename~=-1)
    
    data=load(filename);
    data_min=load(filename_min);
    data_max=load(filename_max);
    
    [pathstr,namei,ext] = fileparts(filename);
    name{i}=namei;
    
    k{i}=data(:,1);
    k_min{i}=data_min(:,1);
    k_max{i}=data_max(:,1);
    
    sdl{i}=data(:,2);
    sdl_min{i}=data_min(:,2);
    sdl_max{i}=data_max(:,2);
    
    filename = fgetl(in);
    filename_min = fgetl(in_min);
    filename_max = fgetl(in_max);
    
    i=i+1;
    
end

%% Plotting power spectra

figure;
set(gcf, 'Units','centimeters', 'Position',im_size)
set(gcf, 'PaperPositionMode','auto')
set(gca, 'FontSize',fntsize);
hold on;box on;grid on;

xlabel('Frequency [cycles/km]','FontSize',fntsize);
ylabel('Topography power [km^{2}]','FontSize',fntsize);

ccj= jet(numel(sdl));

for i=1:numel(sdl);   
    h(i)=plot(k{i},sdl{i},'Color',ccj(i,:),'LineWidth',2,...
     'MarkerSize',1,'Marker','o');
%     plot(k_min{i},sdl_min{i},'--','Color',ccj(i,:),'LineWidth',1);
%     plot(k_max{i},sdl_max{i},'--','Color',ccj(i,:),'LineWidth',1);
end

legend(h,PlanetName,'FontSize',fntsize);
PrintWhite([fig_folder 'Fig_IcyTopoSpec.jpg']);









