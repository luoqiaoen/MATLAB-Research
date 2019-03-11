clear
tic

DimArrayX= (2/0.2); % 8um/200nm
DimArrayY= (130/0.2); %130um/200nm
% scatt1 = -11.753 -1.2596i;
scatt1 = 2; % refractive index = 2


% SingleArrangement = [9,32,35,36,37,39,48,49,51];

% MeshSizeUp = 'wl/';
% MeshSizeDown = num2str(50);
% MeshMaxSize = strcat(MeshSizeUp,MeshSizeDown);



% function []=ISStatisticalAnalysisRandomNanostructure(DimArrayY, scatt1, SingleArrangement, MeshMaxSize)
% function []=ISStatisticalAnalysisRandomNanostructure(DimArrayY, scatt1)

%% start simulation
%fprintf(1, 'Start Simulation\n');
import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/home/min/a/luo132/COMSOL/');

model.modelNode.create('mod1');

model.geom.create('geom1', 2);

model.mesh.create('mesh1', 'geom1');

model.physics.create('emw', 'ElectromagneticWaves', 'geom1');

model.study.create('std1');
model.study('std1').feature.create('freq', 'Frequency');
model.study('std1').feature('freq').activate('emw', true);

%% set length unit
model.geom('geom1').lengthUnit([native2unicode(hex2dec('00b5'), 'Cp1252') 'm']);

%% set parameters
model.param.set('wl', '600e-9 [m]');
model.param.descr('wl', 'wavelength');
model.param.set('c', '3e8 [m/s]');
model.param.descr('c', 'speed of light');
model.param.set('f', 'c/wl'); %3e14
model.param.descr('f', 'frequency');
model.param.set('rr', '0.025');
model.param.descr('rr', 'radius of single scatterer');
model.param.set('theta_i', '0');
model.param.descr('theta_i', 'angle of incidence');
model.param.set('eps0', '8.854187817*1e-12 [F/m]');
model.param.set('mu0', '4*pi*1e-7 [N/A^2]');
model.param.set('eta0','sqrt(mu0/eps0)');
model.param.set('W_slab','130 [um]');
model.param.set('l_th','2 [um]');
model.param.set('p_side','200 [nm]')

%% create geometry
lambda=0.6; % wl = red light
sizeNum = 0.2;
model.geom('geom1').feature.create('r1', 'Rectangle'); %free space
model.geom('geom1').feature('r1').set('type', 'solid');
model.geom('geom1').feature('r1').set('base', 'corner');
model.geom('geom1').feature('r1').set('pos', {'-5*wl' '-5*wl'});
% model.geom('geom1').feature('r1').set('size', {'l_th + 10*wl' 'W_slab + 10*wl'});
model.geom('geom1').feature('r1').set('size', {'l_th + 10*wl' 'W_slab'});
model.geom('geom1').feature.create('r2', 'Rectangle'); % PML
model.geom('geom1').feature('r2').set('type', 'solid');
model.geom('geom1').feature('r2').set('base', 'corner');
model.geom('geom1').feature('r2').set('pos', {'l_th + 5*wl' '-6*wl'});
% model.geom('geom1').feature('r2').set('size', {'wl' 'W_slab + 12*wl'});
model.geom('geom1').feature('r2').set('size', {'wl' 'W_slab'});
model.geom('geom1').feature.create('r3', 'Rectangle'); % PML
model.geom('geom1').feature('r3').set('type', 'solid');
model.geom('geom1').feature('r3').set('base', 'corner');
model.geom('geom1').feature('r3').set('pos', {'-6*wl' '-6*wl'});
% model.geom('geom1').feature('r3').set('size', {'wl' 'W_slab + 12*wl'});
model.geom('geom1').feature('r3').set('size', {'wl' 'W_slab'});
model.geom('geom1').feature.create('r4', 'Rectangle'); % PML
model.geom('geom1').feature('r4').set('type', 'solid');
model.geom('geom1').feature('r4').set('base', 'corner');
% model.geom('geom1').feature('r4').set('pos', {'-6*wl' 'W_slab +5*wl'});
model.geom('geom1').feature('r4').set('pos', {'-6*wl' 'W_slab'});
model.geom('geom1').feature('r4').set('size', {'l_th +12*wl' 'wl'});
model.geom('geom1').feature.create('r5', 'Rectangle'); % PML
model.geom('geom1').feature('r5').set('type', 'solid');
model.geom('geom1').feature('r5').set('base', 'corner');
model.geom('geom1').feature('r5').set('pos', {'-6*wl' '-6*wl'});
model.geom('geom1').feature('r5').set('size', {'l_th +12*wl' 'wl'});


% model.geom('geom1').feature.create('sq5', 'Square');
% model.geom('geom1').feature('sq5').set('type', 'solid');
% model.geom('geom1').feature('sq5').set('base', 'corner');
% model.geom('geom1').feature('sq5').set('pos', {'-1' '-1'});
% model.geom('geom1').feature('sq5').set('size', '0.2');
% model.geom('geom1').run('sq5');

%% Parameters

p=0.5; % filling factor
q=1-p;
StartPointComsol=7; % geometry labels
NumberOfZeros=DimArrayX*DimArrayY*q;
NumberOfOnes=DimArrayX*DimArrayY*p;
R=[repmat(0,1,ceil(NumberOfZeros)) repmat(1,1,floor(NumberOfOnes))];
R=R(randperm(length(R)));
a=linspace(3,DimArrayX*DimArrayY+2,DimArrayX*DimArrayY);
MaterialsWithScatterAsMaterial=a.*R;
%nameMaterialsWithScatterAsMaterial=[int2str(MaterialsWithScatterAsMaterial)];
MaterialsWithScatterAsMaterialWithoutZero=zeros;
AnzahlInNeuemVektor=1;
%% Create Irregular Structure
% R_design = [0 0 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 1 0 0 1 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1];
R_design = R;
R_num = 1;
% FilRadius = num2str(sizeNum/10);
sqtags = cell(1,DimArrayX*DimArrayY);
for j=1:1:DimArrayX
    
    for i=1:1:DimArrayY
        if R_design(R_num) == 1
            name=['sq' int2str(j*100000) int2str(i)];
    
            posx=[num2str((j-1)*sizeNum + sizeNum/2)];
            posy=[num2str((i-1)*sizeNum + sizeNum/2)];
    
            model.geom('geom1').feature.create(name, 'Square');
            model.geom('geom1').feature(name).set('createselection', 'on');
            model.geom('geom1').feature(name).set('type', 'solid');
            model.geom('geom1').feature(name).set('base', 'center'); % sq5 or name
            model.geom('geom1').feature(name).set('pos', {posx posy});
            model.geom('geom1').feature(name).set('size', 'p_side');
            
%             model.geom('geom1').feature.create('scatterer', 'Selection');
%             filname = ['fil' int2str(j*1000) int2str(i)];
%     
%             model.geom('geom1').create(filname, 'Fillet');
%             model.geom('geom1').feature(filname).selection('point').set(name, [1 2 3 4]);
%             model.geom('geom1').feature(filname).set('radius', FilRadius);
            sqtags{R_num} = name;
        end
        R_num = R_num + 1;
    end
end

model.geom('geom1').runAll;
model.geom('geom1').run;
sqtags_filled = cell(1,DimArrayX*DimArrayY/2);
for m=1:1:DimArrayX*DimArrayY
    if ~isempty(sqtags{m})
        sqtags_filled{AnzahlInNeuemVektor}=sqtags{m};
        AnzahlInNeuemVektor=AnzahlInNeuemVektor+1;
    end
end

for i = 1: 1: length(sqtags_filled)
    composite_tag = strcat('geom1_',sqtags_filled{i},'_dom');
    temp = mphgetselection(model.selection(composite_tag));
    sqtags_array(i) = temp.entities;
end

%% create material
%fprintf(1, 'Adding Material\n'); 
model.material.create('mat1');
model.material('mat1').name('Air');
model.material('mat1').set('family', 'air');
model.material('mat1').propertyGroup('def').set('relpermeability', '1');
model.material('mat1').propertyGroup('def').set('relpermittivity', '1');
model.material('mat1').propertyGroup('def').set('electricconductivity', '0[S/m]');
model.material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.material('mat1').propertyGroup('RefractiveIndex').set('n', '1');
model.material('mat1').set('family', 'air');
model.material.create('mat2');
model.material('mat2').name('Silicon');
model.material('mat2').set('family', 'custom');
model.material('mat2').propertyGroup('def').set('relpermeability', '1');
model.material('mat2').propertyGroup('def').set('electricconductivity', '0e-12[S/m]');
%model.material('mat2').propertyGroup('def').set('relpermittivity', '12.11');
%model.material('mat2').propertyGroup('def').set('relpermittivity', '1.156');
model.material('mat2').propertyGroup('def').set('relpermittivity', num2str(scatt1));
%model.material('mat2').propertyGroup.create('RefractiveIndex', 'Refractive index');
%model.material('mat2').propertyGroup('RefractiveIndex').set('n', '3.48');
%model.material('mat2').propertyGroup('RefractiveIndex').set('n', '1.075');

for m=1:1:DimArrayX*DimArrayY
    if MaterialsWithScatterAsMaterial(m) ~= 0
        MaterialsWithScatterAsMaterialWithoutZero(AnzahlInNeuemVektor)=MaterialsWithScatterAsMaterial(m);
        AnzahlInNeuemVektor=AnzahlInNeuemVektor+1;
    end
end



%model.material('mat2').selection.set([14 15 17 19 21 22 23 24 27 28 31 33 35 40 47 48 50 51 52 53 56 58 59 61 62 65 66 67 69 70 71 74 76 77 79 80 81 83 87 88 89 90 91 92 94 95 96 97 98 99 101 103 104 105]);
model.material('mat2').selection.set(sqtags_array);
% SingleArrangement = linspace(7, floor(NumberOfOnes)+6, floor(NumberOfOnes));
% 
% model.material('mat2').selection.set(SingleArrangement*2-13);

%% fillet the corners

%% Perfectly Matched Layer
model.coordSystem.create('pml1', 'geom1', 'PML');
model.coordSystem('pml1').selection.set([1 2 3 4 6 (floor(NumberOfOnes)+7) (floor(NumberOfOnes)+8) (floor(NumberOfOnes)+9)]);

%% Physics
model.physics('emw').prop('BackgroundField').set('SolveFor', 'scatteredField');
model.physics('emw').prop('components').set('components', 'inplane');
model.physics('emw').prop('BackgroundField').set('Eb', {'0' 'exp(-i*x*cos(theta_i)*2*pi*f/c + -i*y*sin(theta_i)*2*pi*f/c)' '0'});
% model.physics('emw').prop('BackgroundField').set('Eb', {'0' 'exp(-i*x*cos(theta_i)*2*pi*f/c + -i*y*sin(theta_i)*2*pi*f/c)*sqrt(2*eta0)' '0'});
%% Mesh

%% Set up study

%% Create Plot

%% Postprocessing

%% Save the Model
%fprintf(1, 'Saving\n');

model.save('C:\Users\luoqi\Dropbox\Workspace\COMSOL\thickness_2.mph');
%model.name('Test1.mph');

out = model;


toc