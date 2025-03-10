clc;
clear;
close all;

% TODO
% 5. convergence
% 2. Animation?
% 3. Hall Sensor Simulink!
% 4. Control on Simulink??
% 6. ENERGY OF THE WHOLE SYSTEM AND TRY TO FIND THE LEAST ENERGY POINT?

%% Constants and Other Definitions

mu0 = 4*pi*1e-7; % Vacuum permeability
d = 0.042; % Length of the dipole
w = 0.008; % Dipole width
thick = 0.01; % Dipole Thickness
M = 0.005; % Mass of the dipole
I = (M*d*d)/12; % Moment of inertia calculation 
n=3; %number of magnets

%from datasheet
Br= 1.2; % N42 value
m_datasheet = Br*(d*w*thick)/mu0;

%% Variables and Matrices Setup

% Intial Values
% positions and angles are with respect to origin aka unmoving but rotating
% first dipole
theta = zeros(n,1); % Angles
pos = zeros(n,2); % Magnet centre Position x,y

% For later calculation
r=zeros(n,n,3); % Magnet Position vectors between 2 magent centres nxn
% 3rd layer is the norm
m=zeros(n,2); % Dipole moment vectors
B = zeros(n,n); % Magnet Field Vectors due to ther dipoles
tau = zeros(n,1); % Torque vectors
alpha = zeros(n,1); % Angular acc

T = 2; % Total time of simulation
ts = linspace(0, T, 200); % Time vector -> decrease time step to converge
tracking=zeros(n,2,length(ts)+1);

%% Initial values setup

for i=1:n
    disp(["Enter details for ",i,"th dipole"]);
    theta(i,1) = input("Enter initial angle");
    pos(i,1)=input("Enter initial x position");
    pos(i,2)=input("Enter initial y position");
end

%% le for loop that killed me : only calcs. 
% get data first, do figure later

bmax=0;
tau_max=0;

for t=1:length(ts)
    tracking(:,:,t)=pos(:,:);

    for k=1:n % for this time step calculate all the dipole's m_vectors
        m(k,1)=m_datasheet*cos(theta(k));
        m(k,2)=m_datasheet*sin(theta(k));
    end
    for i=1:n
        disp(["i=",i]);
        xp=pos(i,1); yp=pos(i,2);
        disp([xp, yp]);
        m1_vec=[m(i,1),m(i,2),0]; %current dipole m to be used in torque calc
        net_tau =0;
       for j=1:n %inner loop to calculate net toruque and field
           if(j==i)
               continue;
           end
            m2_vec=[m(j,1), m(j,2)];
           % r position vector from dipole at j to the i dipole
           r(i,j,1)=xp-pos(j,1); r(i,j,2)=yp-pos(j,2); 
           r_vec = [r(i,j,1),r(i,j,2)];
           disp(r_vec);
           norm(r_vec)
    
            % the field
            k=(mu0 / (4*pi));
            denom = ( norm(r_vec))^5
            dot_pord = 3*dot(m2_vec,r_vec)*r_vec;
            other= m2_vec*(norm(r_vec))^2;
            B_vec = k*(dot_pord-other)/denom
            
           
            B_vec=[B_vec,0];
            b=norm(B_vec);
            if b>bmax
                bmax=b;
            end

             %total torque from all fields from all magnets
            net_tau = net_tau + cross(m1_vec, B_vec);
       end
       if (norm(net_tau)>tau_max)
           tau_max=norm(net_tau);
       end
        tau(i,1)=net_tau(3) % net toruq
        alpha(i,1)=tau(i,1)/I % angular acc
        rad = norm([xp, yp]-[0,0]) %distance from origin
        del_theta = alpha(i,1)*rad%net angle rotated by
        theta(i,1)=theta(i,1)+del_theta; %new angle for later
        
        rot_mat = [cos(del_theta), -sin(del_theta); 
            sin(del_theta), cos(del_theta)];
        new = rot_mat*[xp; yp]
        pos(i,1)=new(1);
        pos(i,2)=new(2);
        tracking(i,:,t)=new;
    end
    
end

%% quick plotting
figure;
hold on;
axis equal;
xlim([-1 1]*0.2); % Adjust axis limits based on your system size
ylim([-1 1]*0.2);
xlabel('X Position (m)');
ylabel('Y Position (m)');


x(:)=tracking(2, 1, :);
y(:)=tracking(2, 2, :);
plot(x,y , 'Color', rand(1, 3)); % Trajectories with random colors

%% Magnetostatic Model

frame=[3;4;[0.2; 0.2; -0.2; -0.2]; 0.2; -0.2; -0.2; 0.2];

%creating the geometry
magnet1=[3;4;[0.021; 0.021; -0.021; -0.021; 
    0.004;   -0.004;  -0.004;  0.004]];

magnet2=[3;4;[0.171; 0.171; 0.129; 0.129; 
    0.004;-0.004; -0.004; 0.004]];

gd=[frame, magnet1, magnet2];
sf="frame+magnet1+magnet2";
ns=char('frame','magnet1','magnet2');
ns=ns';
dl=decsg(gd, sf, ns);

figure;
pdegplot(dl, FaceLabels="on");
%%
model = femodel("AnalysisType","magnetostatic", "Geometry",dl);
model.VacuumPermeability=1.256e-6;
model.MaterialProperties(2)=materialProperties(RelativePermeability=1.05);
model.MaterialProperties(3)=materialProperties(RelativePermeability=1.05);
model.MaterialProperties(1)=materialProperties(RelativePermeability=1);
Mag=1.3;
dir1 = [1;0];
dir2 = [1;0];
model.FaceLoad(2) = faceLoad(Magnetization=Mag*dir1);
model.FaceLoad(3) = faceLoad(Magnetization=Mag*dir2);

model = generateMesh(model, Hface={[2 3], 0.001}, Hmax=0.01, Hgrad=2);
pdemesh(model);
%%
R = solve(model);
Bmag = sqrt(R.MagneticFluxDensity.Bx.^2 + R.MagneticFluxDensity.By.^2);
Hmag = sqrt(R.MagneticField.Hx.^2 + R.MagneticField.Hy.^2);

figure
pdeplot(R.Mesh,XYData=Bmag, ...
               FlowData=[R.MagneticFluxDensity.Bx ...
                         R.MagneticFluxDensity.By])
