clc;
clear;
close all;

% TODO
% 1. forgot the system mechanics

%% Constants and Other Definitions

mu0 = 4*pi*1e-7; % Vacuum permeability
d = 0.020; % Length of the dipole
w = 0.003; % Dipole width
thick = 0.05; % Dipole Thickness
M = 0.005; % Mass of the dipole
I = (M*d*d)/12; % Moment of inertia calculation 
n=3; %number of magnets

%from datasheet
Br= 1.42; % N52 value
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
energy=zeros(n,1); % clculating energy for each dipole
alpha = zeros(n,1); % Angular acc
w=zeros(n,1);%angular velocity

T = 1000; % Total time of simulation
ts = linspace(0, T, 200); % Time vector -> decrease time step to converge
tracking=zeros(n,2,length(ts)+1); %tracking trajecotry of dipoles

%% Initial values setup

for i=1:n
    disp(["Enter details for ",i,"th dipole"]);
    theta(i,1) = input("Enter initial angle");
    pos(i,1)=input("Enter initial x position");
    pos(i,2)=input("Enter initial y position");
end

%% magnetic field calculayions 
% get data first, do figure later

bmax=0;
tau_max=0;
energy_min=1;

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
        net_energy =0;
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
            % force due to m2 on m1
            B_vec = k*(dot_pord-other)/denom
            
           
            B_vec=[B_vec,0];
            b=norm(B_vec);
            if b>bmax
                bmax=b;
            end

             %total torque from all fields from all magnets
             % torque due to field by m2 on m1
            net_tau = net_tau + cross(m1_vec, B_vec);
            net_energy = net_energy+dot(m1_vec, B_vec)
       end
       if (norm(net_tau)>tau_max)
           tau_max=norm(net_tau); %for self tracking
       end
       
        tau(i,1)=net_tau(3) % net toruq
        energy(i,1)=net_energy; % total energy of the system
        alpha(i,1)=tau(i,1)/I % angular acc
        w(i,1)=w(i,1)+alpha(i,1)*0.01; %angular vel for this time step
        del_theta = w(i,1)*0.01; %change in angle for this time step
        rad = norm([xp, yp]-[0,0]) %distance from origin
        del_moved = del_theta*rad%net angle rotated by
        theta(i,1)=theta(i,1)+del_theta; %new angle for later
        
        rot_mat = [cos(del_theta), -sin(del_theta); 
            sin(del_theta), cos(del_theta)]; % claculating the new position
        new = rot_mat*[xp; yp]
        pos(i,1)=new(1);
        pos(i,2)=new(2);
        tracking(i,:,t)=new;
    end
    sum_e=sum(energy);
    if (sum_e<energy_min)
        disp("new low energy reached");
        disp(["timestep is: ", t]);
        disp(tracking(:,:,t));
        energy_min=sum_e
        input("press enter to continue");
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


x2(:)=tracking(2, 1, :);
y2(:)=tracking(2, 2, :);

x1(:)=tracking(1, 1, :);
y1(:)=tracking(1, 2, :);
x3(:)=tracking(3, 1, :);
y3(:)=tracking(3, 2, :);

plot(x2,y2,'DisplayName','Magnet2'); % Trajectories with random colors
hold on
plot(x1, y1,'DisplayName','Magnet1');
hold on
plot(x3,y3,'DisplayName','Magnet3');
