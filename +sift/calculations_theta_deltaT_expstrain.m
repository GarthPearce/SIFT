%Find critical invariant J (J1 or J2)

%Define all the model inputs
%coordtrans = 0;
%nu12=0.25;
%nu23 = 0.45;
%a = 0;
%d = 0;
%e = 0;
%f = 0;
alpha = 2.87e-5;


%Inputs and Parameters
expstrain_step = 0.00001;
J1crit = 0.0316;
J2crit = 0.1368;


%Fibre angle values to test
Theta_min=5;
Theta_max=90;
Theta_step=5;
%Theta_total=floor((Theta_max-Theta_min)/Theta_step + 1);

%deltaT values to test
deltaT_min=-200;
deltaT_max=200;
deltaT_step=50;
%deltaT_total=floor((deltaT_max-deltaT_min)/(deltaT_step)+1);


%i is which combination of Theta and deltaT is being run (iteration at each
%point)
i=0;

%results table 
% results_header(1,1)='iteration';
% results_header(1,2)='Theta';
% results_header(1,3)='deltaT';
% results_header(1,4)='J1crit_expstrain';
% results_header(1,5)='J1crit_J1';
% results_header(1,6)='J1crit_J2';
% results_header(1,7)='J2crit_expstrain';
% results_header(1,8)='J2crit_J1';
% results_header(1,9)='J2crit_J2';


for x=Theta_min:Theta_step:Theta_max %fibre angle
    for y=deltaT_min:deltaT_step:deltaT_max %temp

        Theta=x;
        deltaT=y;
        i=i+1
        
                if Theta==5
                theta_strain= [3.40E+00	-4.17E-01	-1.30E+00	-9.59E+00	4.24E-07	1.94E-06];
                elseif Theta==10
                theta_strain=[8.53E+00	2.31E+00	-4.79E+00	-4.27E+01	-6.36E-07	9.21E-07];
                elseif Theta==15
                theta_strain=[1.20E+01	1.16E+01	-1.05E+01	-8.98E+01	-7.78E-07	1.60E-07];
                elseif Theta==20
                theta_strain=[1.35E+01	2.77E+01	-1.86E+01	-1.39E+02	8.03E-07	6.80E-07];
                elseif Theta==25
                theta_strain=[1.37E+01	4.99E+01	-2.88E+01	-1.85E+02	-2.25E-06	8.92E-07];
                elseif Theta==30
                theta_strain=[1.30E+01	7.78E+01	-4.13E+01	-2.27E+02	-5.36E-07	-1.10E-06];
                elseif Theta==35
                theta_strain=[1.18E+01	1.11E+02	-5.59E+01	-2.62E+02	1.51E-07	4.94E-07];
                elseif Theta==40
                theta_strain=[1.01E+01	1.49E+02	-7.25E+01	-2.90E+02	2.99E-07	2.37E-07];
                elseif Theta==45
                theta_strain=[8.00E+00	1.91E+02	-9.07E+01	-3.10E+02	2.56E-07	-6.17E-07];
                elseif Theta==50
                theta_strain=[5.68E+00	2.36E+02	-1.10E+02	-3.20E+02	7.73E-07	-1.15E-06];
                elseif Theta==55
                theta_strain=[3.19E+00	2.83E+02	-1.31E+02	-3.19E+02	-9.72E-07	8.75E-07];
                elseif Theta==60
                theta_strain=[6.35E-01	3.30E+02	-1.51E+02	-3.06E+02	-9.32E-07	-2.03E-07];
                elseif Theta==65
                theta_strain=[-1.86E+00	3.76E+02	-1.71E+02	-2.81E+02	-4.81E-07	4.27E-07];
                elseif Theta==70
                theta_strain=[-4.17E+00	4.19E+02	-1.90E+02	-2.44E+02	5.07E-07	7.30E-07];
                elseif Theta==75
                theta_strain=[-6.16E+00	4.55E+02	-2.06E+02	-1.95E+02	6.21E-07	1.40E-07];
                elseif Theta==80
                theta_strain=[-7.70E+00	4.83E+02	-2.18E+02	-1.36E+02	-3.22E-07	2.21E-06];
                elseif Theta==85
                theta_strain=[-8.67E+00	5.01E+02	-2.25E+02	-6.99E+01	-1.07E-06	6.51E-07];
                elseif Theta==90 
                theta_strain=[-9.01E+00	5.07E+02	-2.28E+02	1.56E-06	3.77E-07	3.61E-07];
                else
                theta_strain='error';
            end
        
        results(i,1)=i;
        results(i,2)=Theta;
        results(i,3)=deltaT;
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Find expstrain that produces J1 crticial invariant
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %reset initial expstrain
        expstrain = 0;

        %While loop check and counter
        J1crit_invariant = 0;
        J1count = 0;

        while J1crit_invariant==0

            %define what expstrain to use
           % [strain_vector] = ansys_strains(Theta,theta_strain);
            expstrain = expstrain_step*J1count;

            %calculate J1max and J2max
            [J1max,J2max] = onset_theory(theta_strain,expstrain,deltaT);

            %Check if J1 critical invariant has been reached
            if J1max >= J1crit
                J1crit_invariant = 1;
            else
                J1crit_invariant = 0;
            end

            %Increase counter
            J1count = J1count + 1;

        end

        results(i,4)=expstrain;
        results(i,5)=J1max;
        results(i,6)=J2max;

        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Find expstrain that produces J2 crticial invariant
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %reset initial expstrain
        expstrain = 0;

        %While loop check and counter
        J2crit_invariant = 0;
        J2count = 0;

        while J2crit_invariant==0

            %define what expstrain to use
            %[theta_strain] = ansys_strains(Theta,theta_strain);
            expstrain = expstrain_step*J2count;

            %calculate J1max and J2max
            [J1max,J2max] = onset_theory(theta_strain,expstrain,deltaT);

            %Check if J1 critical invariant has been reached
            if J2max >= J2crit
                J2crit_invariant = 1;
            else
                J2crit_invariant = 0;
            end

            %Increase counter
            J2count = J2count + 1;

        end


        results(i,7)=expstrain;
        results(i,8)=J1max;
        results(i,9)=J2max;

        
    end
        
        
        end

   %%%%%%%%%%%%%%
%% PLOTTING %%
%%%%%%%%%%%%%%
            figure
            plot(results(:,5),results(:,6));
            hold on
            plot(results(:,8),results(:,9));
            title('J1 vs J2')
            xlabel('J1')
            ylabel('J2')

            
            figure
            scatter3((results(:,2)),(results(:,3)),(results(:,4)));
            hold on
            scatter3(results(:,2),results(:,3),results(:,7));        
            title('Graph of experimental strain, DeltaT and Theta')
            xlabel('Theta')
            ylabel('DeltaT')
            zlabel('Experimental Strain')


            figure
            tri = delaunay(results(:,2),results(:,3));
            plot(results(:,2),results(:,3),'.');
            h=trisurf(tri,results(:,2),results(:,3),results(:,4));
            hold on
            k=trisurf(tri,results(:,2),results(:,3),results(:,7));
            axis vis3d
            lighting phong
            colorbar EastOutside
            title('Surface plot of experimental strain, DeltaT and Theta')
            xlabel('Theta')
            ylabel('DeltaT')
            zlabel('Experimental Strain')

                
             
  
      





