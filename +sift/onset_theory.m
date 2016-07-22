function [J1max,J2max] = onset_theory(strain_vector,expstrainval,deltaTval)
%Calculate J1max and J2max for a given set of inputs

%Do the coordinate transofrmation to the global strain tensor (if required)
%if coordtrans==1
    %T300_CYCOM970 = strain_ct;
%else
    %T300_CYCOM970 = [-nu12 expstrain (-nu23*expstrain) 0 0 0];
%end
T300_CYCOM970=(strain_vector)*expstrainval;

%import the MME matrices (6 by 7 by 8 matrix) and the rotation matrices
load('MME_T300_CYCOM970.mat');
load('rot_t.mat');


%Multiply strain tensor by each rotation matrix (j=4 rotations) and MME (i=8 layers)
%Add thermal vector multiplied by deltaT(Tcure-Ttest) to each result
for i=1:8
    for j=1:4
        
        %Define local strain tensor based on: global strain tensor, MME, rotations, and deltaT
        MST.T300=STtoMST3(T300_CYCOM970,MME_T300_CYCOM970(:,:,i),deltaTval,rot_t(:,:,j));
        
        %Obtain J1 and J2 from MST
        J1.T300(i,j)=STtoJ1(MST.T300);
        J2.T300(i,j)=STtoJ2_new(MST.T300);
        
    end
end


%Return the maximum J1 and J2 across all rotation matrices and MMEs
  J1array=max(J1.T300);
  J2array=max(J2.T300);
    
  J1max=max(J1array);
  J2max=max(J2array);
  




