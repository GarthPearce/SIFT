function MST=STtoMST3(ST,MME,deltaT,rot)

ST3=rot*ST'; %transpose strain tensor and multiply with the current rotation
% sizing=size(ST3);
MME_dT=zeros(size(ST3)); %allocate array matrix

for i = 1:size(ST3,2)
    MME_dT(:,i)=MME(:,7); %to set the A
end
MST=MME(:,1:6)*ST3+deltaT*MME_dT;
% MST=MME(:,1:6)*ST3;