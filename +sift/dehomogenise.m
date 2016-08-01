function [invariant] = dehomogenise(strain,deltaT,MME,nMatrixIP,nFibreIP,nArrays,nAngles)
%SIFT_user(strain,deltaT,MME,nMatrixIP,nFibreIP,nArrays,nAngles,crit,results
%dehomogenise Summary of this function goes here
%   Detailed explanation goes here
    modStrain = zeros(6,max(nMatrixIP,nFibreIP),nArrays,nAngles);
    for iPosition = 1:nMatrixIP
        for iArray = 1:nArrays
            for iAngle = 1:nAngles
                modStrain(:,iPosition,iArray,iAngle) = MME(:,1:6,iPosition,iArray)*strain
                modStrain(:,iPosition,iArray,iAngle) = modStrain(:,iPosition,iArray,iAngle)+deltaT*MME(:,7,iPosition,iArray)
                matrixDil(:,iPosition,iArray,iAngle) = sift.eDil(modStrain(:,iPosition, iArray, iAngle))
                matrixDis(:,iPosition,iArray,iAngle) = sift.eDis(modStrain(:,iPosition, iArray, iAngle))
            end
        end
    end
    invariant(1) = max(matrixDil);
    invariant(2) = max(matrixDis);
end

