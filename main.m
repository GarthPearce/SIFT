nMatrixIP = 4;
nFibreIP = 3;
nIP = nMatrixIP + nFibreIP;
nArrays = 2;
nAngles = 12;

matrixDilCrit = 0.025;
deltaT = 140;

MME = csvread('MME_T800s_3900-2_0.60.csv');
MME = reshape(MME,[6,7,nIP,nArrays]);

unitStrainTensor = [0.0 1. -0.45 0. 0. 0.]';
increment = 0.0001;
extension = 0;
matrixDil = 0;

while matrixDil < matrixDilCrit
    strain = extension*unitStrainTensor;
    result = sift.dehomogenise(strain,deltaT,MME,nMatrixIP,nFibreIP,nArrays,nAngles);
    matrixDil = result.invariant.matrixDil;
    matrixDis = result.invariant.matrixDis;
    extension = extension + increment;
end

extension
result.invariant.matrixDil
result.invariant.matrixDis
result.invariant.fibreDis