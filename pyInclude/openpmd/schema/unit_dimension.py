from __future__ import annotations

from types import SimpleNamespace


DIMENSIONLESS = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
LENGTH = (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
AREA = (2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
CROSS_SECTION = AREA
TIME = (0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0)
INV_LENGTH = (-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
INV_VOLUME = (-3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
PHOTON_FLUX = (-2.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0)
RATE = (0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0)

unitDimension = SimpleNamespace(
    dimensionless=DIMENSIONLESS,
    length=LENGTH,
    wavelength=LENGTH,
    lambda_=LENGTH,
    area=AREA,
    crossSection=CROSS_SECTION,
    time=TIME,
    tFluo=TIME,
    inverseLength=INV_LENGTH,
    inverseVolume=INV_VOLUME,
    photonFlux=PHOTON_FLUX,
    rate=RATE,

    numberOfPoints=DIMENSIONLESS,
    numberOfTriangles=DIMENSIONLESS,
    numberOfLevels=DIMENSIONLESS,
    thickness=LENGTH,
    nTot=INV_VOLUME,
    crystalTFluo=TIME,
    claddingNumber=DIMENSIONLESS,
    claddingAbsorption=INV_LENGTH,
    minRaysPerSample=DIMENSIONLESS,
    maxRaysPerSample=DIMENSIONLESS,
    mseThreshold=DIMENSIONLESS,
    repetitions=DIMENSIONLESS,
    adaptiveSteps=DIMENSIONLESS,
    useReflections=DIMENSIONLESS,
    spectralResolution=DIMENSIONLESS,
    lambdaResolution=DIMENSIONLESS,
    monochromatic=DIMENSIONLESS,
    maxSigmaAbsorption=CROSS_SECTION,
    maxSigmaEmission=CROSS_SECTION,
    backend=DIMENSIONLESS,
    maxGpus=DIMENSIONLESS,
    parallelMode=DIMENSIONLESS,
    minSampleRange=DIMENSIONLESS,
    maxSampleRange=DIMENSIONLESS,
    rngSeed=DIMENSIONLESS,

    position=LENGTH,
    points=LENGTH,
    connectivity=DIMENSIONLESS,
    neighbors=DIMENSIONLESS,
    forbiddenEdges=DIMENSIONLESS,
    normalPoints=DIMENSIONLESS,
    center=LENGTH,
    normal=DIMENSIONLESS,
    surface=AREA,
    claddingGroup=DIMENSIONLESS,
    refractiveIndex=DIMENSIONLESS,
    reflectivity=DIMENSIONLESS,
    betaVolume=DIMENSIONLESS,
    pointBeta=DIMENSIONLESS,
    phiAse=PHOTON_FLUX,
    mse=DIMENSIONLESS,
    totalRays=DIMENSIONLESS,
    dndtAse=RATE,
    cellCenterX=LENGTH,
    cellCenterY=LENGTH,
    cellNormalX=DIMENSIONLESS,
    cellNormalY=DIMENSIONLESS,
    claddingCellType=DIMENSIONLESS,
    lambdaAbsorption=LENGTH,
    lambdaEmission=LENGTH,
    sigmaAbsorption=CROSS_SECTION,
    sigmaEmission=CROSS_SECTION,

    canonicalPoints=LENGTH,
    canonicalConnectivity=DIMENSIONLESS,
    canonicalOffsets=DIMENSIONLESS,
    canonicalCellTypes=DIMENSIONLESS,
    corePoints=LENGTH,
    coreCellsConnectivity=DIMENSIONLESS,
    coreCellsOffsets=DIMENSIONLESS,
    coreCellsTypes=DIMENSIONLESS,
)
