## To support project 20160905-mie_AuNPs+shells
## See "supporting information - peak analysis" for original code development
# from header import *

from header_import_data import *

# Shared peakutils settings
threshold = 0.25
wlExpectednm = 525
dwlExpectednm = 35
wlScatnm = 400
x0base = wlExpectednm
nPeaks = 1
wlMin0nm, wlMax0nm = 400, 800
wlMinnm, wlMaxnm = wlExpectednm-dwlExpectednm, wlExpectednm+dwlExpectednm

# Core spectrum, without baseline correction
print('LSPR peak wavelengths')
print('=====================')
print('Bare cores only')
print('  * looking for one peak at', wlExpectednm, 
      'nm \n    - in the range [', wlMinnm, ',', wlMaxnm, '] nm')
print('  * peaks detected: \n    wl-nm (interp),\t [wl-nm,  ext]')

lsprPeakWLnmCoreList = []
for sampleNum in range(2):
    sampleType = (['L-Core_D=','X-Core_D='][sampleNum]) + str(int(round(2*coreRnms[sampleNum]))) + '-nm'
    print('   - Core #',sampleNum,'(',
          sampleType,
          ')')
    subdataArray = spectraCores[sampleNum]
    selection = ((subdataArray[:,0] >= wlMin0nm) 
                 & (subdataArray[:,0] <= wlMax0nm)
                 )
    subdataArray = subdataArray[selection]
    minDist = int(subdataArray.size/(nPeaks+1))
    x, y = subdataArray[:,0], subdataArray[:,1]
    indexes = peakutils.indexes(y, thres=threshold, min_dist=minDist)
    xPeaks = x[indexes]
    anArray = (np.array([xPeaks, indexes])).T
    selection = (anArray[:,0] > wlMinnm) & (anArray[:,0] < wlMaxnm)
    indexes = (anArray[selection][:,1]).tolist()
    xPeaks = x[indexes]
    if (len(indexes) > 0): 
        result = subdataArray[indexes[0]]
        xPeaks = peakutils.interpolate(x, y, indexes)
        print('   ', round(100*xPeaks[0])/100,',\t\t',  
              (np.round(100*result)/100).tolist())
        result = [xPeaks[0], coreRnms[sampleNum]]
        lsprPeakWLnmCoreList.append(result)

print("  * output array \"lsprPeakWLnmCoreList\" format: \n    [ peak-wl-nm, Rcore-nm]\n")

# Core spectrum baseline correction, for later subtraction correction
print('  * power law baseline fit')
print('    - of the form ybase = k*(x/x0)^n')
print('    - using log-log/linear fitting with x0 =',x0base)
print('  * looking for one peak at', wlExpectednm, 
      'nm \n    - in the range [', wlMinnm, ',', wlMaxnm, '] nm')
print('    - baseline correction range [', wlMin0nm, ',', wlMax0nm, '] nm')
print('  * Peaks detected: \n      wl-nm (interp),\t [wl-nm, ext], \t (orig wl-nm),\t gamma Ref')

lsprPowrPeakWLnmCoreList = []
lsprPowrBaselineCoreList = []
lsprPowrBaselineCorenkList = []
lsprPowrGammaCoreRefList = []

for sampleNum in range(2):
    sampleType = (['L-Core_D=','X-Core_D='][sampleNum]) + str(int(round(2*coreRnms[sampleNum]))) + '-nm'
    print('    - Core #',sampleNum,'(', sampleType, ')', end="")
    
    subdataArray = spectraCores[sampleNum]
    selection = ((subdataArray[:,0] >= wlMin0nm) 
                 & (subdataArray[:,0] <= wlMax0nm)
                 )
    subdataArray = subdataArray[selection]
    minDist = int(subdataArray.size/(nPeaks+1))
    
    x, y = subdataArray[:,0], subdataArray[:,1]
    yOrig = y
    logx, logy = np.log(x/x0base), np.log(y)
    dlogx = min(np.diff(logx))
    loglogx = np.arange(min(logx), max(logx) + dlogx, dlogx)
    loglogy = np.interp(loglogx, logx, logy)
    
    loglogyLinBase = peakutils.baseline(loglogy, 1)
    logyLinBase = np.interp(logx, loglogx, loglogyLinBase)
    loglogfit = np.polyfit(logx, logyLinBase, 1)
    basePowrn, baseCoeffk = loglogfit[0], np.exp(loglogfit[1])
    print('\t{baseline n =', '{0:.2f}'.format(basePowrn),', k =', '{0:.1f}'.format(baseCoeffk),'}')
    
    yPowrBase = (np.exp(logyLinBase))
    y = y - yPowrBase
    
    indexes = peakutils.indexes(y, thres=threshold, min_dist=minDist)
    xPeaks = x[indexes]
    anArray = (np.array([xPeaks, indexes])).T
    selection = (anArray[:,0] > wlMinnm) & (anArray[:,0] < wlMaxnm)
    indexes = (anArray[selection][:,1]).tolist()
    xPeaks = x[indexes]
    if (len(indexes) > 0): 
        result = subdataArray[indexes[0]]
        xPeaks = peakutils.interpolate(x, y, indexes)

        xPeak = xPeaks[0]
        xPeakOrig = lsprPeakWLnmCoreList[sampleNum][0]
        select = (np.abs(x-xPeak) < 0.55)
        subdataArrayPeak = subdataArray[select][0]
        select = (np.abs(x-wlScatnm) < 0.55)
        subdataArrayScat = subdataArray[select][0]
        gammaCoreRef = subdataArrayPeak[1]/subdataArrayScat[1]

        print('     ', round(100*xPeaks[0])/100,',\t',  
              (np.round(100*subdataArrayPeak[:])/100).tolist(),
              '\t (',round(100*xPeakOrig)/100,')',
              '\t',(np.round(100*subdataArrayScat[:])/100).tolist(),
              '\t {0:.3f} '.format(gammaCoreRef))
        
        result = [xPeaks[0], lsprPeakWLnmCoreList[sampleNum][0], coreRnms[sampleNum]]
        lsprPowrPeakWLnmCoreList.append(result)

        resultGamma = [result[0], subdataArrayPeak[1], subdataArrayScat[1], 
                       gammaCoreRef, result[-1]]
        lsprPowrGammaCoreRefList.append(resultGamma)

        lsprPowrBaselineCoreList.append([x, yOrig, y[indexes][0], y, yPowrBase])
        lsprPowrBaselineCorenkList.append([x0base, basePowrn, baseCoeffk])

print('  * variable: \t lsprPowrPeakWLnmCoreList')
print("    - output array format: \n    [ peak-wl-nm, orig-peak-wl-nm, Rcore-nm ]")
print('  * variable: \t lsprPowrBaselineCoreList')
print("    - output array format: \n    [ wl-nm-data, ext-data, peak-corr-ext",
      "ext-corrected-data, baseline-correction-data]")
print('    - wavelength range: [',lsprPowrBaselineCoreList[0][0][0],',', 
      lsprPowrBaselineCoreList[0][0][-1],'] x (',len(lsprPowrBaselineCoreList[0][0]), 'values)')
print('  * variable: \t lsprPowrBaselineCorenkList')
print("    - output array format: \n    [ x0, n, k ]")
print('  * variable: \t lsprPowrGammaCoreRefList')
print("    - output array format: \n    [ peak-wl-nm, ext-peak, ext-"+str(wlScatnm)+"-nm, gamma-core-Ref, Rcore-nm ]")


# Core-shell spectrum, without baseline correction
print('\nLSPR peak wavelengths')
print('=====================')
print('Core-shell nanoparticles')
print('  * without baseline correction (estimates for large shells)')
print('  * looking for one peak at', wlExpectednm, 
      'nm \n    - in the range [', wlMinnm, ',', wlMaxnm, '] nm')
print('  * Peaks detected: \n    wl-nm,\t [wl-nm, ext, T-C, Rh-nm], \t (orig wl-nm)',
     '\t [wl-nm, ext, T-C, Rh-nm]', '\t gamma')

lsprPeakWLnmList = []
lsprGammaList = []
x = []

def indexesInRange(x, y):
    indexes = peakutils.indexes(y, thres=threshold, min_dist=minDist)
    xPeaks = x[indexes]
    anArray = (np.array([xPeaks, indexes])).T
    selection = (anArray[:,0] > wlMinnm) & (anArray[:,0] < wlMaxnm)
    indexes = (anArray[selection][:,1]).tolist()
    return indexes


for sampleNum in range(6):
    sampleType = (['L-Thn','L-Med','L-Thk','X-Thn','X-Med','X-Thk'][sampleNum])
    print('   - Sample #',sampleNum,'(', sampleType, ')')
    lsprPeakWLnmSubList = []
    indexT = 0
    coreNum = 0
    if sampleNum >= 3: coreNum = 1
    gammaCoreRef = lsprPowrGammaCoreRefList[coreNum][-2]
    
    for T in spectraTemps[sampleNum]:
        selection = ((dataArray[:,-1] == sampleNum) 
                     & (dataArray[:,-3] == T)
                     & (dataArray[:,0] >= wlMin0nm) 
                     & (dataArray[:,0] <= wlMax0nm)
                     )
        subdataArray = dataArray[selection]
        minDist = int(subdataArray.size/(nPeaks+1))
        x, y = subdataArray[:,0], subdataArray[:,1]
        indexes = peakutils.indexes(y, thres=threshold, min_dist=minDist)
        xPeaks = x[indexes]
        anArray = (np.array([xPeaks, indexes])).T
        selection = (anArray[:,0] > wlMinnm) & (anArray[:,0] < wlMaxnm)
        indexes = (anArray[selection][:,1]).tolist()
        xPeaks = x[indexes]

        if (len(indexes) > 0): 
            result = subdataArray[indexes[0]]
            Rhnm = result[-2]
            xPeaks = peakutils.interpolate(x, y, indexes)
            indexT = indexT + 1
            
            xPeak = xPeaks[0]
            
            select = (np.abs(x-xPeak) < 0.51)
            subdataArrayPeak = subdataArray[select][0]
            select = (np.abs(x-wlScatnm) < 0.51)
            subdataArrayScat = subdataArray[select][0]
            gamma = subdataArrayPeak[1]/subdataArrayScat[1]
            gamma = gamma/gammaCoreRef
            
            print('     ', round(100*xPeaks[0])/100,',\t',  
                  (np.round(100*subdataArrayPeak[0:-1])/100).tolist(),
                  '\t',(np.round(100*subdataArrayScat[0:-1])/100).tolist(),
                  '\t {0:.3f} '.format(gamma))
            
            result = [xPeaks[0], gamma, result[-3], result[-2], int(result[-1])]
            lsprPeakWLnmList.append(result)

            resultGamma = [result[0], subdataArrayPeak[1], subdataArrayScat[1], 
                           gamma, result[-3], result[-2], int(result[-1])]
            lsprGammaList.append(resultGamma)

            result = [xPeaks[0], result[-3], result[-2]]
            lsprPeakWLnmSubList.append(result)
            
        else: # using negative curvature method
            x1 = 0.25*x[:-2] + 0.5*x[1:-1] + 0.25*x[2:]
            x1 = np.concatenate((x1, 0.5*x[:-1] + 0.5*x[1:]))
            y1 = 0.25*y[:-2] + 0.5*y[1:-1] + 0.25*y[2:]
            y1 = np.concatenate((y1, 0.5*y[:-1] + 0.5*y[1:]))
            sort = np.argsort(x1)
            x1, y1 = x1[sort], y1[sort]
            
            dy = np.gradient(y1, np.diff(x1)[0])
            ddy = np.gradient(dy, np.diff(x1)[0])
            select = ((ddy <= 0) & (x1 > 500) & (x1 < 535))
            x2, y2 = x1[select], y1[select]
            xPeaks = [np.mean(x2)]
            indexes = np.arange(len(x))[np.argsort(np.abs(x - xPeaks[0]))][[0]] 
            
            result = subdataArray[indexes[0]]
            Rhnm = result[-2]
            xPeaks = peakutils.interpolate(x, y, indexes)
            indexT = indexT + 1
            
            xPeak = xPeaks[0]
            
            select = np.abs(x-xPeak).argmin()
            subdataArrayPeak = subdataArray[select]
            select = np.abs(x-xPeak).argmin()
            subdataArrayScat = subdataArray[select]
            gamma = subdataArrayPeak[1]/subdataArrayScat[1]
            gamma = gamma/gammaCoreRef
            
            print('     ', round(100*xPeaks[0])/100,',\t',  
                  (np.round(100*subdataArrayPeak[0:-1])/100).tolist(),
                  '\t',(np.round(100*subdataArrayScat[0:-1])/100).tolist(),
                  '\t {0:.3f} '.format(gamma))
            
            result = [xPeaks[0], gamma, result[-3], result[-2], int(result[-1])]
            lsprPeakWLnmList.append(result)

            resultGamma = [result[0], subdataArrayPeak[1], subdataArrayScat[1], 
                           gamma, result[-3], result[-2], int(result[-1])]
            lsprGammaList.append(resultGamma)

            result = [xPeaks[0], result[-3], result[-2]]
            lsprPeakWLnmSubList.append(result)

            
            
lsprPeakWLnmArray = np.asarray(lsprPeakWLnmList)
lsprScatGammaArray = np.asarray(lsprGammaList)

print('  * variable: \t lsprPeakWLnmArray')
print("    - output array format: \n      [ peak-wl-nm, gamma, temp-C, Rh-nm, sample-# ]")
print("      sample-# = [0,...,5] \n      for [lin-Thin, lin-Med., lin-Long, X-Thin, X-Med., X-Long]")
print('  * variable: \t lsprScatGammaArray')
print("    - output array format: \n      [ peak-wl-nm, ext-peak, ext-"+str(wlScatnm)+"-nm, gamma, temp-C, Rh-nm, sample-# ]")
print("      sample-# = [0,...,5] \n      for [lin-Thin, lin-Med., lin-Long, X-Thin, X-Med., X-Long]")


# Core-shell spectrum, with baseline power-law correction
print('LSPR peak wavelengths')
print('=====================')
print('Core-shell nanoparticles')
print('  * power law baseline fit')
print('    - of the form ybase = k*(x/x0)^n')
print('    - using log-log/linear fitting with x0 =',x0base)
print('  * looking for one peak at', wlExpectednm, 
      'nm \n    - in the range [', wlMinnm, ',', wlMaxnm, '] nm')
print('    - baseline correction range [', wlMin0nm, ',', wlMax0nm, '] nm')
print('  * Peaks detected: \n    wl-nm,\t [wl-nm, ext, T-C, Rh-nm], \t (orig wl-nm)',
     '\t [wl-nm, ext, T-C, Rh-nm]', '\t gamma')

lsprPowrPeakWLnmList = []
lsprPowrBaselineList = []
lsprPowrBaselinenkList = []
lsprPowrBaselineExpvsScaPrmList = []
lsprPowrGammaList = []
x = []

def indexesInRange(x, y):
    indexes = peakutils.indexes(y, thres=threshold, min_dist=minDist)
    xPeaks = x[indexes]
    anArray = (np.array([xPeaks, indexes])).T
    selection = (anArray[:,0] > wlMinnm) & (anArray[:,0] < wlMaxnm)
    indexes = (anArray[selection][:,1]).tolist()
    return indexes


for sampleNum in range(6):
    sampleType = (['L-Thn','L-Med','L-Thk','X-Thn','X-Med','X-Thk'][sampleNum])
    print('   - Sample #',sampleNum,'(', sampleType, ')')
    lsprPeakWLnmSubList = []
    indexT = 0
    coreNum = 0
    if sampleNum >= 3: coreNum = 1
    gammaCoreRef = lsprPowrGammaCoreRefList[coreNum][-2]
    
    for T in spectraTemps[sampleNum]:
        selection = ((dataArray[:,-1] == sampleNum) 
                     & (dataArray[:,-3] == T)
                     & (dataArray[:,0] >= wlMin0nm) 
                     & (dataArray[:,0] <= wlMax0nm)
                     )
        subdataArray = dataArray[selection]
        minDist = int(subdataArray.size/(nPeaks+1))
        
        x, y = subdataArray[:,0], subdataArray[:,1]
        yOrig = y
        logx, logy = np.log(x/x0base), np.log(y)
        dlogx = min(np.diff(logx))
        loglogx = np.arange(min(logx), max(logx) + dlogx, dlogx)
        loglogy = np.interp(loglogx, logx, logy - min(logy))

        loglogyLinBase = peakutils.baseline(loglogy, 1)
        logyLinBase = np.interp(logx, loglogx, loglogyLinBase) + min(logy)
        loglogfit = np.polyfit(logx, logyLinBase, 1)
        basePowrn0, baseCoeffk0 = loglogfit[0], np.exp(loglogfit[1])
    
        yPowrBase = (np.exp(logyLinBase))
        y = y - yPowrBase
        indexes = indexesInRange(x, y)
        # if no peaks:
        if (len(indexes) == 0): 
            selection = (x >= wlMinnm) & (x <= wlMaxnm)
            xzoom, yzoom = x[selection], y[selection]
            indexeszoom = peakutils.indexes(yzoom, thres=threshold, min_dist=minDist)
            indexes = indexeszoom + (wlMinnm - wlMin0nm)
        xPeaks = x[indexes]
                
        selection = ( (lsprPowrBaselineCoreList[coreNum][0] >= wlMin0nm)
                     & (lsprPowrBaselineCoreList[coreNum][0] <= wlMax0nm)
                     )
        yCorePowrBase = np.asarray(lsprPowrBaselineCoreList[coreNum][-1])[selection]
        yCorePeak = lsprPowrBaselineCoreList[coreNum][2]
        scaling = 1*y[indexes][0]/yCorePeak
        y = y + yCorePowrBase*scaling
        
        # fit [(base) - (core base)] to power law
        yPowrBaseTotal = yPowrBase - yCorePowrBase*scaling
        logy = np.log(yPowrBaseTotal)
        if min(yPowrBaseTotal) > 0: 
            logy = np.log(yPowrBaseTotal)
            loglogfit = np.polyfit(logx, logy, 1)
        else: loglogfit = [0, 0]
        basePowrn, baseCoeffk = loglogfit[0], np.exp(loglogfit[1])
        
        indexes = indexesInRange(x, y)
        # if no peaks:
        if (len(indexes) == 0): 
            selection = (x >= wlMinnm) & (x <= wlMaxnm)
            xzoom, yzoom = x[selection], y[selection]
            indexeszoom = peakutils.indexes(yzoom, thres=threshold, min_dist=minDist)
            indexes = indexeszoom + (wlMinnm - wlMin0nm)
        xPeaks = x[indexes]

        # iterate 
        for i in range(20):
            oldScaling = scaling
            scaling = 1*(y[indexes[0]] - yCorePowrBase[indexes[0]]*scaling)/yCorePeak
            if ( abs(1-oldScaling/scaling) < 0.000001 ): break
            
            y = yOrig - yPowrBase + yCorePowrBase*scaling
            indexes = indexesInRange(x, y)
            # if no peaks:
            if (len(indexes) == 0): 
                selection = (x >= wlMinnm) & (x <= wlMaxnm)
                xzoom, yzoom = x[selection], y[selection]
                indexeszoom = peakutils.indexes(yzoom, thres=threshold, min_dist=minDist)
                indexes = indexeszoom + (wlMinnm - wlMin0nm)

        if (len(indexes) > 0): 
            result = subdataArray[indexes[0]]
            Rhnm = result[-2]
            xPeaks = peakutils.interpolate(x, y, indexes)
            indexT = indexT + 1
            
            xPeak = xPeaks[0]
            
            select = (np.abs(x-xPeak) < 0.51)
            subdataArrayPeak = subdataArray[select][0]
            select = (np.abs(x-wlScatnm) < 0.51)
            subdataArrayScat = subdataArray[select][0]
            gamma = subdataArrayPeak[1]/subdataArrayScat[1]
            gamma = gamma/gammaCoreRef
            
            print('     ', round(100*xPeaks[0])/100,',\t',  
                  (np.round(100*subdataArrayPeak[0:-1])/100).tolist(),
                  '\t',(np.round(100*subdataArrayScat[0:-1])/100).tolist(),
                  '\t {0:.3f} '.format(gamma))
            
            result = [xPeaks[0], gamma, result[-3], result[-2], int(result[-1])]
            lsprPowrPeakWLnmList.append(result)

            resultGamma = [result[0], subdataArrayPeak[1], subdataArrayScat[1], 
                           gamma, result[-3], result[-2], int(result[-1])]
            lsprPowrGammaList.append(resultGamma)

            lsprPowrBaselineList.append([x, yOrig, scaling, y, yPowrBase, 
                                         yPowrBase-yCorePowrBase*scaling, Rhnm, sampleNum])
            lsprPowrBaselinenkList.append([x0base, basePowrn0, baseCoeffk0])
            result = [xPeaks[0], result[-3], result[-2]]
            lsprPeakWLnmSubList.append(result)
            
            scatParam = 2*pi*Rhnm/(xPeaks[0]/1.35)
            lsprPowrBaselineExpvsScaPrmList.append([scatParam, basePowrn])
            
lsprPowrPeakWLnmArray = np.asarray(lsprPowrPeakWLnmList)
lsprPowrScatGammaArray = np.asarray(lsprPowrGammaList)

print('  * variable: \t lsprPowrPeakWLnmArray')
print("    - output array format: \n      [ peak-wl-nm, gamma, temp-C, Rh-nm, sample-# ]")
print("      sample-# = [0,...,5] \n      for [lin-Thin, lin-Med., lin-Long, X-Thin, X-Med., X-Long]")
print('  * variable: \t lsprPowrBaselineList')
print("    - output array format: \n",
      "      [ wl-nm-data, ext-orig-data, ext-scaling-factor, ext-corrected-data,",
      "baseline-correction-no-core-data, baseline-correction-data, Rh-nm, sample-#]")
print('    - wavelength range: [',x[0],',',x[-1],'] x (',len(x), 'values)')
print('  * variable: \t lsprPowrBaselinenkList')
print("    - output array format: \n      [ x0, n0, k0 ]")
print('  * variable: \t lsprPowrBaselineExpvsScaPrmList')
print("    - output array format: \n      [ ScatParam, n-exponent ]")
print('  * variable: \t lsprPowrScatGammaArray')
print("    - output array format: \n      [ peak-wl-nm, ext-peak, ext-"+str(wlScatnm)+"-nm, gamma, temp-C, Rh-nm, sample-# ]")
print("      sample-# = [0,...,5] \n      for [lin-Thin, lin-Med., lin-Long, X-Thin, X-Med., X-Long]")
