## To support project 20160905-mie_AuNPs+shells
## See "supporting information - data analysis" for original code development
# from header import *

from header import *

# Linear PNIPAM shells
print("Import and organize data")
print("========================")
print("Linear PNIPAM: ")
folderName = 'data_20160905'
fileNames = []
print("Importing data in folder",folderName,"...")
for file in os.listdir(folderName):
    #if file.endswith('.dat'): 
        fileNames.append(file)
spectraNames = []
spectraTemps = []
fileNames = np.sort(fileNames)

def importData(fileName):
    print("  *",fileName)
    headerLines = 0
    if (fileName == 'Au-linPNIPAM-long-hydr-radius.dat'
        or fileName == 'Au-linPNIPAM-medium-hydr-radius.dat'
        or fileName == 'Au-linPNIPAM-short-hydr-radius.dat'): 
        headerLines = 3
    elif (fileName == 'Au-linPNIPAM-long.dat'
          or fileName == 'Au-linPNIPAM-medium.dat'
          or fileName == 'Au-linPNIPAM-short.dat'): 
        headerLines = 2
        columnNames = np.genfromtxt(folderName+'/'+fileName, delimiter='\t', dtype=str, comments='\\', skip_header=0, max_rows=1)[1:-1]
        spectraNames.append(columnNames.tolist())
        temps = []
        for name in columnNames:
            temps.append(float(name[-5:-3]))
        spectraTemps.append(temps)
        #print(columnNames)
    else: headerLines = 2
    tempData = np.genfromtxt(folderName+'/'+fileName, delimiter='\t', skip_header=headerLines)
    return tempData
data = list(map(importData,fileNames))

# Last column indicates sample shell thickness:
# 2 = long (thick)
# 1 = medium
# 0 = short (thin)
print("Oganizing data ...")
dataList = []
spectraLRhnms = []
for i in range(3):
    sampleNum = 2-i
    dlsData = data[2*i]
    spectraData0 = data[2*i+1]
    temps = spectraTemps[i]
    hydroRs = []
    for T in temps:
        aList = list(zip(dlsData[:,0], dlsData[:,1], np.abs(dlsData[:,0] - T)))
        dtype = [('T',float), ('Rh',float), ('dT',float)]
        anArray = np.array(aList, dtype=dtype)
        anArray = np.sort(anArray, order='dT')[0:3]
        subset = []
        for pt in anArray:
            subset.append(pt[1])
        rh = np.mean(subset)
        hydroRs.append(rh)
    for ii in range(len(temps)):
        for iii in range(len(spectraData0)):
            dataPoint = [ spectraData0[iii,0], spectraData0[iii,ii+1], temps[ii], hydroRs[ii].tolist(), sampleNum ]
            dataList.append(dataPoint)
            spectraLRhnms.append(hydroRs)
dataLArray = np.asarray(dataList)
spectraLTemps = spectraTemps
print("  * extrapolate hydrodynamic radius from DLS to UV-Vis temps")
print("  * organized as [ wl-nm, abs, temp, Rh, sample-# ]")
print("  * sample 2 = thick shell, 1 = medium, 0 = thin")

# Au cores only
spectrumCore = data[-1]
print("  * core alone sorted separately")

print("Done.\n")

# Cross-linked PNIPAM shells
print("Import and organize data")
print("========================")
print("X-linked PNIPAM: ")
folderName = 'data_20160906'
fileNames = []
print("Importing data in folder",folderName,"...")
for file in os.listdir(folderName):
    #if file.endswith('.dat'): 
        fileNames.append(file)
spectraXTemps = []
fileNames = np.sort(fileNames)

def importData(fileName):
    print("  *",fileName)
    headerLines = 3
    if (fileName == 'temperature-calibration.dat'): 
        headerLines = 2
    if (fileName == 'Au-xPNIPAM-long.dat'
          or fileName == 'Au-xPNIPAM-medium.dat'
          or fileName == 'Au-xPNIPAM-short.dat'): 
        temps = np.genfromtxt(folderName+'/'+fileName, delimiter='\t', skip_header=2, max_rows=1)[1:-1]
        temps = temps.tolist()
        spectraXTemps.append(temps)
        #print(temps)

    tempData = np.genfromtxt(folderName+'/'+fileName, delimiter='\t', skip_header=headerLines)
    # print(tempData)
    return tempData
data = list(map(importData,fileNames))

tempConvert = interpolate.interp1d(data[-1][:,0], data[-1][:,1], kind="cubic")
aListOfLists = spectraXTemps
spectraXTemps = []
for i in range(3):
    aList = []
    for ii in range(len(aListOfLists[i])):
        newTemp = (tempConvert(aListOfLists[i][ii])).tolist()
        aList.append(newTemp)
    spectraXTemps.append(aList)

# Last column indicates sample shell thickness:
# 5 = long (thick)
# 4 = medium
# 3 = short (thin)
print("Oganizing data ...")
dataXList = []
spectraXRhnms = []
for i in range(3):
    sampleNum = int(3 + 2-i)
    dlsData = data[1+2*i]
    spectraData0 = data[1+2*i+1]
    temps = spectraXTemps[i]
    hydroRs = []
    for T in temps:
        aList = list(zip(dlsData[:,0], dlsData[:,1], np.abs(dlsData[:,0] - T)))
        dtype = [('T',float), ('Rh',float), ('dT',float)]
        anArray = np.array(aList, dtype=dtype)
        anArray = np.sort(anArray, order='dT')[0:3]
        subset = []
        for pt in anArray:
            subset.append(pt[1])
        rh = np.mean(subset)
        hydroRs.append(rh)
    for ii in range(len(temps)):
        for iii in range(len(spectraData0)):
            dataPoint = [ spectraData0[iii,0], spectraData0[iii,ii+1], temps[ii], hydroRs[ii].tolist(), sampleNum ]
            dataXList.append(dataPoint)
            spectraXRhnms.append(hydroRs)
dataXArray = np.asarray(dataXList)
print("  * extrapolate hydrodynamic radius from DLS to UV-Vis temps")
print("  * organized as [ wl-nm, abs, temp, Rh, sample-# ]")
print("  * sample 5 = thick shell, 4 = medium, 3 = thin")

# Au cores only
spectrumXCore = data[0]
print("  * core alone sorted separately")

print("Done.\n")


# Combined data sets
print('Combining data sets...')
dataList = dataLArray.tolist()
for i in range(len(dataXArray)): dataList.append(dataXArray[i])
dataArray = np.asarray(dataList)

spectraCores = np.array([spectrumCore, spectrumXCore])
coreRnms = [17/2, 13.2/2]

spectraTemps = []
spectraRhnms = []
for sampleNum in range(6):
    selection = (dataArray[:,-1] == sampleNum)
    temps = np.unique( (dataArray[selection])[:,-3] )
    Rhnms = []
    for T in temps:
        selection = (dataArray[:,-1] == sampleNum) & (dataArray[:,-3] == T)
        RhOfT = np.unique( (dataArray[selection])[:,-2] ) 
        #print('  ', T, RhOfT[0])
        Rhnms.append( RhOfT[0] )
    Rhnms = np.asarray(Rhnms)
    print(sampleNum,':\t', temps.size, 'temperatures')
    spectraTemps.append(temps)
    spectraRhnms.append(Rhnms)
spectraTemps = np.asarray(spectraTemps)
spectraRhnms = np.asarray(spectraRhnms)


print("\nImport and organize data")
print("========================")
print("X-linked PNIPAM: ")
folderName = 'data_20160914'
fileNames = []
print("Importing data in folder",folderName,"...")
for file in os.listdir(folderName):
    #if file.endswith('.dat'): 
        fileNames.append(file)
slsNames = []
slsTemps = []
fileNames = np.sort(fileNames)
def importData(fileName):
    print("  *",fileName)
    headerLines = 3

    columnNames = np.genfromtxt(folderName+'/'+fileName, 
                                delimiter='\t', dtype=str, 
                                comments='\\', skip_header=0, max_rows=1)[0:-1]
    slsNames.append(columnNames.tolist())
    temp = np.genfromtxt(folderName+'/'+fileName, delimiter='\t', dtype=str, 
                                comments='\\', skip_header=2, max_rows=1)[2]
    temp = float(temp)
    slsTemps.append(temp)
    #print(columnNames)
    temporaryData = np.genfromtxt(folderName+'/'+fileName, delimiter='\t', skip_header=headerLines)
    return temporaryData
slsData = list(map(importData,fileNames))
print('Columns:',slsNames[2])
print('Collapsed-state reference samples 4-5')
slsIvsqRefSample4 = slsData[3][:,[1,2]]
slsIvsqRefSample5 = slsData[1][:,[1,2]]

# Print out variables for easy reference

print('\n\nVariables:')
print('  * indexes:')
print('    - \"coreNum\" (2): ')
print('      (0=Lin),\t(1=X-L)')
print('    - \"sampleNum\" (6): ')
print('      (0=L-Thn),\t(1=L-Med),\t(2=L-Thk),')
print('      (3=X-Thn),\t(4=X-Med),\t5=(X-Thk)')
print('  * \"spectraCores\" = core spectra')
print('    - type:\t',type(spectraCores))
print('    - size:\t',np.shape(np.asarray(spectraCores)))
print('    - [coreNum]  x  [wl-index]  x  [ wl-nm, ext-au ]')
print('  * \"coreRnms\" = core radiuses')
print('    - type:\t',type(coreRnms))
print('    - size:\t',np.shape(np.asarray(coreRnms)))
print('    - [coreNum] x [core-radius-TEM-nm]')
print('  * \"dataArray\" = core-shell spectra')
print('    - type:\t',type(dataArray))
print('    - size:\t',np.shape(dataArray))
print('    - [data-index]  x  [ wl-nm, ext-au, temp-C, Rh-nm, sampleNum ]')
print('  * \"spectraTemps\" = core-shell temperatures')
print('    - type:\t',type(spectraTemps))
print('    - size:\t',np.shape(spectraTemps))
print('    - [sampleNum]  x  [ temp-C ]')
print('  * \"spectraRhnms\" = core-shell hydrodynamic radiuses')
print('    - type:\t',type(spectraTemps))
print('    - size:\t',np.shape(spectraTemps))
print('    - [sampleNum]  x  [ Rh-nm ]')
print('  * \"slsData\" = all 5 SLS data sets, sorted as input above')
print('    - type:\t',type(slsData))
print('    - size:\t',np.shape(slsData))
print('    - [ theta-deg, q-nm, I(q), error-I(q) ]')
print('  * \"slsIvsqRefSample4\" = collapsed state for sample 4 (X2-med)')
print('    - type:\t',type(slsIvsqRefSample4))
print('    - size:\t',np.shape(slsIvsqRefSample4))
print('    - [ q-nm, I(q) ]')
print('  * \"slsIvsqRefSample5\" = collapsed state for sample 5 (X3-thk)')
print('    - type:\t',type(slsIvsqRefSample5))
print('    - size:\t',np.shape(slsIvsqRefSample5))
print('    - [ q-nm, I(q) ]')
print("Done.\n")

