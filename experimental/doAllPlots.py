import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

utilspath=os.path.join("..","util")
sys.path = [utilspath]+sys.path
import myLandlabTools as tools

doPlots = True

# The dash list is a even size list that gives the ink on, ink off in pixels.
dashList = [(3,0),(1,1),(3,3),(7,7),(7,2),(5,2,2,2,5,2),(2,2,7,2,2,2),(3,3,2,2),(5,2,10,2),(10,5,20,5)] 

stormDict = {
    'high':[30., 17100.],
    'low':[15., 17100.],
    'intense': [50., 900.],
    'intense1': [50., 3600.],
    'intense2': [100., 900.],
    '1995-11-18':[15.6, 17100.],
    '1997-03-11':[34.8, 4500.],
    '2009-02-10':[18.0, 7200.],
    '2010-12-08':[16.8, 5400.],
    '2012-11-24':[16.8, 3900],
}

storm = 'intense1'

theFolder = 'valleyTopoIntense1'

df = tools.filesToDf(theFolder)
if df.shape[0]  < 1:
    print('dataframe has zero length, check theFolder={}'.format(theFolder))
    exit()

# print(df.head())
grad = []
for filename in df['filename']:
    Slope,ManningN,HCond,InfilDepth = tools.decodeLLFilename(filename)
    # get the slope value from the string
    gradient = Slope.strip('valley')
    grad.append(int(gradient))
df['grad'] = grad

# print(df)
slopes = df['Slope'].unique()
grad = np.sort(df['grad'].unique())
infd = np.sort(df['InfilDepth'].unique())
hcond = np.sort(df['HCond'].unique())
manning = np.sort(df['ManningN'].unique())
print('slopes=\n{}\n'.format(slopes))
print('gradients=\n{}\n'.format(grad))
print('infd=\n{}\n'.format(infd))
print('hcond=\n{}\n'.format(hcond))
print('manning=\n{}\n'.format(manning))

# check folders and create if not exist
for theOutPath in ['Hydro_ManningN','Flow_ManningN','Hydro_HydConduct','Flow_HydConduct',
    'Hydro_SoilInfDepth','Hydro_Slope','Flow_SlopeManningN']:
    p = os.path.join(theFolder,theOutPath)
    if not os.path.exists(p):
        os.makedirs(p)

if doPlots:
    # -------------- variation in manningN - Hydrograph ---------------------------------

    for HydC in hcond:  
        for InfilD in infd:
            for thisSlope in slopes:
                
                dfPlot = df[(df['Slope'] == thisSlope) & (df['HCond'] == HydC)  & (df['InfilDepth'] == InfilD)]

                plt.figure(figsize = (15,6))
                theLegend = []

                for m,filename in enumerate(dfPlot['filename']):
                    Slope,ManningN,HCond,InfilDepth = tools.decodeLLFilename(filename)
                    hydrograph_time, discharge_at_outlet, h, d = tools.readLLData(filename)

                    plt.plot(hydrograph_time, discharge_at_outlet,dashes=dashList[m])
                    theLegend.append(ManningN)

                plt.legend(theLegend)   
                plt.xlabel('Time (h)')
                plt.ylabel('Discharge (m$^3$/s)')
                # plt.title('Outlet Hydrograph for various values of Manning n \nStorm {} mm for {} s\ndem with {}, initial water infiltration depth = {} m, soil hydraulic conductivity = {} m/s'.format(
                #     stormDict[storm][0],stormDict[storm][1],thisSlope,InfilD,HydC))

                plt.savefig(theFolder+'/Hydro_ManningN/Hydro_varyingManning_{}_Inf{}_Hyd{}.png' .format(thisSlope, InfilD, HydC),dpi=300)

    # ----------------- variation in manningN - Max Discharge & Cum flow --------------------------


    for HydC in hcond:  
        for InfilD in infd:

            legend = []
            plt.figure(figsize = (15,10))

            m = 0 
            for gradient in grad:
                maxDischarge = []
                sumDischarge = []
                manning = []


                dfPlot = df[(df['grad'] == gradient) & (df['HCond'] == HydC)  & (df['InfilDepth'] == InfilD)]
                thisSlope = dfPlot['Slope'].unique()
                legend.append(thisSlope[0])

                for filename in dfPlot['filename']:
                    Slope,ManningN,HCond,InfilDepth = tools.decodeLLFilename(filename)
                    hydrograph_time, discharge_at_outlet, h, d = tools.readLLData(filename)

                    maxDischarge.append(max(np.abs(discharge_at_outlet)))
                    sumDischarge.append(sum(np.abs(discharge_at_outlet)[1:]*np.diff(np.array(hydrograph_time))*3600.))
                    manning.append(ManningN)

                plt.subplot(211)
                plt.plot(manning, maxDischarge,dashes=dashList[m])
                plt.subplot(212)
                plt.plot(manning, sumDischarge,dashes=dashList[m])
                m = m + 1

            plt.subplot(211)
            plt.legend(legend)
            plt.xlabel('Manning n')
            plt.ylabel('Maximum Discharge (m$^3$/s)')
            # plt.title('Discharge for various values of Manning n and terrain slope \nStorm {} mm for {} s\ninitial water infiltration depth = {} m, soil hydraulic conductivity = {} m/s'.format(
            #     stormDict[storm][0],stormDict[storm][1],InfilD,HydC))

            plt.subplot(212)
            plt.legend(legend)
            plt.xlabel('Manning n')
            plt.ylabel('Cumulative Discharge (m$^3$)')
            plt.savefig(theFolder+'/Flow_ManningN/MaxAndCumDisch_varyingManningSlope_Inf{}_Hyd{}.png' .format(InfilD, HydC),dpi=300)

    # variation in HydConduct - Hydrograph

    for ManningN in manning:  
        for InfilDepth in infd:
            for Slope in slopes:
 
                plt.figure(figsize = (15,6))
                theLegend = []

                for m,HCond in enumerate(hcond):

                    dfPlot = df[(df['ManningN'] == ManningN) & (df['Slope'] == Slope) & (df['HCond'] == HCond)  & (df['InfilDepth'] == InfilDepth)]
                    filename = dfPlot['filename'].unique()[0]

                    hydrograph_time, discharge_at_outlet, h, d = tools.readLLData(filename)

                    plt.plot(hydrograph_time, discharge_at_outlet,dashes=dashList[m])
                    theLegend.append('{} m/s'.format(HCond))

                plt.legend(theLegend)   
                plt.xlabel('Time (h)')
                plt.ylabel('Discharge (m$^3$/s)')
                # plt.title('Outlet Hydrograph for various values of soil hydraulic conductivity (m/s) \nStorm {} mm for {} s\ndem with valley{}, initial infiltration depth = {} m, Manning n = {}'.format(
                #     stormDict[storm][0],stormDict[storm][1],int(gradient),InfilD,Manning))

                plt.savefig(theFolder+'/Hydro_HydConduct/Hydro_varyingHydCond_{}_manning{}_inf{}.png' .format(Slope, ManningN, InfilDepth),dpi=300)
                
    # variation in HydConduct - Max Discharge & Cum flow

    for Manning in manning:  
        for InfilD in infd:

            legend = []
            plt.figure(figsize = (15,10))

            m = 0
            for gradient in grad:
                maxDischarge = []
                sumDischarge = []
                hydCon = []

                dfPlot = df[(df['grad'] == gradient) & (df['ManningN'] == Manning)  & (df['InfilDepth'] == InfilD)]
                thisSlope = dfPlot['Slope'].unique()
                legend.append(thisSlope[0])
               
                for filename in dfPlot['filename']:
                    Slope,ManningN,HCond,InfilDepth = tools.decodeLLFilename(filename)
                    hydrograph_time, discharge_at_outlet, h, d = tools.readLLData(filename)

                    maxDischarge.append(max(np.abs(discharge_at_outlet)))
                    sumDischarge.append(sum(np.abs(discharge_at_outlet)[1:]*np.diff(np.array(hydrograph_time))*3600.))
                    hydCon.append(HCond)

                dfNow = pd.DataFrame({'hydCon':hydCon, 'sumDischarge':sumDischarge, 'maxDischarge':maxDischarge })
                dfNow = dfNow.sort_values(by=['hydCon'])

                plt.subplot(211)
                plt.plot(dfNow['hydCon'], dfNow['maxDischarge'],dashes=dashList[m])
                plt.subplot(212)
                plt.plot(dfNow['hydCon'], dfNow['sumDischarge'],dashes=dashList[m])
                m = m + 1


            plt.subplot(211)
            plt.legend(legend)
            plt.xlabel('Hydraulic Conductivity (m/s)')
            plt.ylabel('Maximum Discharge (m$^3$/s)')
            # plt.title('Discharge for various values of soil hydraulic conductivity and terrain slope\nStorm {} mm for {} s\ninitial water infiltration depth = {} m, Manning n = {}'.format(
            #     stormDict[storm][0],stormDict[storm][1],InfilD,Manning))

            plt.subplot(212)
            plt.legend(legend)
            plt.xlabel('Hydraulic Conductivity (m/s)')
            plt.ylabel('Cumulative Discharge (m$^3$)')

            plt.savefig(theFolder+'/Flow_HydConduct/MaxAndCumDisch_varyingHydCondSlope_manning{}_inf{}.png' .format(Manning, InfilD),dpi=300)

    # variation in Slope - Hydrograph

    if 'valley' in theFolder:
        
        for Manning in manning:  
            for HydC in hcond:
                for InfilD in infd:
                    
                    dfPlot = df[(df['InfilDepth'] == InfilD) & (df['ManningN'] == Manning)  & (df['HCond'] == HydC)]

                    plt.figure(figsize = (15,6))
                    theLegend = []

                    for m,gradient in enumerate(grad):

                        dfPlot = df[(df['ManningN'] == Manning) & (df['grad'] == gradient) & (df['HCond'] == HydC)  & (df['InfilDepth'] == InfilD)]
                        filename = dfPlot['filename'].unique()[0]
                        
                        Slope,ManningN,HCond,InfilDepth = tools.decodeLLFilename(filename)
                        hydrograph_time, discharge_at_outlet, h, d = tools.readLLData(filename)

                        plt.plot(hydrograph_time, discharge_at_outlet,dashes=dashList[m])
                        theLegend.append(Slope)

                    plt.legend(theLegend)   
                    plt.xlabel('Time (h)')
                    plt.ylabel('Discharge (m$^3$/s)')
                    # plt.title('Outlet Hydrograph for various terrain slopes \nStorm {} mm for {} s\ninitial water infiltration depth = {} m, soil hydraulic conductivity = {} m/s, Manning n = {}'.format(
                    #     stormDict[storm][0],stormDict[storm][1],InfilD,HydC,Manning))

                    plt.savefig(theFolder+'/Hydro_Slope/Hydro_varyingSlope_inf{}_manning{}_hyd{}.png' .format(InfilD, Manning, HydC),dpi=300)

        # variation in Slope - Max Discharge & Cum flow

        for HydC in hcond:
            for InfilD in infd:

                dfPlot = df[(df['HCond'] == HydC)  & (df['InfilDepth'] == InfilD)]
                slope = []

                # for each slope get the gradient value and add to data frame
                for i, filename in enumerate(dfPlot['filename']):

                    Slope,ManningN,HCond,InfilDepth = tools.decodeLLFilename(filename)

                    # get the slope value from the string
                    #gradient = Slope.strip('nslope')
                    gradient = Slope.strip('valley')

                    grad = int(gradient)
                    slope.append(grad)

                dfPlot['gradient'] = slope

                # sort dataframe on gradient
                dfPlot.sort_values('gradient',inplace=True)

                # prepare for data extraction and plotting
                plt.figure(figsize = (15,10))
                legend = []

                m = 0

                for Manning in manning:

                    dfX = dfPlot[dfPlot['ManningN'] == Manning]
                    maxDischarge = []
                    sumDischarge = []
                    slope = []
                    legend.append(Manning)

                    for filename, grad in zip(dfX['filename'],dfX['gradient']):

                        Slope,ManningN,HCond,InfilDepth = tools.decodeLLFilename(filename)
                        hydrograph_time, discharge_at_outlet, h, d = tools.readLLData(filename)

                        # exclude flat topo
                        if grad > 0:

                            maxDischarge.append(max(np.abs(discharge_at_outlet)))
                            sumDischarge.append(sum(np.abs(discharge_at_outlet)[1:]*np.diff(np.array(hydrograph_time))*3600.))
                            slope.append(grad)


                    plt.subplot(211)
                    plt.plot(slope, maxDischarge,dashes=dashList[m])
                    plt.subplot(212)
                    plt.plot(slope, sumDischarge,dashes=dashList[m])
                    m = m + 1

                plt.subplot(211)
                plt.xlabel('Terrain Slope (degrees)')
                plt.ylabel('Maximum Discharge (m$^3$/s)')
                # plt.title('Discharge for various values of terrain slope and Manning n\nStorm {} mm for {} s\ninitial water infiltration depth = {} m, soil hydraulic conductivity = {} m/s'.format(
                #     stormDict[storm][0],stormDict[storm][1],InfilD,HydC))
                plt.legend(legend)

                plt.subplot(212)
                plt.xlabel('Terrain Slope (degrees)')
                plt.ylabel('Cumulative Discharge (m$^3$)')
                plt.legend(legend)

                plt.savefig(theFolder+'/Flow_SlopeManningN/MaxAndCumDisch_varyingSlopeManningN_HydCond_{}_inf{}.png' .format(HydC, InfilD),dpi=300)
                    
