from root_pandas import read_root
import ROOT
from ROOT import TFile, TTree
import sys
import numpy as np
import pandas as pd


cols = ['PrimaryParticleEnergy','image_E1', 'image_E2', 'image_DRT_S', 'image_DRT_C', 'image_TT', 'theta_seed', 'phi_seed']

#for i in xrange(20*20): #change 20 to gridSize to change grid size 
    #cols.append('pixel_%i' % i)
    
#df = read_root(sys.argv[1], columns=cols, flatten=['image_E1', 'image_E2', 'image_DRT_S', 'image_DRT_C', 'image_TT'])
#df = read_root(sys.argv[1], columns=cols, flatten=['image_*'])
df = read_root(sys.argv[1], columns=cols)

inputFile = TFile(sys.argv[1], "READ")
myTree    = inputFile.Get("CNN")

#o_PrimaryParticleEnergy = []
#o_PrimaryParticleMomentum = []
o_PrimaryParticleName = []
#o_theta_seed = []
#o_phi_seed = []

#o_image_TT = []
#o_image_E1 = []
#o_image_E2 = []
#o_image_DRT_C = []
#o_image_DRT_S = []


counter = 0
for iEntry in myTree:         
    #if (counter <10) : 
    if (True) :         
        #o_PrimaryParticleEnergy.append(iEntry.PrimaryParticleEnergy)
        o_PrimaryParticleName.append(str(iEntry.PrimaryParticleName))
        #tt_var = []        
        #e1_var = []
        #e2_var = []
        #drtc_var = []
        #drts_var = []
        #if (counter%100 == 0): print ("event: ", counter)
        #for i in range(len(iEntry.image_TT)):
            #tt_var.append(float(iEntry.image_TT[i]))
            #e1_var.append(float(iEntry.image_E1[i]))
            #e2_var.append(float(iEntry.image_E2[i]))
            #drtc_var.append(float(iEntry.image_DRT_C[i]))
            #drts_var.append(float(iEntry.image_DRT_S[i]))

        #o_image_TT.append(tt_var)
        #o_image_E1.append(e1_var)
        #o_image_E2.append(e2_var)
        #o_image_DRT_C.append(drtc_var)
        #o_image_DRT_S.append(drts_var)        
        ##print (iEntry.PrimaryParticleName)
    #counter+=1


    
#print("read ", counter, " events")

#data = {'PrimaryParticleEnergy' : o_PrimaryParticleEnergy,
        #'PrimaryParticleName'   : o_PrimaryParticleName,
        #'image_TT'              : o_image_TT,
        #'image_E1'              : o_image_E1,
        #'image_E2'              : o_image_E2,
        #'image_DRT_C'           : o_image_DRT_C,
        #'image_DRT_S'           : o_image_DRT_S
       #}

       
#print("created data ")

#df = pd.DataFrame(data)#, columns = ['PrimaryParticleEnergy', 'PrimaryParticleName'])#, 'image_E1', 'image_E2'])
df['PrimaryParticleName'] = o_PrimaryParticleName
print("created dataframe")

print (df)

#df['GenDeltaR'] = df['GenDeltaR'].str[0]
#df.to_hdf("pixelTrain.h5",key='df',mode='w',encoding='utf-8')
#df.to_hdf(sys.argv[2],key='df',mode='w',encoding='utf-8')

df.to_hdf(sys.argv[2], key='df',mode='w',encoding='utf-8')

print("\nroot file: ", sys.argv[1], " --> converted to: ", sys.argv[2])
#print("\nroot file: ", inputFileName, " --> converted to: ", outputFile)        
print ("done")
