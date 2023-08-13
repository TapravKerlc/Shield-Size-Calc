# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 23:17:06 2022

@author: Zverina
"""
import random
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from sklearn import linear_model
import pandas as pd
# from mendeleev import element
import chemicals


#########fotoni##########

energije = np.linspace(1000, 1000000000, 41) #upoštevamo energije od 1 KeV do 1 GeV (čeprav je to čisti overkill, pomembne energije so do nekaj MeV)
kovine = [13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 55, 56, 57, 72, 73, 74, 75, 76, 77, 78, 79, 81, 82] #vrstna števila kovin
skupnipresek = []
print(len(kovine))

#tukaj lahko dodamo lanth kovine
lanth = (range(58,72))
for elem in lanth:
    kovine.append(elem)
#tukaj lahko dodamo actind kovine    
actind = ((range(90, 104)))
for elem in actind:
    kovine.append(elem)
kovine.sort()
#print(kovine)

maxim = 0
kater = 0
tocke = []
def fotoefekt (z, E):
    presek = z**4/(E**3) #enačba za presek pri fotoefektu
    return presek

def compton (z, E):
    presek = z #enačba za presek pri comptonu
    return presek

def tvorbapar (z, E):
    presek = z**2 * math.log(E) #enačba za presek pri tvorbi parov
    return presek

# print(kovine[89])

for j in range(69):
    kovinaizgube = []
    skupni = 0
    #print(kovine[j])
    for i in energije: #gremo čez naše željene energije
        
        gama = i/0.511
        beta=math.sqrt(1-(1/gama)**2) #relativistična enačba za beta, oziroma hitrost/energijo delca

        simbol = chemicals.periodic_table[kovine[j]].symbol
        A=chemicals.periodic_table[simbol].MW
        z = kovine[j]
        izgube = z/(A*beta**2) 
        kovinaizgube.append(izgube)
        skupni = skupni + izgube
        
    lm = linear_model.LinearRegression() #initialize linear model, da predvidevamo črto odziva LOR
    y = np.array(kovinaizgube).reshape(-1, 1)
    x = np.array(energije).reshape(-1, 1)
    lm.fit(x, y)
    lm.predict(x)
    tocke.append(lm.score(x, y)) #dodajamo rezultat na naš seznam za vsako kovino z
    #print(lm.score(x,y))
    if maxim < (lm.score(x, y)): #navaden if za preverjanje če je maximum presežen
        maxim = lm.score(x, y)
        #print(maxim)
        kater = kovine[j] #našli smo boljšo kovino
    

    
    # if skupni < minim:
    #     minim = skupni
    #     kater = j
    # skupniizgube.append(skupni)
print(kater, maxim)
    
#y = skupniizgube #dependent
#x = kovine #independent
# x = sm.add_constant(x) #dodamo konst
# lm = sm.OLS(y,x).fit() #fit
# #lm.predict(x)
# print(lm.summary())

# lm = linear_model.LinearRegression()
# y = np.array(skupniizgube).reshape(-1, 1)
# x = np.array(kovine).reshape(-1, 1)
# lm.fit(x, y)
# lm.predict(x)
# print(lm.score(x, y))



data = pd.DataFrame([kovine, tocke]).T
data.columns = ['kovine', 'tocke']

sns.lmplot(x="kovine", y="tocke", data=data, order=1) #za risanje zahtevnejšega plota 
plt.title("$R^2$ linearne regresije za kovine")
plt.ylabel('linearnost energijskega odziva')
plt.xlabel('Kovine')
plt.show

########težji delci, alfa, beta sevanje#######
#TODO
# def izgubeTezki(rho, Z, A, E): #gostota, masno, atomsko in energija delca
#     izgubeT = rho * Z / (A*E)
    
    
#for i in energije: