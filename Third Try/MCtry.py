import numpy as np
import matplotlib.pyplot as plt
import xcom
import math
import csv
import periodictable
#import os

data1 = 0
data1 = xcom.calculate_cross_section(1)
#print(data1)
material = xcom.Material([1]) # Create hydrogen
data2 = xcom.calculate_attenuation(material)
#print(data2)

#thanks https://plotdigitizer.com/app


def correct_data(energies, data, correction_csv):
    """
    Corrects the data using correction factors from a CSV file.

    Args:
    - energies (list): A list of energies used to measure the data.
    - data (list): A list containing the uncorrected data.
    - correction_csv (str): Path to the CSV file containing correction factors.
    
    Returns:
    - list: A list containing the corrected data.
    """
    # Load correction factors from CSV file
    correction_factors = {}
    with open(correction_csv, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header row if present
        for row in reader:
            energy, correction_factor = map(float, row)
            correction_factors[energy] = pow(correction_factor, 0.3)
    
    # Divide energies by 1000 to convert from keV to eV
    energies = [energy / 1000 for energy in energies]

    # Interpolate correction factors for each energy in ListByE
    corrected_data = []
    for energy, measurement in zip(energies, data):
        # Check if energy is within the range of available correction factors
        if energy < min(correction_factors.keys()):
            # Extrapolate below the minimum energy using the first two available data points
            energy1 = min(correction_factors.keys())
            energy2 = sorted(correction_factors.keys())[1]
            factor1 = correction_factors[energy1]
            factor2 = correction_factors[energy2]
            corrected_measurement = np.interp(energy, [energy1, energy2], [factor1, factor2]) * measurement
        elif energy > max(correction_factors.keys()):
            # Extrapolate above the maximum energy using the last two available data points
            energy2 = max(correction_factors.keys())
            energy1 = sorted(correction_factors.keys())[-2]
            factor2 = correction_factors[energy2]
            factor1 = correction_factors[energy1]
            corrected_measurement = np.interp(energy, [energy1, energy2], [factor1, factor2]) * measurement
        else:
            # Interpolate correction factor using nearest energies
            nearest_energies = sorted(correction_factors.keys(), key=lambda x: abs(x - energy))[:2]
            energy1, energy2 = nearest_energies
            factor1 = correction_factors[energy1]
            factor2 = correction_factors[energy2]
            corrected_measurement = np.interp(energy, [energy1, energy2], [factor1, factor2]) * measurement
        
        corrected_data.append(corrected_measurement)
    
    return corrected_data


def compton(entryE, začetnKot): #funkcija ki izračuna nov kot in energijo po comptonski interakciji
    alfa = entryE/1.022 #1.022= 2mc^2
    epsilon0 = 1/(1+2*alfa)
    alfa1 = (1-epsilon0**2)/2
    alfa2 = -math.log(epsilon0)    
    bully = True
    
    ####Sledi ful zakomplicirana funkcija za določanje epsilona, skupej s preverjanjem, da je prava rešitev (bully = True pomeni napačna rešitev)
    
    while bully:
        xi1 = np.random.random(1)
        xi2 = np.random.random(1)
        if xi1 < alfa1/(alfa1+alfa2):
            epsilonČ = math.sqrt(epsilon0**2+2*(alfa1+alfa2)*xi1)
        else:
            epsilonČ = epsilon0*math.exp((alfa1+alfa2)*xi1-alfa1)
            
        if xi2 < (1-(epsilonČ*(1-epsilonČ**2))/((1+epsilonČ**2)*(alfa*epsilonČ)**2)): #preverjamo če rešitev ustreza
            bully = False
    
    novKot = math.acos(1-((1-epsilonČ)/(alfa*epsilonČ))) + začetnKot
    exitE = epsilonČ*entryE
    global odloženaE
    odloženaE = odloženaE+(1-epsilonČ)*entryE
    #print(exitE, novKot)
    return exitE, novKot


# def min_coefficient_of_variation(data):
#     """
#     Find the sublist with the lowest coefficient of variation in a list of lists of lists.

#     Args:
#     - data (list): A list containing lists of lists of numeric data.

#     Returns:
#     - tuple: A tuple containing the position (index) of the sublist with the lowest coefficient of variation 
#              and the coefficient of variation itself.
#     """

#     min_cv = float('inf')
#     min_cv_index = None

#     for idx1, sublist1 in enumerate(data):
#         for idx2, sublist2 in enumerate(sublist1):
#             mean = np.mean(sublist2)
#             std_dev = np.std(sublist2)
#             if mean != 0:
#                 cv = std_dev / mean
#                 if cv < min_cv:
#                     min_cv = cv
#                     min_cv_index = (idx1, idx2)
    
#     return min_cv_index, min_cv

def min_linear_sublist(data):
    """
    Find the sublist with the highest Pearson correlation coefficient in a list of materials, sublists of thickness,sublist of passed by energies.

    - V data (ki ga pridobimo iz filteredListByMat) je struktura organizirana hierarhično po:

    materialu (prvi nivo: idx1) (al, ti, fe, pb),
    debelini materiala (drugi nivo: idx2) (1 do 10 cm),
    energijah za vsak poskus (notranji nivo: sublist2).

    Ko funkcija min_linear_sublist vrne rezultat, imamo:

        best[0]: indeks materiala (idx1) v filteredListByMat,
        best[1]: indeks debeline (idx2) za ta material.

    Vrednost rekord pa predstavlja najvišji dosežen korelacijski koeficient, ki smo ga našli.
    """

    max_corr = -1  # Initialize maximum correlation to a low value
    max_corr_index = None

    for idx1, sublist1 in enumerate(data): #by material
        for idx2, sublist2 in enumerate(sublist1): #by thickness
            # Calculate Pearson correlation coefficient
            corr_coef = np.corrcoef(sublist2[:-1], sublist2[1:])[0, 1]  # Assuming data is 1D
            if abs(corr_coef) > max_corr:
                max_corr = abs(corr_coef)
                max_corr_index = (idx1, idx2)
    
    return max_corr_index, max_corr
#to find the sublist that is most linear instead of having the lowest coefficient of variation, you can use a measure of linearity such as the Pearson correlation coefficient.

odloženaE = 0

def plot(data):
    x = data["energy"]/1e6
    # plt.plot(x, data["coherent"], label="Coherent")
    # plt.plot(x, data["incoherent"], label="Incoherent")
    # plt.plot(x, data["photoelectric"], label="Photoelectric")
    # plt.plot(x, data["pair_atom"], label="Pair atom")
    # plt.plot(x, data["pair_electron"], label="Pair electron")
    plt.plot(x, data["total"], label="Total")
    plt.xlabel("MeV")
    plt.ylabel("barn/atom")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.show()
#print(data1["energy"]/1e6)
#plot(data1)

data = xcom.calculate_cross_section(82)
#plot(data)
#print(data["Material"])



barn_to_cm2 = 1e-24  # 1 barn = 1e-24 cm²
energies = data["energy"]#/1e6 #energije v eV
lower_bound = 50000  # 50 keV in eV
upper_bound = 3000000  # 3 MeV in eV
debel = np.linspace(0.1, 10, 99) #v cm
materiali = [13, 22, 26, 82] #aluminij, titan, železo, svinec
imenaMat = ["Aluminij", "Titan", "Železo", "Svinec"]

cena_na_kg = {
    "Aluminij": 2.614,
    "Titan": 6.600,
    "Železo": 0.461,
    "Svinec": 2.048
}

energies_in_range = energies[(energies >= lower_bound) & (energies <= upper_bound)]
energies = energies_in_range

counterPobegli = 0 #koliko jih pride skozi
print(energies)

dolzinaD = len(debel)
dolzinaE = len(energies)
dolzinaM = len(materiali)
ListByMat = [None]*dolzinaM
ListByD = [None]*dolzinaD
newListByD = [None]*dolzinaD
ListByMasa = [None]*dolzinaD


tries = 100
detector_radius = 9.53 / 10  # Pretvori v cm
detector_length = 45.2 / 10  # Pretvori v cm
max_weight = 2000  # 2 kg v gramih

for stev, mat in enumerate(materiali):  # index and atomic number
    data = xcom.calculate_cross_section(mat)  # calculate cross-section for material
    element = periodictable.elements[mat]  # get element properties
    print(imenaMat[stev])
    
    for stetje, debelina in enumerate(debel):  # index and thickness in cm
        # Izračun prostornine ovoja
        outer_radius = detector_radius + debelina  # Zunanji radij ovoja
        shield_volume = math.pi * (outer_radius**2 - detector_radius**2) * detector_length  # Prostornina ovoja v cm³

        # Izračun mase ovoja
        shield_mass = shield_volume * element.density  # Masa ovoja v g

        # Preveri, ali masa presega 1 kg
        if shield_mass > max_weight:
            #print(f"Debelina {debelina} cm za material {imenaMat[stev]} je preskočena zaradi teže {shield_mass/1000:.2f} kg.")
            continue  # Preskoči to debelino, če masa presega 1 kg

        print(stetje / dolzinaD * 100, "%", debelina, " cm")
        ListByE = [None] * dolzinaE
        skip_debelina = False  # Flag to skip this thickness

        for zaporedni, E1 in enumerate(energies):
            counterPobegli = 0  # Reset the escape counter for each energy
            exit_loops = False  # Flag to break out of nested loops early

            for x in range(tries):
                currentE = E1
                currentKot = 0
                dRoba = debelina

                while currentE > 0:
                    dRoba = math.cos(currentKot) * dRoba

                    # Incoherent cross-section
                    currentCSpathBarn = data["incoherent"][zaporedni]
                    currentCSpathCm2 = currentCSpathBarn * barn_to_cm2  # Convert barns to cm²
                    currentCSpath = 1 / (currentCSpathCm2 * (6.022e23 * element.density / element.mass))

                    # Photoelectric cross-section
                    currentFEpathBarn = data["photoelectric"][zaporedni]
                    currentFEpathCm2 = currentFEpathBarn * barn_to_cm2  # Convert barns to cm²
                    currentFEpath = 1 / (currentFEpathCm2 * (6.022e23 * element.density / element.mass))

                    dPE = -math.log(np.random.random(1)) / currentFEpath
                    dCS = -math.log(np.random.random(1)) / currentCSpath

                    if (dRoba <= dCS) and (dRoba <= dPE):  # escape
                        counterPobegli += 1
                        #print(f"Pobegli count za to energijo: {counterPobegli}")
                        currentE = 0

                        # Check if counterPobegli exceeds 95% of tries
                        if counterPobegli > 0.95 * tries:
                            exit_loops = True  # Break the loops early
                            #print(f"Napaka! counterPobegli ({counterPobegli}) je več kot dovoljenih 95%.")
                            break  # Exit inner while loop

                    elif (dCS <= dPE):  # Compton effect
                        currentE, currentKot = compton(currentE, currentKot)
                    else:  # Photoelectric effect
                        currentE = 0

                if exit_loops:
                    break  # Exit the outer 'for x in range(tries)' loop

            if counterPobegli < 0.05 * tries:
                skip_debelina = True  # Mark this thickness to be skipped
                break  # Break out of the energy loop

            # If not skipping, store the result for this energy
            ListByE[zaporedni] = counterPobegli

        # If debelina is flagged to be skipped, continue to the next thickness
        if skip_debelina:
            continue

        # If the debelina is valid (not skipped), apply corrections
        Corrected = correct_data(energies, ListByE, "zp1201EkeV.csv")
        ListByD[stetje] = Corrected
        ListByMasa[stetje] = shield_mass

    # Ensure that only valid (non-None) entries are passed to ListByMat
    ListByMat[stev] = [d for d in ListByD if d is not None]

# Before calling min_linear_sublist, ensure to filter None from ListByMat
filteredListByMat = [sublist for sublist in ListByMat if sublist is not None]

# Now call min_linear_sublist
best, rekord = min_linear_sublist(filteredListByMat)



photon_energies = {
    "Tehnecij 99m": 140e3,   # 140 keV
    "GM detektor tipični zgornji rob": 500e3,  # 500 keV
    "Cs137": 662e3    
}

colors = {
    "Cs137": 'red',
    "Tehnecij 99m": 'blue',
    "GM detektor tipični zgornji rob": 'green'
}



if best is not None:
    optimalna_debelina_cm = (best[1] / 10)+1  # Konvertirajte debelino v cm
    masa_scita = ListByMasa[best[1]]  # Masa ščita v g
    ime_materiala = imenaMat[best[0]]
    cena_materiala = cena_na_kg.get(ime_materiala, "N/A")  # Cena na kg
    cena_materiala = cena_materiala * masa_scita/1000
    cena_materiala = round(cena_materiala, 3)
    cena_materiala = str(cena_materiala).replace(".", ",")

    print("Optimalna debelina v cm:", optimalna_debelina_cm)
    print("Masa tega ščita v g:", masa_scita)
    print(rekord, ime_materiala, ListByD[best[1]])
    plt.plot(
    energies,
    [value / 100 for value in ListByD[best[1]]],
    label=f"Ščit iz materiala: {ime_materiala}, debeline: {optimalna_debelina_cm} cm\nCena: {cena_materiala} $")
    for label, energy in photon_energies.items():
        plt.axvline(x=energy, color=colors[label], linestyle='--', label=f"{label}: {energy / 1e3} keV")
    plt.xlabel("eV")
    plt.ylabel("Relativni odziv")
    plt.xscale("log")
    plt.legend()
    plt.show()
else:
    print("Nismo našli primerne.")
""" 
# Seznami za shranjevanje mase in cene za vsak material
mase_za_vsak_material = {material: [] for material in imenaMat}
cene_za_vsak_material = {material: [] for material in imenaMat}

# Zanka skozi vse materiale in debeline, da napolnimo podatke za maso in ceno
for stev, mat in enumerate(materiali):
    element = periodictable.elements[mat]
    ime_materiala = imenaMat[stev]
    for debelina in debel:
        # Izračun prostornine ovoja
        outer_radius = detector_radius + debelina
        shield_volume = math.pi * (outer_radius**2 - detector_radius**2) * detector_length
        
        # Izračun mase ovoja
        shield_mass = shield_volume * element.density  # Masa ovoja v g
        if shield_mass > max_weight:
            continue  # Preskoči debelino, če masa presega 1 kg
        
        # Shranjevanje mase in cene
        mase_za_vsak_material[ime_materiala].append(shield_mass / 1000)  # v kg
        cena_materiala = cena_na_kg[ime_materiala] * (shield_mass / 1000)
        cene_za_vsak_material[ime_materiala].append(cena_materiala)

# Risanje grafa za maso
plt.figure(figsize=(12, 6))
for ime_materiala in imenaMat:
    plt.plot(debel[:len(mase_za_vsak_material[ime_materiala])], mase_za_vsak_material[ime_materiala],
             label=f"{ime_materiala}")
plt.xlabel("Debelina ščita (cm)")
plt.ylabel("Masa ščita (kg)")
plt.title("Naraščanje mase z debelino ščita za različne materiale")
plt.legend()
plt.grid()
plt.show()

# Risanje grafa za ceno
plt.figure(figsize=(12, 6))
for ime_materiala in imenaMat:
    plt.plot(debel[:len(cene_za_vsak_material[ime_materiala])], cene_za_vsak_material[ime_materiala],
             label=f"{ime_materiala}")
plt.xlabel("Debelina ščita (cm)")
plt.ylabel("Cena ščita ($)")
plt.title("Naraščanje cene z debelino ščita za različne materiale")
plt.legend()
plt.grid()
plt.show()
 """
