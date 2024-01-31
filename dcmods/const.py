import math


def rp(field_strength, agent='Primovist'): # Gadoxetate
# Szomolanyi P, Rohrer M, Frenzel T, Noebauer-Huhmann IM, Jost G, Endrikat J, Trattnig S, Pietsch H. Comparison of the Relaxivities of Macrocyclic Gadolinium-Based Contrast Agents in Human Plasma at 1.5, 3, and 7 T, and Blood at 3 T. Invest Radiol. 2019 Sep;54(9):559-564. doi: 10.1097/RLI.0000000000000577.
    field = math.floor(field_strength)
    if agent=='Primovist':
        if field == 1.5: return 8.1     # relaxivity of blood in Hz/mM
        if field == 3.0: return 6.4     # relaxivity of blood in Hz/mM
        if field == 4.0: return 6.4     # relaxivity of blood in Hz/mM
        if field == 7.0: return 6.2     # relaxivity of blood in Hz/mM
        if field == 9.0: return 6.1     # relaxivity of blood in Hz/mM 
    if agent=='Clariscan':
        if field == 3.0: return 2.72 
    if agent=='Dotarem':
        if field == 3.0: return 2.72 
    msg = 'No relaxivity data for ' + agent + ' at ' + str(field_strength) + ' T.'
    raise ValueError(msg)

def rh(field_strength):
    field = math.floor(field_strength)
    if field == 1.5: return 14.6    # relaxivity of hepatocytes in Hz/mM
    if field == 3.0: return 9.8     # relaxivity of hepatocytes in Hz/mM
    if field == 4.0: return 7.6     # relaxivity of hepatocytes in Hz/mM
    if field == 7.0: return 6.0     # relaxivity of hepatocytes in Hz/mM
    if field == 9.0: return 6.1     # relaxivity of hepatocytes in Hz/mM

def R1_blood(field_strength=3.0, Hct=0.45):
    field = math.floor(field_strength)
    if field == 1.5: return 1000.0 / 1480.0    # aorta R1 in 1/sec 
    if field == 3.0: return 0.52 * Hct + 0.38  # Lu MRM 2004 

def R1_liver(field_strength=3.0):
    field = math.floor(field_strength)
    if field == 1.5: return 1000.0/602.0     # liver R1 in 1/sec (Waterton 2021)
    if field == 3.0: return 1000.0/752.0     # liver R1 in 1/sec (Waterton 2021)
    if field == 4.0: return 1.281     # liver R1 in 1/sec (Changed from 1.285 on 06/08/2020)
    if field == 7.0: return 1.109     # liver R1 in 1/sec (Changed from 0.8350 on 06/08/2020)
    if field == 9.0: return 0.920     # per sec - liver R1 (https://doi.org/10.1007/s10334-021-00928-x)

def R1_kidney(field_strength=3.0):
    # Reference values average over cortext and medulla from Cox et al
    # https://academic.oup.com/ndt/article/33/suppl_2/ii41/5078406
    field = math.floor(field_strength)
    if field == 1.5: return 1000.0/((1024+1272)/2)
    if field == 3.0: return 1000.0/((1399+1685)/2)