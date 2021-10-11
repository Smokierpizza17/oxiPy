from oxi import *

questionsList = '''Cl2
Cl -
Na
Na +
O2
N2
Al 3+
H2O
NO3 -
NO2 -
NO2
Cr2O7 2-
KCl
NH3
CaH2
SO4 2-
Na2O2
SiO2
CaCl2
PO4 3-
MnO2
FeO
Fe2O3
H2O2
CaO
H2S
H2(SO4)
(NH4)Cl
K3(PO4)
H(NO3)
K(NO2)'''

if questionsList:
    for query in questionsList.split("\n"):
        if query != "":
            print(wrapper(query.strip()))
            print("")

counter = 1

while True:
    query = input("question %s > " % counter)
    if query != "":
        print("answer %s:" % counter)
        print(wrapper(query.strip()))
        print("")
    counter += 1