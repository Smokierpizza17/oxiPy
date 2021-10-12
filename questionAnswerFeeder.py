import sys
from oxi import *

generateSubgroups = False
interactive = False

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
K(NO2)
Ca(NO3)2'''

questionsList = '''OH -
O2 2-

OCN -
SCN -

PO4 3-
PO3 3-

SO4 2-
SO3 2-
S2O3 2-

NO3 -
NO2 -
NO -

CrO4 2-
Cr2O7 2-
AsO4 3-
MnO4 -

CO3 2-
C2O4 2-'''

if generateSubgroups:
    for query in questionsList.split("\n"):
        if query != "":
            rawElements = getOxiNumbers(query.strip())
            newElements = ((*rawElements[0],), rawElements[1])
            print("\"%s\": " % query.split(" ")[0], end="")
            print(newElements, end="")
            print(",")
        else:
            print("")
    # sys.exit()

if questionsList:
    for query in questionsList.split("\n"):
        if query != "":
            # print(wrapper(query.strip()))
            print(printResult(*getOxiNumbers(query)))
            print("")

if interactive:
    counter = 1
    while True:
        query = input("question %s > " % counter)
        if query != "":
            print("answer %s:" % counter)
            print(wrapper(query.strip()))
            print("")
        counter += 1