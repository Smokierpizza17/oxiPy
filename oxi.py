'''takes a chemical structure (in SMILES or just simple plaintext) and computes
oxidisation states of the individual atoms.'''
import re, sys

NONMETALELEMENTS = ["H", "He", "C", "N", "O", "F", "P", "S", "Cl", "Se", "Br",
                    "I", "At", "Ts"]

AUTOSETTABLEOXINUMBERS = {"F": -1, "O": -2,
                          "Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1, "Fr": 1,
                          "Be": 2, "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2, "Ra": 2,
                          "Sc": 3, "Y": 3, "Zr": 4, "Hf": 4, "Ta": 5, "Tc": 7,
                          "Zn": 2, "Cd": 2, "B": 3, "Al": 3, "Ga": 3, "In": 3,
                          "Ge": 4, "Kr": 2, "Rn": 2,
                          "La": 3, "Ac": 3, "Th": 4, "Nd": 3, "Pm": 3, "Gd": 3,
                          "Dy": 3, "Ho": 3, "Es": 3, "Er": 3, "Fm": 3, "Md": 3,
                          "Lu": 3, "Lr": 3}

MOLECULARIONS = {"CO3": ((['C', 1, 4], ['O', 3, -2]), -2),
                 "SO3": ((['S', 1, 4], ['O', 3, -2]), -2),
                 "SO4": ((['S', 1, 6], ['O', 4, -2]), -2),
                 "NO2": ((['N', 1, 3], ['O', 2, -2]), -1),
                 "NO3": ((['N', 1, 5], ['O', 3, -2]), -1),
                 "PO3": ((['P', 1, 3], ['O', 3, -2]), -3),
                 "PO4": ((['P', 1, 5], ['O', 4, -2]), -3),
                 "CN": ((['C', 1, 2], ['N', 1, -3]), -1),
                 "SCN": ((['S', 1, -2], ['C', 1, 4], ['N', 1, -3]), -1),
                 "OH": ((['O', 1, -2], ['H', 1, 1]), -1),
                 "NH4": ((['N', 1, -3], ['H', 4, 1]), 1)}

overallChargesR = re.compile(r"[0-9]?[+,-]")
elementGroupsR = re.compile(r"[A-Z][a-z]?[0-9]*|(?:\().*(?:\))[0-9]*")
elementNamesR = re.compile(r"[A-Z][a-z]?")
elementCountsR = re.compile(r"[0-9]+")

subgroupNamesR = re.compile(r"(?<=\().*(?=\))")
subgroupCountsR = re.compile(r"(?<=\))[0-9]")

capitalLettersR = re.compile(r"[A-Z]")
numbersR = re.compile(r"[0-9]")

subscript = str.maketrans("0123456789+-", "₀₁₂₃₄₅₆₇₈₉₊₋")
superscript = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")


def spaceinator(oxiNumberStringlet, elementStringlet, delimiter):
    '''prepends and appends spaces so that things line up'''
    numberIndex = re.search(numbersR, oxiNumberStringlet).start()
    letterIndex = re.search(capitalLettersR, elementStringlet).start()
    while numberIndex != letterIndex:
        numberIndex = re.search(numbersR, oxiNumberStringlet).start()
        letterIndex = re.search(capitalLettersR, elementStringlet).start()
        if numberIndex < letterIndex:
            oxiNumberStringlet = delimiter + oxiNumberStringlet
        if numberIndex > letterIndex:
            elementStringlet = delimiter + elementStringlet

    while len(oxiNumberStringlet) < len(elementStringlet):
        oxiNumberStringlet += delimiter
    while len(oxiNumberStringlet) > len(elementStringlet):
        elementStringlet += delimiter
    
    return oxiNumberStringlet, elementStringlet


def oxiNumbersMissing(elements):
    '''returns list of indices where oxiNumbers missing
    (num of elements[n][2] == None)'''
    indices = []
    for index, element in enumerate(elements):
        if element[2] is None:
            indices.append(index)
    return indices


def interpretElements(rawString, oxiFiller=None):
    '''get formatted elements without oxi Numbers (fill with oxiFiller)'''
    try:
        overallChargeStr = overallChargesR.findall(rawString)[-1]

        if overallChargeStr == "0":
            overallCharge = 0
        if overallChargeStr == "+":
            overallCharge = 1
        elif overallChargeStr == "-":
            overallCharge = -1
        else:
            overallCharge = int(overallChargeStr[::-1])  # reverse (2+ -> +2)
    except IndexError:
        overallCharge = 0

    # split up rawString into format [[name, count, oxiNumber], ...]
    # also identify subgroups and solve separately
    elementGroups = elementGroupsR.findall(rawString)
    elements = []

    for elementGroup in elementGroups:
        if "(" in elementGroup:
            elementName = subgroupNamesR.findall(elementGroup)[0]
            elementCount = subgroupCountsR.findall(elementGroup)
        else:
            elementName = elementNamesR.findall(elementGroup)[0]
            elementCount = elementCountsR.findall(elementGroup)

        if len(elementCount) == 0:
            elementCount = 1
        else:
            elementCount = int(elementCount[0])

        if elementName in MOLECULARIONS.keys():
            elements.append((elementName, elementCount, oxiFiller))
        else:
            elements.append([elementName, elementCount, oxiFiller])

    return elements, overallCharge


def getOxiNumbers(rawString):
    elements, overallCharge = interpretElements(rawString)

    # if only one atom given, solve immediately and return
    if len(elements) == 1:
        elements[0][2] = overallCharge  # oxi number automatically charge for
        # elemental compounds
        return elements, overallCharge

    # indentify subgroups and solve from lookup table
    for pos, element in enumerate(elements):
        if element[0] in MOLECULARIONS.keys():
            subgroupName = element[0]
            replacementElement = [MOLECULARIONS[subgroupName][0], element[1],
                                  MOLECULARIONS[subgroupName][1]]

            elements[pos] = replacementElement

    # assign +1 or -1 to Hydrogen based on neighbouring atoms
    for pos, element in enumerate(elements):
        if type(element[2]) is int:
            continue

        if element[0] == "H":
            # -1 if with metal, +1 if with nonmetal, priority -1
            try:
                if elements[pos + 1][0] not in NONMETALELEMENTS:
                    element[2] = -1
                elif elements[pos + 1][0] in NONMETALELEMENTS:
                    element[2] = 1
            except IndexError:
                pass
            try:
                if elements[pos - 1][0] in NONMETALELEMENTS:
                    element[2] = 1
                elif elements[pos - 1][0] not in NONMETALELEMENTS:
                    element[2] = -1
            except IndexError:
                pass
        if len(oxiNumbersMissing(elements)) <= 1:
            break

    # iterate and find auto-settable atoms, in order of importance
    for autoElement, autoOxiNumber in AUTOSETTABLEOXINUMBERS.items():
        for element in elements:
            if type(element[2]) == int:
                continue
            if element[0] == autoElement:
                element[2] = autoOxiNumber
            if len(oxiNumbersMissing(elements)) <= 1:
                break
        if len(oxiNumbersMissing(elements)) <= 1:
            break

    # compute missing oxiNumber from sum of others and overall charge
    missingIndices = oxiNumbersMissing(elements)
    if len(missingIndices) == 1:
        missingIndex = missingIndices[0]
        sumOfOxiNumbers = 0
        for element in elements:
            if type(element[2]) == int:
                sumOfOxiNumbers += element[1] * element[2]
        lastOxiNumber = 0 - sumOfOxiNumbers + overallCharge
        elements[missingIndex][2] = lastOxiNumber // elements[missingIndex][1]
        elements[missingIndex][2] = int(elements[missingIndex][2])
    return elements, overallCharge


def printResult(elements, overallCharge=0, delimiter=" ", passUp=False, verbose=False):
    '''prettyprints result with subscript indices and overhead oxiNumbers'''
    oxiNumbers = []
    elementStringlets = []
    if verbose:
        print(elements)

    # split up elements into parts of the full string (Stringlets™) that each
    # have one oxiNumber
    for element in elements:
        if type(element[0]) is tuple:
            newoxiNumbers, newElementStringlets = printResult(element[0],
                                                              0, " ", True)
            newElementStringlets[0] = "(" + newElementStringlets[0]  # add
            # starting bracket
            newElementStringlets[-1] += ")"  # add ending bracket
            oxiNumbers += newoxiNumbers
            elementStringlets += newElementStringlets
            continue

        # get oxiNumbers
        if element[2] is not None:
            oxiNumbers.append(str(element[2]))
        else:
            oxiNumbers.append("?")
        if element[1] >= 2:
            elementString = element[0] + str(element[1])
        else:
            elementString = element[0]
        elementStringlets.append(elementString)
    if passUp:
        return oxiNumbers, elementStringlets
    else:
        oxiNumberString = ""
        elementString = ""

        # extend either string so that length is equal
        for (oxiNumber, elementStringlet) in zip(oxiNumbers, elementStringlets):
            elementStringlet += delimiter
            paddedOxiNumber, paddedElement = spaceinator(oxiNumber, elementStringlet, delimiter)
            oxiNumberString += paddedOxiNumber
            elementString += paddedElement

        oxiNumberString = oxiNumberString.rstrip()
        elementString = elementString.rstrip()

        # create ChargeString
        chargeString = ""

        if overallCharge == -1:
            chargeString = "-"
        elif overallCharge == 1:
            chargeString = "+"

        elif overallCharge < 0:
            chargeString = str(abs(overallCharge)) + "-"
        elif overallCharge > 0:
            chargeString = str(abs(overallCharge)) + "+"

        chargeString = chargeString.translate(superscript)

        upperString = oxiNumberString
        lowerString = elementString.translate(subscript) + chargeString

        return upperString + "\n" + lowerString


def wrapper(rawString, oxiFiller=None, printDelimiter=" "):
    outputText = printResult(*getOxiNumbers(rawString), printDelimiter)
    return outputText


if __name__ == "__main__":
    query = ' '.join(sys.argv[1:])
    if query != "":
        print(printResult(*getOxiNumbers(query.strip())))
