import re

element_mw_list = {
    'h': 1.00784,
    'he': 4.002602,
    'li': 6.941,
    'be': 9.0121831,
    'b': 10.811,
    'c': 12.0107,
    'n': 14.0067,
    'o': 15.9994,
    'f': 18.998403163,
    'ne': 20.1797,
    'na': 22.98976928,
    'mg': 24.3050,
    'al': 26.9815385,
    'si': 28.0855,
    'p': 30.973761998,
    's': 32.065,
    'cl': 35.453,
    'ar': 39.948,
    'k': 39.0983,
    'ca': 40.078,
    'sc': 44.955908,
    'ti': 47.867,
    'v': 50.9415,
    'cr': 51.9961,
    'mn': 54.938044,
    'fe': 55.845,
    'co': 58.933194,
    'ni': 58.6934,
    'cu': 63.546,
    'zn': 65.38,
    'ga': 69.723,
    'ge': 72.630,
    'as': 74.921595,
    'se': 78.971,
    'br': 79.904,
    'kr': 83.798,
    'rb': 85.4678,
    'sr': 87.62,
    'y': 88.90584,
    'zr': 91.224,
    'nb': 92.90637,
    'mo': 95.95,
    'tc': 98.0,
    'ru': 101.07,
    'rh': 102.90550,
    'pd': 106.42,
    'ag': 107.8682,
    'cd': 112.414,
    'in': 114.818,
    'sn': 118.710,
    'sb': 121.760,
    'te': 127.60,
    'i': 126.90447,
    'xe': 131.293,
    'cs': 132.90545196,
    'ba': 137.327,
    'la': 138.90547,
    'ce': 140.116,
    'pr': 140.90766,
    'nd': 144.242,
    'pm': 145.0,
    'sm': 150.36,
    'eu': 151.964,
    'gd': 157.25,
    'tb': 158.92535,
    'dy': 162.500,
    'ho': 164.93033,
    'er': 167.259,
    'tm': 168.93422,
    'yb': 173.054,
    'lu': 174.9668,
    'hf': 178.49,
    'ta': 180.94788,
    'w': 183.84,
    're': 186.207,
    'os': 190.23,
    'ir': 192.217,
    'pt': 195.084,
    'au': 196.966569,
    'hg': 200.592,
    'tl': 204.38,
    'pb': 207.2,
    'bi': 208.98040,
    'po': 209.0,
    'at': 210.0,
    'rn': 222.0,
    'fr': 223.0,
    'ra': 226.0,
    'ac': 227.0,
    'th': 232.0377,
    'pa': 231.03588,
    'u': 238.02891,
    'np': 237.0,
    'pu': 244.0,
    'am': 243.0,
    'cm': 247.0,
    'bk': 247.0,
    'cf': 251.0,
    'es': 252.0,
    'fm': 257.0,
    'md': 258.0,
    'no': 259.0,
    'lr': 262.0,
    'rf': 267.0,
    'db': 270.0,
    'sg': 271.0,
    'bh': 270.0,
    'hs': 277.0,
    'mt': 278.0,
    'ds': 281.0,
    'rg': 282.0,
    'cn': 285.0,
    'nh': 286.0,
    'fl': 289.0,
    'mc': 290.0,
    'lv': 293.0,
    'ts': 294.0,
    'og': 294.0
}
indefinite_entry = True



def format_formula(compound: str) -> list[str]:
    if compound is None or not compound.strip():
        print("[ERROR]: An invalid compound was entered.")
        return None
    
    sorted_elements = sorted(element_mw_list.keys(), key=len, reverse=False)
    element_pattern = "|".join(sorted_elements)
    
    raw_formatted_compound = re.findall(f"({element_pattern})|(\\d+)", compound.lower())
    raw_formatted_compound = [item for group in raw_formatted_compound for item in group if item]
    
    formatted_compound: list[str] = []
    for i in range(len(raw_formatted_compound)):
        entry = raw_formatted_compound[i]
        formatted_compound.append(entry)

        if entry.isalpha():
            if i + 1 == len(raw_formatted_compound) or raw_formatted_compound[i + 1].isalpha():
                formatted_compound.append('1')
    
    if formatted_compound is None:
        print("[ERROR]: Could not correctly format compound.")
        return None  
    return formatted_compound


def calculate_formula_mw(compound: list[str]) -> float:
    if compound is None:
        print("[ERROR]: Can't proceed due to incorrectly formatted compound.")
        return

    total_mw: float = 0.000
    for i in range(len(compound)):
        for elm, mw, in element_mw_list.items():
            if compound[i] == elm:
                total_mw = total_mw + (mw * int(compound[i+1]))

    if total_mw <= 0.0:
        print("[ERROR]: Failed to calculate molecular weight of the choosen compound.")
    return total_mw



if __name__ == "__main__":
    while indefinite_entry:
        user_compound = str(input("Enter a compound's molecular formula: "))

        compound_formula = format_formula(user_compound)
        compound_weight = calculate_formula_mw(compound_formula)

        if compound_weight is not None and compound_weight > 0.0:
            print("Molecular weight of ", user_compound.upper(), ":  ", round(compound_weight, 5), " g/mol\n")
