import re
import numpy
import sympy

indefinite_entry = True


def separate_chemical_equation(equation: str) -> dict[str]:
    if equation is None or not equation.strip():
        print("[ERROR]: An invalid equation was entered.")
        return None

    reactants: dict[str] = {}
    product: dict[str] = {}
    
    try:
        raw_split_equation = (equation.upper()).split("->")
        reactants = [item.strip() for item in raw_split_equation[0].split("+")]
        products = [item.strip() for item in raw_split_equation[1].split("+")]
    except:
        print("[ERROR] An invalid equation was entered.")
        return None

    if reactants is None or products is None:
        print("[ERROR]: Could not separate the chemical reactants or products.")
        return None
    return reactants, products


def format_chemical_equation(rxn_segment: dict[str]) -> list:
    if rxn_segment is None:
        print("[ERROR]: An invalid equation half was entered.")
        return None
    formatted_equation: list = []
    
    for part in rxn_segment:
        chemical_match = re.match(r'^(\d+)?\s*(.*)', part.strip())
        if not chemical_match:
            continue
            
        coeffecient = int(chemical_match.group(1)) if chemical_match.group(1) else 1
        chemical_formula = chemical_match.group(2)
        
        atom_matches = re.findall(r'([A-Z][a-z]*)(\d*)', chemical_formula)
        organized_atoms: list = []
        
        for symbol, subscript in atom_matches:
            count = int(subscript) if subscript else 1
            total_count = count * coeffecient
            
            organized_atoms.append({"coeffecient": total_count, "atom": symbol.lower()})
            
        if organized_atoms is not None:
            formatted_equation.append(organized_atoms)


    if formatted_equation is None:
        print("[ERROR]: Could not format a chemical equation half.")
        return None  
    return formatted_equation


def construct_equation_matrix(rxn_reactants: list, rxn_products: list):
    if rxn_reactants is None or rxn_products is None:
        print("[ERROR]: The reaction products or reactants cannot be used.")
        return None
    
    combined_atoms = set()
    for part in rxn_reactants + rxn_products:
        for entry in part:
            combined_atoms.add(entry['atom'])

    sorted_atoms = sorted(list(combined_atoms))
    atom_index = {atom: i for i, atom in enumerate(sorted_atoms)}

    elements_present = len(sorted_atoms)
    element_amount = len(rxn_reactants) + len(rxn_products)

    equation_matrix = sympy.zeros(elements_present, element_amount)
    current_col = 0
    for molecule_list in rxn_reactants:
        for entry in molecule_list:
            row_idx = atom_index[entry['atom']]
            equation_matrix[row_idx, current_col] += entry['coeffecient']
        current_col += 1

    for molecule_list in rxn_products:
        for entry in molecule_list:
            row_idx = atom_index[entry['atom']]
            equation_matrix[row_idx, current_col] -= entry['coeffecient']
        current_col += 1


    if sorted_atoms is None or equation_matrix is None:
        print("[ERROR]: Could not format a chemical equation matrix.")
        return None  
    return equation_matrix


def solve_chemical_equation(chemical_matrix):
    if chemical_matrix is None:
        print("[ERROR]: An invalid chemical equation matrix was passed.")
        return None

    null_space = chemical_matrix.nullspace()
    balanced_coeffecients = 0

    if null_space:
        solution = null_space[0]
        denominator = sympy.lcm([val.q for val in solution])
        balanced_coeffecients = solution * denominator

        
    if balanced_coeffecients is None or balanced_coeffecients == 0:
        print("[WARNING]: No possible solutions were found.")
        return None
    return balanced_coeffecients

def format_solved_chemical_equation(final_coeffecients, raw_reactants: list, raw_products: list):
    all_molecules = raw_reactants + raw_products
    formatted_parts: list[str] = []
    for i in range(len(all_molecules)):
        display_coeffecient = str(final_coeffecients[i]) if int(final_coeffecients[i]) > 1 else ""
        formatted_parts.append(f"{display_coeffecient}{all_molecules[i]}")

    num_reactants = len(raw_reactants)
    reactants_side = " + ".join(formatted_parts[:num_reactants])
    products_side = " + ".join(formatted_parts[num_reactants:])

    final_equation = f"{reactants_side} -> {products_side}"
    
    if final_equation is None:
        print("[ERROR]: Could not format final chemical equation.")
        return None
    return final_equation


if __name__ == "__main__":
    while indefinite_entry:
        print("Equation format: A + B -> C + D")
        user_equation = str(input("Enter a chemical equation: "))

        separated_equations = separate_chemical_equation(user_equation)
        raw_reactants = separated_equations[0]
        raw_products = separated_equations[1]

        formatted_reactants = format_chemical_equation(raw_reactants)
        formatted_products = format_chemical_equation(raw_products)

        raw_chemical_matrix = construct_equation_matrix(formatted_reactants, formatted_products)
        balanced_chemical_matrix = solve_chemical_equation(raw_chemical_matrix)
        solved_user_equation = format_solved_chemical_equation(balanced_chemical_matrix, raw_reactants, raw_products)

        print("Solved Chemical Equation:  ", solved_user_equation, "\n")
