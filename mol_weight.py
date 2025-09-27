from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem

def calculate_molecular_properties(smiles):
    """
    Calculates various physicochemical properties for a given SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        dict: A dictionary containing the calculated properties, or None if the SMILES is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print(f"Error: Invalid SMILES string: {smiles}")
        return None

    properties = {}

    # Basic Descriptors
    properties['SMILES'] = smiles
    properties['Molecular Formula'] = Chem.rdMolDescriptors.CalcMolFormula(mol)
    properties['Molecular Weight (Monoisotopic)'] = Descriptors.ExactMolWt(mol)
    properties['Molecular Weight (Average)'] = Descriptors.MolWt(mol)
    properties['LogP (Octanol-Water Partition Coefficient)'] = Descriptors.MolLogP(mol)
    properties['TPSA (Topological Polar Surface Area)'] = Descriptors.descriptors.CalcTPSA(mol)
    properties['Num H Donors'] = Descriptors.NumHDonors(mol)
    properties['Num H Acceptors'] = Descriptors.NumHAcceptors(mol)
    properties['Num Rotatable Bonds'] = Descriptors.NumRotatableBonds(mol)
    properties['Num Rings'] = Descriptors.RingCount(mol)
    properties['Num Aromatic Rings'] = Descriptors.NumAromaticRings(mol)
    properties['Num Saturated Rings'] = Descriptors.NumSaturatedRings(mol)
    properties['Num Heteroatoms'] = Descriptors.NumHeteroatoms(mol)

    # Vina related properties (simplified examples)
    # These are often used in drug discovery for lead optimization
    properties['QED (Quantitative Estimation of Drug-likeness)'] = Descriptors.rdMolDescriptors.CalcQED(mol)
    properties['SAscore (Synthetic Accessibility Score)'] = Descriptors.rdMolDescriptors.CalcSAscore(mol)

    # More advanced descriptors (examples, RDKit has many more)
    properties['Fraction CSP3'] = Descriptors.FractionCSP3(mol)
    properties['Num Amide Bonds'] = Descriptors.NumAmideBonds(mol)

    # Generate 2D coordinates for visualization (optional)
    AllChem.Compute2DCoords(mol)
    properties['2D Coordinates Generated'] = True # Just to indicate it was attempted

    return properties

def display_properties(properties):
    """
    Prints the calculated molecular properties in a formatted way.

    Args:
        properties (dict): A dictionary of molecular properties.
    """
    if properties:
        print("\n--- Molecular Properties ---")
        for key, value in properties.items():
            print(f"{key}: {value}")
        print("----------------------------")
    else:
        print("No properties to display (SMILES might have been invalid).")

if name == "main":
    # Example usage:
    smiles_molecules = [
        "CCO",        # Ethanol
        "C1=CC=C(C=C1)C(=O)O", # Benzoic Acid
        "CC(=O)Oc1ccccc1C(=O)O", # Aspirin
        "C(C(=O)O)N", # Glycine
        "INVALID_SMILES" # An intentionally invalid SMILES
    ]

    for smiles in smiles_molecules:
        print(f"\nProcessing SMILES: {smiles}")
        mol_properties = calculate_molecular_properties(smiles)
        display_properties(mol_properties)

        # You can now access individual properties
        if mol_properties:
            print(f"  Molecular Weight of {smiles}: {mol_properties['Molecular Weight (Monoisotopic)']:.2f}")

    print("\n--- Visualization Example (requires RDKit drawing capabilities and Matplotlib) ---")
    try:
        from rdkit.Chem.Draw import MolToImage
        from PIL import Image
        import matplotlib.pyplot as plt

        smiles_to_draw = "CC(=O)Oc1ccccc1C(=O)O" # Aspirin
        mol_draw = Chem.MolFromSmiles(smiles_to_draw)
        if mol_draw:
            img = MolToImage(mol_draw)
            plt.imshow(img)
            plt.title(f"2D
 Structure of {smiles_to_draw}")
            plt.axis('off')
            plt.show()
            # If you want to save it
            # img.save("aspirin.png")
            print(f"Generated 2D structure for {smiles_to_draw} and displayed/saved.")
        else:
            print(f"Could not generate 2D structure for {smiles_to_draw}")
    except ImportError:
        print("Skipping visualization example: RDKit drawing or Matplotlib not installed.")
        print("To enable, install: pip install rdkit-pypi matplotlib pillow")
