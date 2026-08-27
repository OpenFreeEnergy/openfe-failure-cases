import argparse
from gufe import ProteinComponent, ProteinMembraneComponent, Transformation, LigandAtomMapping, SmallMoleculeComponent
from gufe.visualization.mapping_visualization import draw_mapping
import numpy as np
from openff.units import unit as offunit
from openff.units import Quantity
from openmm import unit as omm_unit
from rdkit import Chem
import pathlib
import tempfile
from rdkit.Chem import Draw


def run_posebusters_validation(transformation: Transformation) -> dict[str, list[str]]:
    """
    Run Posebusters validation on the provided transformation.

    Parameters
    ----------
    transformation : Transformation
        The transformation to validate.

    Notes
    -----
    This function will attempt to run the Posebusters validation on the provided transformation.
    If Posebusters is not installed, this function will print a warning and skip the validation.

    """
    try:
        from posebusters import PoseBusters
        print("posebusters_avilable: True\n")
    except ImportError:
        print("posebusters_avilable: False\n")
        return {}

    # get the receptor if we have one
    pcs = transformation.stateA.get_components_of_type(ProteinComponent)
    receptor = None
    temp_receptor_path: pathlib.Path | None = None
    if pcs:
        # use the file path input as water clash detection is missed if we use rdkit.
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
            temp_receptor_path = pathlib.Path(tmp.name)
        pcs[0].to_pdb_file(temp_receptor_path.as_posix())
        receptor = temp_receptor_path.as_posix()
        # use the dock as we have no ground truth for the ligands
        buster_method = "dock"
    else:
        # if we have no receptor just bust the ligands
        buster_method = "mol"

    try:
        buster = PoseBusters(config=buster_method)

        pb_fails = {}

        # get a unique set of smcs across the end states this should work for all protocols
        smcs = {*transformation.stateA.get_components_of_type(SmallMoleculeComponent),
                *transformation.stateB.get_components_of_type(SmallMoleculeComponent)}
        for smc in smcs:
            df = buster.bust(mol_pred=smc.to_rdkit(), mol_cond=receptor, mol_true=None)
            pb_fails[smc.name] = []
            for _, row in df.iterrows():
                data = row.to_dict()
                for key, value in data.items():
                    if not value:
                        pb_fails[smc.name].append(key)
        return pb_fails
    finally:
        if temp_receptor_path is not None:
            temp_receptor_path.unlink(missing_ok=True)


def check_for_bond_breaks(mapping: LigandAtomMapping):
    """
    Checks for bond breaking/forming atom mappings in the provided LigandAtomMapping.

    Parameters
    ----------
    mapping : LigandAtomMapping
        The mapping to check

    Notes
    -----
    Bond breaking is detected by checking that if atoms are bonded in one state and mapped, that their mapped
    counterparts are also bonded in the other state.

    Raises
    ------
    ValueError
        If any bonds would be broken/introduced via the provided mapping.
    """
    # generate a list of bonds in the end states
    mol_a_bonds = {
        frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        for bond in mapping.componentA.to_rdkit().GetBonds()
    }
    mol_b_bonds = {
        frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        for bond in mapping.componentB.to_rdkit().GetBonds()
    }

    # build the list of directions to check
    checks = (
        # Bond in A is broken in B -> Broken
        ("A", "B", mol_a_bonds, mol_b_bonds, mapping.componentA_to_componentB),
        # Bond in B is missing in A -> Introduced
        ("B", "A", mol_b_bonds, mol_a_bonds, mapping.componentB_to_componentA),
    )

    for state_from, state_to, bonds_from, bonds_to, atom_map in checks:
        for atom_a, atom_b in bonds_from:
            if atom_a in atom_map and atom_b in atom_map:
                mapped_bond = frozenset((atom_map[atom_a], atom_map[atom_b]))
                if mapped_bond not in bonds_to:
                    # if bonded atoms are mapped but not bonded in the other endstate raise an error
                    raise ValueError(
                        f"Bond {atom_a}-{atom_b} in component{state_from} is broken in component{state_to} via the provided mapping."
                    )

def check_empty_mapping(mapping: LigandAtomMapping):
    if not mapping.componentA_to_componentB:
        raise ValueError("No atoms are mapped between the two alchemical components.")


def check_minimum_mapping(mapping: LigandAtomMapping):
    """
    Validates the provided mapping meets the minimum number of mapped heavy atoms (4) rule based on the ``mncar_score`` in lomap.

    Parameters
    ----------
    mapping : LigandAtomMapping
        The mapping between transforming components.

    Raises
    ------
    ValueError
        * If the atom mapping has less than 4 mapped heavy atoms and each component has more than 6 heavy atoms.

    Notes
    -----
    The two components must contain at least 4 mapped heavy atoms to pass this validation.
    If the components contain fewer than 6 heavy atoms in total then only a single mapped heavy atom is required.

    See Also
    --------
    lomap.gufe_bindings.scorers.mncar_score
    """
    atom_map = mapping.componentA_to_componentB
    mol_a = mapping.componentA.to_rdkit()
    mol_b = mapping.componentB.to_rdkit()

    mapped_heavy_atoms = dict()
    for atom_idx_a, atom_idx_b in atom_map.items():
        atom_a = mol_a.GetAtomWithIdx(atom_idx_a)
        atom_b = mol_b.GetAtomWithIdx(atom_idx_b)

        if atom_a.GetAtomicNum() != 1 and atom_b.GetAtomicNum() != 1:
            mapped_heavy_atoms[atom_idx_a] = atom_idx_b

    num_heavy_mol_a = mapping.componentA.to_rdkit().GetNumHeavyAtoms()
    num_heavy_mol_b = mapping.componentB.to_rdkit().GetNumHeavyAtoms()

    num_mapped_heavy_atoms = len(mapped_heavy_atoms)

    passed = (
        (num_mapped_heavy_atoms >= 4)
        or (num_heavy_mol_a < 6 and num_mapped_heavy_atoms >= 1)
        or (num_heavy_mol_b < 6 and num_mapped_heavy_atoms >= 1)
    )

    if not passed:
        raise ValueError(
            f"Number of mapped heavy atoms is less than 4: {num_mapped_heavy_atoms} which is required for this protocol."
        )


def check_for_missing_residues(protein: ProteinComponent | ProteinMembraneComponent, peptide_bond_cutoff: float = 3.0 * offunit.angstrom, box_vectors: Quantity | None = None):
    """
    Check for missing residues or capping groups by detecting large peptide bond distances.

    Parameters
    ----------
    peptide_bond_cutoff : Quantity
        The cutoff used to detect large peptide bond distances with units.
    box_vectors : Quantity, default None
        Periodic box vectors with units of length, compatible with
        nanometers. Must be a (3, 3) array in reduced form.

    Notes
    -----
    * We only check single bonded C-N inter-residue distances for this check, as these are the most indicative of missing residues or capping groups.
    * PBC can optionaly be used during the distance calculation.

    Raises
    ------
    ComponentValidationError
        If any inter-residue peptide C-N bonds are found to be longer than the specified cutoff,
        indicating likely uncapped or missing residues.
    """
    rd_mol = protein._rdkit
    conf = np.asarray(protein.to_openmm_positions().value_in_unit(omm_unit.nanometer))

    bond_threshold = peptide_bond_cutoff.m_as(offunit.nanometer)

    # Look for peptide bonds, these always involve C (residue i) - N (residue i+1)
    # We walk through atoms first to avoid performance issues with rd_mol.GetBonds()
    candidates = []

    for atom in rd_mol.GetAtoms():
        mi1 = atom.GetMonomerInfo()

        # Skip if we don't have monomer info or the wrong atom
        if (mi1 is None) or (mi1.GetName().strip() != "C"):
            continue

        for bond in atom.GetBonds():
            # Only want single bonds
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue

            other = bond.GetOtherAtom(atom)
            mi2 = other.GetMonomerInfo()

            if (mi2 is None) or (mi2.GetName().strip() != "N"):
                continue

            if mi1.GetChainId() != mi2.GetChainId():
                continue

            # Check if we have the same residue AND the same insertion code
            if (mi1.GetResidueNumber() == mi2.GetResidueNumber()) and (
                mi1.GetInsertionCode() == mi2.GetInsertionCode()
            ):
                continue

            candidates.append((mi1, mi2, atom.GetIdx(), other.GetIdx()))

    if not candidates:
        return

    if box_vectors is None:
        idx_i = np.array([c[2] for c in candidates])
        idx_j = np.array([c[3] for c in candidates])
        distances = np.linalg.norm(conf[idx_i] - conf[idx_j], axis=1)
    else:
        from openmm.app.internal import compiled
        # Use OpenMM's compiled functions to compute distances with PBC
        periodic_distance_func = compiled.periodicDistance(box_vectors.m_as(offunit.nanometer))
        distances = [periodic_distance_func(conf[c[2]], conf[c[3]]) for c in candidates]

    possible_bad_bonds = [(m1, m2, d * offunit.nanometer) for (m1, m2, _, _), d in zip(candidates, distances) if d > bond_threshold]

    if possible_bad_bonds:
        msg = "\n".join(
            f"{m1.GetChainId()}:{m1.GetResidueName()}{m1.GetResidueNumber()}:{m1.GetName().strip()} - "
            f"{m2.GetChainId()}:{m2.GetResidueName()}{m2.GetResidueNumber()}:{m2.GetName().strip()} = {d.m_as(offunit.angstrom):.2f} A"
            for m1, m2, d in possible_bad_bonds
        )
        raise ValueError(
            "Detected long inter-residue peptide C-N bonds, likely uncapped/missing residues. Check the following bonds:\n"
            + msg
        )


def main(transformation_file: str, write_local_files: bool = False):
    """
    Validate the provided transformation file for common issues.

    The following checks are performed:
    1. Check for missing residues in protein components.
    2. Check for empty atom mappings.
    3. Check for bond breaks in the mapping.
    4. Check for minimum number of mapped heavy atoms (at least 4).
    5. Run posebusters external validation if available.

    The results are writen to a local dir when requested with the same name as the input transformation file and include:
    - `receptor.pdb` the receptor used as input including waters
    - `ligands.sdf` the alchemical end state ligand in there input geometries
    - `mapping.png` the atom mapping for the transformation
    - `errors.log` the errors found by the script which are also printed to the terminal.
    """
    print(f"Validating transformation file: {transformation_file}")
    transformation = Transformation.from_json(transformation_file)
    tf_file = pathlib.Path(transformation_file)

    output_dir = pathlib.Path(f"{tf_file.stem}_validation_outputs") if write_local_files else None
    if output_dir is not None:
        output_dir.mkdir(exist_ok=True, parents=True)

    # get the protocol type
    print(f"Protocol type: {transformation.protocol.__class__.__name__}\n")

    # run the validation checks on the transformation
    # Check protein for missing residues if present
    pcs = transformation.stateA.get_components_of_type(ProteinComponent)

    tf_errors = []

    if pcs:
        protein = pcs[0]
        if isinstance(protein, ProteinMembraneComponent):
            box_vectors = protein.box_vectors
        else:
            box_vectors = None
        try:
            check_for_missing_residues(protein, peptide_bond_cutoff=3.0 * offunit.angstrom, box_vectors=box_vectors)
        except ValueError as e:
            tf_errors.append(f"Missing residues: {e}")

    # check the mapping
    mapping = transformation.mapping
    if mapping is not None:
        if isinstance(mapping, list):
            mapping = mapping[0]
        try:
            check_empty_mapping(mapping)
        except ValueError as e:
            tf_errors.append(f"Empty mapping: {e}")

        try:
            check_for_bond_breaks(mapping)
        except ValueError as e:
            tf_errors.append(f"Bond breaks: {e}")

        try:
            check_minimum_mapping(mapping)
        except ValueError as e:
            tf_errors.append(f"Minimum mapping: {e}")

    pb_fails = run_posebusters_validation(transformation)
    for key, errors in pb_fails.items():
        if errors:
            tf_errors.append(f"Posebusters validation failures for {key}: {', '.join(errors)}")

    if tf_errors:
        # print the errors to the terminal and write them to a log file
        header = f"VALIDATION SUMMARY for: {transformation.name}:"
        message = f"{header}\n" + "\n".join(f"  {error}" for error in tf_errors) + "\n"

        print(message, end="")

        if output_dir is not None:
            with (output_dir / "errors.log").open("w") as out:
                out.write(message)

    else:
        print(f"Transformation {transformation.name} passed all validation checks.")

    # extract the inputs for visualisation
    if output_dir is not None:
        print(f"SCRIPT OUTPUT SUMMARY GENERATED at: {output_dir}")
        print("  - receptor.pdb: the receptor used as input including waters (if present)")
        print("  - ligands.sdf: the alchemical end state ligand in there input geometries")
        print("  - mapping.svg: the atom mapping for the transformation")
        print("  - errors.log: the errors found by the script which are also printed to the terminal.")

        smcs = {*transformation.stateA.get_components_of_type(SmallMoleculeComponent),
                *transformation.stateB.get_components_of_type(SmallMoleculeComponent)}

        supplier = Chem.SDWriter((output_dir / "ligands.sdf").as_posix())
        for smc in smcs:
            supplier.write(smc.to_rdkit())

        pcs = transformation.stateA.get_components_of_type(ProteinComponent)
        if pcs:
            pcs[0].to_pdb_file((output_dir / "receptor.pdb").as_posix())

        # write out the mapping to svg for visualisation
        if mapping is not None:
            d2d = Draw.rdMolDraw2D.MolDraw2DSVG(600, 600, 300, 300)
            svg_text = draw_mapping(mol1_to_mol2=mapping.componentA_to_componentB, mol1=mapping.componentA.to_rdkit(), mol2=mapping.componentB.to_rdkit(), d2d=d2d)
            # Keep only the first complete SVG document, without this there is an error displayed with the mapping
            start = svg_text.find("<svg")
            end = svg_text.find("</svg>", start)
            if start == -1 or end == -1:
                raise ValueError("draw_mapping did not return a valid SVG document.")
            svg_text = svg_text[start:end + len("</svg>")]

            with (output_dir / "mapping.svg").open("w", encoding="utf-8") as out:
                out.write(svg_text)
    else:
        print("SCRIPT OUTPUT SUMMARY: enable local file writing with --write-local-files for additional visualisation of the Transformation.")


def cli():
    """Command line interface for validating a transformation file."""
    parser = argparse.ArgumentParser(
        description="Validate a transformation file for common issues."
    )
    parser.add_argument(
        "transformation_file",
        type=str,
        help="Path to the transformation JSON file to validate.",
    )
    parser.add_argument(
        "--write-local-files",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Write receptor.pdb, ligands.sdf, mapping.svg, and errors.log locally for visualisation and debugging, by default this is disabled.",
    )
    args = parser.parse_args()
    main(args.transformation_file, write_local_files=args.write_local_files)


if __name__ == "__main__":
    cli()
