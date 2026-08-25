# openfe-failure-cases

This repository collects OpenFE simulation failure cases so we can better understand what fails, why it fails, and under what conditions. The goal is to support protocol improvements by sharing reproducible examples, diagnostics, and context around crashes, `NaN` results, and large repeat-to-repeat variation.

## What counts as a failure

In this repository, "failure" means a simulation that:
- crashes
- returns `NaN`
- shows large variation between repeats, typically around `~5 kcal/mol` or more

This does **not** mean a result is scientifically incorrect just because it differs from a reference value.

## Submission

To make a submission, please fork this repository and create a pull request with your case(s), and fill in the submission template. Multiple cases can be submitted in a single pull request, but please submit each case in its own folder containing the transformation JSON, logs and any other output files.

## What to include

Please submit cases with as much context as possible, including:
- the OpenFE transformation JSON
- environment details
- hardware details
- any diagnostics already run
- For `NaN` failures the xml file of the last state before the crash, stored in `nan_error_logs`

## Example submission

You can find an example of a failure case submission here: https://github.com/OpenFreeEnergy/openfe-failure-cases/pull/2


## Troubleshooting

### Periodic box size / nonbonded cutoff errors

You may see errors like:

- `openmm.OpenMMException: The periodic box size has decreased to less than twice the nonbonded cutoff.`
- `openmm.OpenMMException: NonbondedForce: The cutoff distance cannot be greater than half the periodic box size.`

These errors have the same cause: the system is too small for the chosen nonbonded cutoff.

#### Possible causes
- insufficient solvent padding
- incorrect periodic box vectors
- a solvated system that does not leave enough space around the solute

#### Suggested fixes
- If the system is solvated by the OpenFE protocol, increase the protocol's `solvent_padding` setting.
- If the system is explicitly solvated, verify that the correct box vectors were supplied.
- Ensure there is sufficient solvent around the system.

We do **not** recommend changing the nonbonded cutoff, since force fields are typically parameterized for these values and changing them may affect accuracy.

### SMIRNOFF force field parameter assignment errors

You may see errors like:

- `openff.interchange.exceptions.UnassignedBondError: BondHandler was not able to find parameters for the following valence terms: - Topology indices (0, 1): names and elements (C1 C), (Si1 Si),`

This usually means the SMIRNOFF-style force field does not have parameters for part of the chemistry in the system. In the example above, the force field was unable to assign bond parameters for bonds involving silicon.

#### Possible causes
- the system contains atoms or functional groups that are not covered by the selected force field

#### Suggested fixes
- check whether the chemistry in the system is expected to be supported by the chosen force field
- consider using a molecule specific force field with custom parameters to guarantee coverage of the chemistry in your system

If the error lists specific bonds or valence terms, those terms are the ones that could not be assigned parameters.

### NaN errors during simulation

You may see errors like:

- `openmm.OpenMMException: Particle coordinate is NaN.`
- `openmmtools.multistate.utils.SimulationNaNError: Propagating replica 0 at state 10 resulted in a NaN!`

For more information on the OpenMM error, see the OpenMM FAQ entry on NaNs:
https://github.com/openmm/openmm/wiki/Frequently-Asked-Questions#nan

These errors usually mean that the simulation became numerically unstable during minimization or propagation.

#### Possible causes
- missing capping groups or other issues in the receptor structure
- a poor atom mapping, including mappings that break bonds, map too few heavy atoms or have a large number of atoms whos hybridization changes
- an initial clash between the ligand and a crystal water, the receptor, or another part of the system that the minimizer could not relax
- a poor input structure for the ligand that the minimizer could not fix

#### Suggested fixes
- run `scripts/validate_transformation.py` on the transformation JSON to try to identify the source of the problem
- inspect the receptor for missing residues or missing capping groups
- review the atom mapping for bond-breaking transformations or mappings with too few heavy atoms, try to improve the mapping by aligning the ligands before generating the mapping
  - consider running the transformation using the [SepTop](https://docs.openfree.energy/en/latest/guide/protocols/septop.html) protocol which does not require an atom mapping and is robust to poor ligand alignment
- check the input ligand and starting pose for clashes or poor geometry, try and relax the ligand in the receptor before running the protocol
- if possible, rebuild or re-prep the input structures before rerunning the protocol

The state of the system and integrator before the error are often saved in a `nan-error-logs` directory, which can help with debugging.