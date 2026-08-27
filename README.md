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

## Local analysis

We provide [scripts](scripts) for local analysis of failure cases. These scripts can be helpful to diagnose private failure cases which cannot be submitted to this public repository for analysis by the OpenFE team.

The scripts help validate transformations and check for common issues, such as bond-breaking RBFE mappings or missing protein residues without caps and should help identify possible causes of failure for further investigation.

The scripts work with any OpenFE alchemical protocol including the `RelativeHybridTopologyProtocol`, `SepTopProtocol` and `AbsoluteBindingProtocol` though the atom mapping analysis is only relevant to the `RelativeHybridTopologyProtocol`.

The scripts are not exhaustive however so if you have a failure case that is not explained by the scripts, please submit it to the repository for further analysis.

## Required environment

These failure cases are intended to be analyzed with an OpenFE environment that includes:

- `openfe` (version 1.8 or newer, ideally the same version used to generate the Transformation JSON)
- `posebusters` (installable via `pip install posebusters`)

### Running a script

Each script takes a transformation JSON file as input. For example:

```bash
python scripts/validate_transformation.py path/to/transformation.json
```