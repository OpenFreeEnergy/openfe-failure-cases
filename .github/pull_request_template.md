# OpenFE Failure Case Submission

## Failure type
- [ ] Crash
- [ ] NaN result
- [ ] Large variation between repeats (~5 kcal/mol or more)
- [ ] Other:

## Notes
Add any extra context that could help reproduce or understand the failure.

- Did you try to run the same transformation using a different OpenFE protocol (SepTop or plainMD)?
- Did you try to run the same transformation using a different atom mapping method?
- Did you try to run the same transformation using a different force field?
- Where automatic protocol restarts enabled (by default this is true of openfe quickrun)? 
- Does the simulation crash immediately or after some progress?
- Do repeats crash at the same point in the simulation?
- How were the input structures prepared? (e.g., docking, co-folding, etc. provide as much detail as possible)

## Attachments
- [ ] Protocol Transformation JSON(s)
- [ ] Logs
- [ ] For `NaN` results the `state.xml` file stored in the `nan-error-logs` folder
- [ ] Environment YAML (conda list output)
- [ ] Hardware details

## Environment details
- **OpenMM accelerator details:** (CPU / GPU, CUDA version, etc.)
- **OS:**
