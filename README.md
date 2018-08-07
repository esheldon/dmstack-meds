# medsdm

Code to create MEDS files from a LSST dmstack interface

## Example

The following command will process a particular tract, patch of DC2 run 1.1p:
```
$ python -m medsdm.maker --config test/config_DC2.yaml /global/projecta/projectdirs/lsst/global/in2p3/Run1.1/output 4638 0,0 i
```
