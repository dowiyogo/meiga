# NA64-mini / Muography Hodoscope Geometry Summary

- number of modules: 2
- lower module:
  - z_min = 0.00 mm
  - z_max = 12.00 mm
  - z_center = 6.00 mm
- upper module:
  - z_min = 98.40 mm
  - z_max = 110.40 mm
  - z_center = 104.40 mm
- center separation = 98.40 mm
- free separation between modules = 86.40 mm
- bars per axis = 16
- bars per module = 32
- total active bars = 64
- pixels per module = 16 x 16 = 256
- subplanes per axis = 2
- bars per subplane = 8
- subplane offset = 2.00 mm
- bar length = 75.00 mm
- bar width = 3.00 mm
- active thickness = 1.00 mm
- X-Y separation D = 5.00 mm between the inner X_B and Y_A layer centers
- lower module reflector = VikuitiESR
- upper module reflector = TiO2Paint
- active material = BC408 or equivalent scintillating polystyrene
- light yield metadata = 12800 photons/MeV
- attenuation length metadata = 2100 mm
- MPPC metadata = Hamamatsu S12572-100P
- dead material mode = active_only
- CheckOverlaps result = PASS: no severe overlap message found

## Channel Convention

`axis = "X"` denotes a channel used to reconstruct X. Its bar is geometrically
long along Y and segmented in X. `axis = "Y"` denotes a channel used to
reconstruct Y; its bar is geometrically long along X and segmented in Y.

## TODO

- Iteration 2: simplified dead material.
- Iteration 3: parameterized optical response.
- Iteration 4: SiPM pulses and electronics.
- Iteration 5: angular reconstruction and F matrix.
