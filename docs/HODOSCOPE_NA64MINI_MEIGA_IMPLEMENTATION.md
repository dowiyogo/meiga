# NA64-mini / Muography Hodoscope in Meiga

This implementation adds an explicit `NA64Mini` mode to `eHodoscope`.

## Geometry

- Two X-Y detector modules are built by one logical hodoscope detector.
- Lower module: `z_min = 0.00 mm`, `z_max = 15.00 mm`, `z_center = 7.50 mm`.
- Upper module: `z_min = 98.40 mm`, `z_max = 113.40 mm`, `z_center = 105.90 mm`.
- Equivalent center separation: `98.40 mm`.
- Free gap between module faces: `83.40 mm`.
- Module thickness: `15.00 mm`.
- Active material mode defaults to `active_only`.

Each module has two axes, `X` and `Y`. Each axis has 16 resolvable channels:

- Subplane A: channels `0..7`.
- Subplane B: channels `8..15`.
- Subplane offset: `2.00 mm`.
- Transverse centres in subplane A: `-15, -11, -7, -3, 1, 5, 9, 13 mm`.
- Transverse centres in subplane B: `-13, -9, -5, -1, 3, 7, 11, 15 mm`.
- Bar pitch: `4.00 mm`.
- Bar width: `3.00 mm`.
- Bar length: `75.00 mm`.
- Active bar thickness: `1.00 mm`.

The four active layer centers inside each module are:

- `X_A`: `z_min + 1.5 mm`.
- `X_B`: `z_min + 4.5 mm`.
- `Y_A`: `z_min + 10.5 mm`.
- `Y_B`: `z_min + 13.5 mm`.

The same-axis subplane centre spacing is `3.0 mm` (the bars have a `2.0 mm`
air gap). The X-block to Y-block centre spacing is `6.0 mm` (the inner bar
faces have a `5.0 mm` air gap). The module-centre separation is `98.4 mm`.

## Axis Convention

The code uses reconstruction-axis naming:

- `axis = "X"` is a channel used to reconstruct X. Its scintillator bar is long
  in Y and segmented in X.
- `axis = "Y"` is a channel used to reconstruct Y. Its scintillator bar is long
  in X and segmented in Y.

This avoids ambiguity between the physical long direction and the reconstructed
coordinate.

## Channels

The sensitive bar IDs are one-based for Geant4, but the hodoscope global channel
IDs are zero-based:

- `LOWER_X_00 ... LOWER_X_15`: global channels `0..15`.
- `LOWER_Y_00 ... LOWER_Y_15`: global channels `16..31`.
- `UPPER_X_00 ... UPPER_X_15`: global channels `32..47`.
- `UPPER_Y_00 ... UPPER_Y_15`: global channels `48..63`.

Total active bars: `64`.

## Pixels

Each module stores a discrete `16 x 16` pixel matrix:

```text
pixel_id = x_index + 16 * y_index
```

The overlap geometry does not create extra pixels. For events with multiple
hits in one axis, the primary pixel uses the hit channel with the largest energy
deposit, while the full channel list remains in the output.

Sub-pixel triangulation in overlap regions is intentionally left as a TODO.

## Materials

- Active material metadata: BC-408 or equivalent scintillating polystyrene.
- Light yield metadata: `12800 photons/MeV`.
- Attenuation length metadata: `2100 mm`.
- MPPC metadata: Hamamatsu `S12572-100P`.
- Lower module reflector metadata: `VikuitiESR`.
- Upper module reflector metadata: `TiO2Paint`.

Reflectors are not mixed into the active scintillator. In the default
`active_only` mode, reflector and dead-material descriptions are metadata only.
The provisional TiO2 paint thickness for later passive material modeling is
`100 um` (ASSUMPTION).

## Scoring

The sensitive detector now stores, per event:

- active hits by bar above `active_hit_threshold_MeV = 0.1`;
- deposited energy per bar;
- first hit time per bar when available;
- `LOWER_X`, `LOWER_Y`, `UPPER_X`, `UPPER_Y` hit flags;
- lower and upper X-Y coincidences;
- full four-plane coincidence;
- lower and upper discrete pixel IDs.

## Diagnostics

Run:

```bash
diagnostics/run_hodoscope_geometry_check.sh
```

This builds `G4HodoscopeSimulator`, runs 1000 events with `CheckOverlaps = true`,
and writes:

- `diagnostics/hodoscope_geometry_check.log`
- `diagnostics/hodoscope_geometry_summary.md`
- `diagnostics/hodoscope_output_summary.md`
- `diagnostics/hodoscope_bar_occupancy.csv`
- `diagnostics/hodoscope_pixel_occupancy_lower.csv`
- `diagnostics/hodoscope_pixel_occupancy_upper.csv`
- `diagnostics/hodoscope_coincidence_summary.md`
- `diagnostics/hodoscope_coincidences.csv`

## TODO

- Iteration 2: simplified dead material.
- Iteration 3: parameterized optical response.
- Iteration 4: SiPM pulses and electronics.
- Iteration 5: angular reconstruction and F matrix.
