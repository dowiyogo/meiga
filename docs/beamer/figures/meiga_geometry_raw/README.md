# Raw Geant4/Meiga Geometry Output

This folder is populated by `docs/beamer/render_meiga_geometry.sh`.

The guaranteed geometry artifact is `g4_00.wrl`. The primary manual renderer
is Castle Model Viewer:

`/home/reriosto/castle-model-viewer/castle-model-viewer`

If `figures/meiga_geometry_view.png` is not generated automatically, render the
WRL with Castle Model Viewer or another VRML-capable viewer and export a
screenshot manually with that name to enable the optional Beamer slide.

For an exact manual camera, save the printed Viewpoint text to:

- `manual_viewpoint_geometry.wrl.txt`
- `manual_viewpoint_trajectories.wrl.txt`

If either file exists, `render_meiga_geometry.sh` uses it preferentially when
building the zoom WRLs.

## Castle Model Viewer

Manual geometry render:

```bash
/home/reriosto/castle-model-viewer/castle-model-viewer \
	g4_00.wrl \
	--geometry 1600x1200 \
	--anti-alias 4 \
	--screenshot 0 ../meiga_geometry_view.png
```

Manual trajectories render:

```bash
/home/reriosto/castle-model-viewer/castle-model-viewer \
	g4_trajectories.wrl \
	--geometry 1600x1200 \
	--anti-alias 4 \
	--screenshot 0 ../meiga_trajectories_view.png
```

## Manual camera workflow

1. Open `/home/reriosto/castle-model-viewer/castle-model-viewer docs/beamer/figures/meiga_geometry_raw/g4_00.wrl`.
2. Adjust the view with the mouse: left button rotates, middle button or
   Shift+left translates, right button or wheel zooms.
3. Use `Clipboard -> Print Current Camera`.
4. Paste the printed Viewpoint into `manual_viewpoint_geometry.wrl.txt`.
5. Repeat the same process for `g4_trajectories.wrl` and
   `manual_viewpoint_trajectories.wrl.txt`.

## AlmaLinux 9 / Blender headless

Suggested packages:

```bash
sudo dnf install -y mesa-libGL mesa-libGLU mesa-dri-drivers mesa-libEGL mesa-libgbm libglvnd-glx xorg-x11-server-Xvfb blender ImageMagick
```

Manual smoke test:

```bash
LIBGL_ALWAYS_SOFTWARE=1 blender --background --version
```

Manual render test:

```bash
LIBGL_ALWAYS_SOFTWARE=1 MESA_LOADER_DRIVER_OVERRIDE=llvmpipe \
xvfb-run -a -s "-screen 0 1600x1200x24" \
blender --background --python docs/beamer/render_wrl_with_blender.py -- \
--input docs/beamer/figures/meiga_geometry_raw/g4_00.wrl \
--output docs/beamer/figures/meiga_geometry_view.png
```
