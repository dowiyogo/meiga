#!/usr/bin/env bash

set -u

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
raw_dir="${script_dir}/figures/meiga_geometry_raw"
figures_dir="${script_dir}/figures"
app_dir="${repo_root}/build/Applications/G4HodoscopeSimulator"
src_app_dir="${repo_root}/src/Applications/G4HodoscopeSimulator"
log_file="${raw_dir}/render_meiga_geometry.log"
environment_file="${raw_dir}/render_environment.md"
status_file="${raw_dir}/render_status.md"
castle_model_viewer_default="/home/reriosto/castle-model-viewer/castle-model-viewer"
castle_model_converter_default="/home/reriosto/castle-model-viewer/castle-model-converter"
castle_model_viewer_exe=""
castle_model_converter_exe=""

mkdir -p "${raw_dir}" "${figures_dir}"

warn() {
	echo "WARNING: $*" | tee -a "${log_file}"
}

status_line() {
	local key="$1"
	local value="$2"
	printf -- '- %s: %s\n' "${key}" "${value}" >>"${status_file}"
}

run_logged() {
	echo "" >>"${log_file}"
	echo "\$ $*" >>"${log_file}"
	"$@" >>"${log_file}" 2>&1
}

find_windows_blender() {
	if [ -n "${WINDOWS_BLENDER_EXE-}" ] && [ -x "${WINDOWS_BLENDER_EXE}" ]; then
		printf '%s\n' "${WINDOWS_BLENDER_EXE}"
		return 0
	fi

	local candidate
	for candidate in \
		/mnt/c/Program\ Files/Blender\ Foundation/*/blender.exe \
		/mnt/c/Program\ Files/Blender*/blender.exe \
		/mnt/c/Users/*/AppData/Local/Programs/Blender\ Foundation/*/blender.exe; do
		if [ -x "${candidate}" ]; then
			printf '%s\n' "${candidate}"
			return 0
		fi
	done

	return 1
}

wsl_to_windows_path() {
	local path="$1"
	if command -v wslpath >/dev/null 2>&1; then
		wslpath -w "${path}"
		return 0
	fi
	printf '%s\n' "${path}"
}

windows_blender_version() {
	local blender_exe="$1"
	"${blender_exe}" --version 2>&1 || true
}

find_castle_model_viewer() {
	if [ -n "${CASTLE_MODEL_VIEWER-}" ] && [ -x "${CASTLE_MODEL_VIEWER}" ]; then
		printf '%s\n' "${CASTLE_MODEL_VIEWER}"
		return 0
	fi

	if [ -x "${castle_model_viewer_default}" ]; then
		printf '%s\n' "${castle_model_viewer_default}"
		return 0
	fi

	return 1
}

find_castle_model_converter() {
	if [ -n "${CASTLE_MODEL_CONVERTER-}" ] && [ -x "${CASTLE_MODEL_CONVERTER}" ]; then
		printf '%s\n' "${CASTLE_MODEL_CONVERTER}"
		return 0
	fi

	if [ -x "${castle_model_converter_default}" ]; then
		printf '%s\n' "${castle_model_converter_default}"
		return 0
	fi

	return 1
}

castle_viewer_help() {
	local viewer_exe="$1"
	"${viewer_exe}" --help 2>&1 || true
}

insert_viewpoint_into_wrl() {
	local input_wrl="$1"
	local output_wrl="$2"
	local viewpoint_text_file="$3"

	python3 - "$input_wrl" "$output_wrl" "$viewpoint_text_file" <<'PY' >>"${log_file}" 2>&1
import sys
from pathlib import Path

input_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
viewpoint_path = Path(sys.argv[3])

source = input_path.read_text(encoding="utf-8")
viewpoint_text = viewpoint_path.read_text(encoding="utf-8").strip()

if not source.startswith("#VRML V2.0 utf8"):
    raise SystemExit("input is not a VRML 2.0 file")

lines = source.splitlines(keepends=True)
insert_at = 1
while insert_at < len(lines) and lines[insert_at].lstrip().startswith("#"):
	insert_at += 1
while insert_at < len(lines) and lines[insert_at].strip() == "":
	insert_at += 1

rendered = "".join(lines[:insert_at])
if not rendered.endswith("\n"):
	rendered += "\n"
rendered += viewpoint_text + "\n\n"
rendered += "".join(lines[insert_at:])
output_path.write_text(rendered, encoding="utf-8")
PY
}

make_viewpoint_wrl() {
	local input_wrl="$1"
	local output_wrl="$2"
	local viewpoint_name="$3"
	local position="$4"
	local orientation="$5"
	local field_of_view="$6"
	local manual_viewpoint_file="$7"

	rm -f "${output_wrl}"
	if [ -n "${manual_viewpoint_file}" ] && [ -s "${manual_viewpoint_file}" ]; then
		insert_viewpoint_into_wrl "${input_wrl}" "${output_wrl}" "${manual_viewpoint_file}"
		return $?
	fi

	python3 "${script_dir}/add_wrl_viewpoint.py" \
		--input "${input_wrl}" \
		--output "${output_wrl}" \
		--name "${viewpoint_name}" \
		--position "${position}" \
		--orientation "${orientation}" \
		--field-of-view "${field_of_view}" >>"${log_file}" 2>&1
}

render_castle_viewpoint_png() {
	local input_wrl="$1"
	local output_png="$2"
	local viewpoint_name="$3"
	local label="$4"

	if [ ! -f "${input_wrl}" ]; then
		warn "${label}: input file not found: ${input_wrl}"
		return 1
	fi

	if render_with_castle_model_viewer "${input_wrl}" "${output_png}" "${viewpoint_name}"; then
		status_line "${label}" "${viewpoint_name}"
		return 0
	fi

	return 1
}

render_castle_viewpoint_png_with_viewpoint() {
	local input_wrl="$1"
	local output_png="$2"
	local viewpoint_name="$3"
	local label="$4"
	local render_wrl="${input_wrl}"

	if render_with_castle_model_viewer "${input_wrl}" "${output_png}" "${viewpoint_name}"; then
		status_line "${label}" "${viewpoint_name}"
		return 0
	fi

	return 1
}

render_castle_wrl_with_viewpoint() {
	local input_wrl="$1"
	local output_png="$2"
	local viewpoint_name="$3"
	local label="$4"
	local rendered=0

	if [ ! -f "${input_wrl}" ]; then
		warn "${label}: input file not found: ${input_wrl}"
		return 1
	fi

	render_with_castle_model_viewer "${input_wrl}" "${output_png}" "${viewpoint_name}" || rendered=1
	return ${rendered}
}

render_zoom_variant() {
	local base_wrl="$1"
	local variant_wrl="$2"
	local variant_png="$3"
	local viewpoint_name="$4"
	local position="$5"
	local orientation="$6"
	local field_of_view="$7"
	local manual_viewpoint_file="$8"
	local label="$9"

	make_viewpoint_wrl "${base_wrl}" "${variant_wrl}" "${viewpoint_name}" "${position}" "${orientation}" "${field_of_view}" "${manual_viewpoint_file}"
	render_with_castle_model_viewer "${variant_wrl}" "${variant_png}" "${viewpoint_name}"
	if file_nonempty "${variant_png}"; then
		status_line "${label}" "éxito (${viewpoint_name})"
		return 0
	fi
	status_line "${label}" "fallo (${viewpoint_name})"
	return 1
}

write_environment() {
	castle_model_viewer_exe="$(find_castle_model_viewer || true)"
	castle_model_converter_exe="$(find_castle_model_converter || true)"
	{
		echo "# Render Environment"
		echo
		echo '```text'
		echo '$ /home/reriosto/castle-model-viewer/castle-model-viewer --help || true'
		/home/reriosto/castle-model-viewer/castle-model-viewer --help || true
		echo
		echo '$ /home/reriosto/castle-model-viewer/castle-model-converter --help || true'
		/home/reriosto/castle-model-viewer/castle-model-converter --help || true
		echo
		echo '$ echo ${CASTLE_MODEL_VIEWER-}'
		echo "${CASTLE_MODEL_VIEWER-}"
		echo
		echo '$ echo ${CASTLE_MODEL_CONVERTER-}'
		echo "${CASTLE_MODEL_CONVERTER-}"
		echo
		echo '$ CASTLE_MODEL_VIEWER_DEFAULT'
		echo "${castle_model_viewer_default}"
		echo
		echo '$ CASTLE_MODEL_CONVERTER_DEFAULT'
		echo "${castle_model_converter_default}"
		echo
		echo '$ CASTLE_MODEL_VIEWER_EXE'
		echo "${castle_model_viewer_exe:-not found}"
		echo
		echo '$ CASTLE_MODEL_CONVERTER_EXE'
		echo "${castle_model_converter_exe:-not found}"
		echo
		echo '$ uname -a'
		uname -a || true
		echo
		echo '$ cat /etc/os-release || true'
		cat /etc/os-release || true
		echo
		echo '$ which blender || true'
		which blender || true
		echo
		echo '$ which wslpath || true'
		which wslpath || true
		echo
		echo '$ ls /mnt/c/Program\ Files/Blender\ Foundation || true'
		ls /mnt/c/Program\ Files/Blender\ Foundation 2>/dev/null || true
		echo
		echo '$ echo ${WINDOWS_BLENDER_EXE-}'
		echo "${WINDOWS_BLENDER_EXE-}"
		echo
		windows_blender_exe="$(find_windows_blender || true)"
		if [ -n "${windows_blender_exe}" ]; then
			echo '$ WINDOWS_BLENDER_EXE=... --version'
			windows_blender_version "${windows_blender_exe}"
		else
			echo '$ WINDOWS_BLENDER_EXE=... --version'
			echo "Windows Blender not detected"
		fi
		echo
		echo '$ blender --version || true'
		blender --version || true
		echo
		echo '$ which xvfb-run || true'
		which xvfb-run || true
		echo
		echo '$ which glxinfo || true'
		which glxinfo || true
		echo
		echo '$ which glxgears || true'
		which glxgears || true
		echo
		echo '$ which python3 || true'
		which python3 || true
		echo
		echo '$ which magick || true'
		which magick || true
		echo
		echo '$ which convert || true'
		which convert || true
		echo
		echo '$ rpm -qa | grep -E "mesa|glvnd|blender|xorg-x11-server-Xvfb|ImageMagick" | sort || true'
		rpm -qa | grep -E "mesa|glvnd|blender|xorg-x11-server-Xvfb|ImageMagick" | sort || true
		echo
		echo '$ echo $DISPLAY'
		echo "${DISPLAY-}"
		echo
		echo '$ echo $LIBGL_ALWAYS_SOFTWARE'
		echo "${LIBGL_ALWAYS_SOFTWARE-}"
		echo
		echo '$ echo $MESA_LOADER_DRIVER_OVERRIDE'
		echo "${MESA_LOADER_DRIVER_OVERRIDE-}"
		echo
		echo '$ LIBGL_ALWAYS_SOFTWARE=1 glxinfo -B || true'
		if command -v glxinfo >/dev/null 2>&1; then
			LIBGL_ALWAYS_SOFTWARE=1 glxinfo -B || true
		else
			echo "glxinfo not found"
		fi
		echo '```'
	} >"${environment_file}" 2>&1
}

write_geometry_config() {
	local output_json="$1"
	local render_file="$2"
	local geovis="$3"
	local trajvis="$4"
	local n_particles="$5"
	local output_data="$6"
	python3 - "${src_app_dir}/G4HodoscopeSimulator.json" "${output_json}" "${render_file}" "${geovis}" "${trajvis}" "${n_particles}" "${output_data}" <<'PY' >>"${log_file}" 2>&1
import json
import sys
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
render_file = sys.argv[3]
geovis = sys.argv[4].lower() == "true"
trajvis = sys.argv[5].lower() == "true"
n_particles = int(sys.argv[6])
output_data = sys.argv[7]

data = json.loads(src.read_text(encoding="utf-8"))
data.setdefault("Input", {})["InputNParticles"] = n_particles
data.setdefault("Output", {})["OutputFile"] = output_data
data["Output"]["CompressOutput"] = False
data["DetectorList"] = "./DetectorList.xml"
data.setdefault("Simulation", {})["GeoVisOn"] = geovis
data["Simulation"]["TrajVisOn"] = trajvis
data["Simulation"]["CheckOverlaps"] = True
data["Simulation"]["RenderFile"] = render_file
data["Simulation"]["Verbosity"] = 0
dst.write_text(json.dumps(data, indent=2), encoding="utf-8")
PY
}

run_meiga_config() {
	local config_file="$1"
	local label="$2"
	if [ ! -x "${app_dir}/G4HodoscopeSimulator" ]; then
		warn "G4HodoscopeSimulator executable not found at ${app_dir}; cannot run ${label}."
		return 1
	fi
	(
		cd "${raw_dir}" || exit 1
		"${app_dir}/G4HodoscopeSimulator" -c "./$(basename "${config_file}")"
	) >>"${log_file}" 2>&1
}

file_nonempty() {
	[ -s "$1" ]
}

try_windows_blender_attempt() {
	local attempt_name="$1"
	local input_file="$2"
	local output_file="$3"
	local blender_script="$4"
	local blender_exe
	local input_windows
	local output_windows
	local script_windows

	if [ ! -f "${input_file}" ]; then
		warn "${attempt_name}: input file not found: ${input_file}"
		status_line "${attempt_name}" "fallo (WRL no existe)"
		return 1
	fi

	blender_exe="$(find_windows_blender || true)"
	if [ -z "${blender_exe}" ]; then
		warn "${attempt_name}: Windows Blender not found"
		status_line "${attempt_name}" "fallo (Windows Blender no disponible)"
		return 1
	fi

	input_windows="$(wsl_to_windows_path "${input_file}")"
	output_windows="$(wsl_to_windows_path "${output_file}")"
	script_windows="$(wsl_to_windows_path "${blender_script}")"

	rm -f "${output_file}"
	echo "Trying ${attempt_name} / Windows Blender: ${blender_exe}" | tee -a "${log_file}"
	"${blender_exe}" --background --python "${script_windows}" -- \
		--input "${input_windows}" \
		--output "${output_windows}" >>"${log_file}" 2>&1 || true

	if file_nonempty "${output_file}"; then
		echo "${attempt_name}: rendered ${output_file}" | tee -a "${log_file}"
		status_line "${attempt_name} / Windows Blender" "éxito"
		return 0
	fi

	warn "Windows Blender was found but could not render WRL."
	status_line "${attempt_name} / Windows Blender" "fallo"
	return 1
}

try_blender_attempt() {
	local attempt_name="$1"
	local input_file="$2"
	local output_file="$3"
	shift 3

	if [ ! -f "${input_file}" ]; then
		warn "${attempt_name}: input file not found: ${input_file}"
		status_line "${attempt_name}" "fallo (WRL no existe)"
		return 1
	fi
	if ! command -v blender >/dev/null 2>&1; then
		warn "${attempt_name}: blender not found"
		status_line "${attempt_name}" "fallo (blender no disponible)"
		return 1
	fi

	rm -f "${output_file}"
	echo "Trying ${attempt_name}" | tee -a "${log_file}"
	"$@" >>"${log_file}" 2>&1

	if file_nonempty "${output_file}"; then
		echo "${attempt_name}: rendered ${output_file}" | tee -a "${log_file}"
		status_line "${attempt_name}" "éxito"
		return 0
	fi

	warn "${attempt_name}: no PNG produced"
	status_line "${attempt_name}" "fallo"
	return 1
}

render_wrl_with_blender_chain() {
	local input_file="$1"
	local output_file="$2"
	local label="$3"
	local script_path="${script_dir}/render_wrl_with_blender.py"

	if try_blender_attempt "${label} / Blender normal" "${input_file}" "${output_file}" \
		blender --background --python "${script_path}" -- \
			--input "${input_file}" \
			--output "${output_file}"; then
		return 0
	fi

	if try_blender_attempt "${label} / Blender software GL" "${input_file}" "${output_file}" \
		env LIBGL_ALWAYS_SOFTWARE=1 \
		blender --background --python "${script_path}" -- \
			--input "${input_file}" \
			--output "${output_file}"; then
		return 0
	fi

	if try_blender_attempt "${label} / Blender llvmpipe" "${input_file}" "${output_file}" \
		env LIBGL_ALWAYS_SOFTWARE=1 MESA_LOADER_DRIVER_OVERRIDE=llvmpipe \
		blender --background --python "${script_path}" -- \
			--input "${input_file}" \
			--output "${output_file}"; then
		return 0
	fi

	if command -v xvfb-run >/dev/null 2>&1; then
		if try_blender_attempt "${label} / Xvfb + llvmpipe" "${input_file}" "${output_file}" \
			env LIBGL_ALWAYS_SOFTWARE=1 MESA_LOADER_DRIVER_OVERRIDE=llvmpipe \
			xvfb-run -a -s "-screen 0 1600x1200x24" \
			blender --background --python "${script_path}" -- \
				--input "${input_file}" \
				--output "${output_file}"; then
			return 0
		fi
	else
		warn "xvfb-run is not available for ${label} fallback."
		status_line "${label} / Xvfb + llvmpipe" "fallo (xvfb-run no disponible)"
	fi

	if try_windows_blender_attempt "${label}" "${input_file}" "${output_file}" "${script_path}"; then
		return 0
	fi

	warn "Blender could not render ${label}; raw WRL remains available."
	return 1
}

render_with_castle_model_viewer() {
	local input_file="$1"
	local output_file="$2"
	local viewpoint_name="${3-}"
	local viewer_exe="${castle_model_viewer_exe}"
	local -a castle_command

	if [ -z "${viewer_exe}" ]; then
		viewer_exe="$(find_castle_model_viewer || true)"
	fi

	if [ -z "${viewer_exe}" ]; then
		warn "Castle Model Viewer not available"
		return 1
	fi

	if [ ! -f "${input_file}" ]; then
		warn "Castle Model Viewer: input file not found: ${input_file}"
		return 1
	fi

	rm -f "${output_file}"
	echo "Trying Castle Model Viewer / input-first: ${viewer_exe}" | tee -a "${log_file}"
	castle_command=("${viewer_exe}" "${input_file}" --geometry 1600x1200 --anti-alias 4 --screenshot 0 "${output_file}")
	if [ -n "${viewpoint_name}" ]; then
		castle_command=("${viewer_exe}" "${input_file}" --viewpoint "${viewpoint_name}" --geometry 1600x1200 --anti-alias 4 --screenshot 0 "${output_file}")
	fi
	"${castle_command[@]}" >>"${log_file}" 2>&1 || true

	if file_nonempty "${output_file}"; then
		echo "Castle Model Viewer: rendered ${output_file}" | tee -a "${log_file}"
		status_line "Castle Model Viewer / input-first" "éxito"
		return 0
	fi

	rm -f "${output_file}"
	echo "Trying Castle Model Viewer / flags-first: ${viewer_exe}" | tee -a "${log_file}"
	castle_command=("${viewer_exe}" --geometry 1600x1200 --anti-alias 4 --screenshot 0 "${output_file}" "${input_file}")
	if [ -n "${viewpoint_name}" ]; then
		castle_command=("${viewer_exe}" --viewpoint "${viewpoint_name}" --geometry 1600x1200 --anti-alias 4 --screenshot 0 "${output_file}" "${input_file}")
	fi
	"${castle_command[@]}" >>"${log_file}" 2>&1 || true

	if file_nonempty "${output_file}"; then
		echo "Castle Model Viewer: rendered ${output_file}" | tee -a "${log_file}"
		status_line "Castle Model Viewer / flags-first" "éxito"
		return 0
	fi

	warn "Castle Model Viewer could not render ${input_file}"
	status_line "Castle Model Viewer" "fallo"
	return 1
}

try_vtk_python_render() {
	local input_file="$1"
	local output_file="$2"

	if [ ! -f "${input_file}" ]; then
		warn "VTK Python: input file not found: ${input_file}"
		status_line "VTK Python" "fallo (WRL no existe)"
		return 1
	fi
	if ! command -v python3 >/dev/null 2>&1; then
		warn "VTK Python: python3 not available"
		status_line "VTK Python" "fallo (python3 no disponible)"
		return 1
	fi
	if ! python3 - <<'PY' >/dev/null 2>&1
try:
    import vtk  # noqa: F401
except Exception:
    raise SystemExit(1)
PY
	then
		warn "VTK Python module not available"
		status_line "VTK Python" "no disponible"
		return 1
	fi

	rm -f "${output_file}"
	python3 - "${input_file}" "${output_file}" >>"${log_file}" 2>&1 <<'PY' || true
import sys
from pathlib import Path

import vtk

input_wrl = Path(sys.argv[1])
output_png = Path(sys.argv[2])

importer = vtk.vtkVRMLImporter()
importer.SetFileName(str(input_wrl))
render_window = vtk.vtkRenderWindow()
importer.SetRenderWindow(render_window)
importer.Read()
render_window.OffScreenRenderingOn()
render_window.SetSize(1600, 1200)
render_window.Render()

renderer = render_window.GetRenderers().GetFirstRenderer()
if renderer is None:
    raise SystemExit(1)

renderer.ResetCamera()
render_window.Render()

window_to_image = vtk.vtkWindowToImageFilter()
window_to_image.SetInput(render_window)
window_to_image.Update()

writer = vtk.vtkPNGWriter()
writer.SetFileName(str(output_png))
writer.SetInputConnection(window_to_image.GetOutputPort())
if writer.Write() != 1:
    raise SystemExit(1)
PY

	if file_nonempty "${output_file}"; then
		echo "VTK Python: rendered ${output_file}" | tee -a "${log_file}"
		status_line "VTK Python" "éxito"
		return 0
	fi

	warn "VTK Python could not render ${input_file}"
	status_line "VTK Python" "fallo"
	return 1
}

render_wrl_with_castle_then_fallbacks() {
	local input_file="$1"
	local output_file="$2"
	local label="$3"

	if render_with_castle_model_viewer "${input_file}" "${output_file}"; then
		return 0
	fi

	if try_vtk_python_render "${input_file}" "${output_file}"; then
		return 0
	fi

	if render_wrl_with_blender_chain "${input_file}" "${output_file}" "${label}"; then
		return 0
	fi

	warn "${label}: Castle Model Viewer, VTK Python and Blender fallbacks failed; raw WRL remains available."
	return 1
}

convert_to_png() {
	local input_file="$1"
	local output_file="$2"
	if command -v magick >/dev/null 2>&1; then
		magick "${input_file}" "${output_file}" >>"${log_file}" 2>&1
		return $?
	fi
	if command -v convert >/dev/null 2>&1; then
		convert "${input_file}" "${output_file}" >>"${log_file}" 2>&1
		return $?
	fi
	return 1
}

try_raytracer() {
	local output_file="${figures_dir}/meiga_geometry_view.png"
	local ray_config="${raw_dir}/G4HodoscopeSimulator_raytracer.json"
	local candidates

	echo "Trying Geant4 RayTracer route" | tee -a "${log_file}"
	write_geometry_config "${ray_config}" "RayTracer" "true" "false" "1" "./raytracer_output.json"
	rm -f "${raw_dir}"/g4RayTracer.* "${raw_dir}"/raytrace.* "${raw_dir}"/RayTracer.*

	if ! run_meiga_config "${ray_config}" "RayTracer"; then
		warn "RayTracer run failed."
		status_line "RayTracer" "fallo/no disponible"
		return 1
	fi

	candidates="$(find "${raw_dir}" -maxdepth 1 \( -iname 'g4RayTracer*.jpeg' -o -iname 'g4RayTracer*.jpg' -o -iname 'g4RayTracer*.png' -o -iname 'raytrace*.jpeg' -o -iname 'raytrace*.jpg' -o -iname 'raytrace*.png' -o -iname 'RayTracer*.jpg' -o -iname 'RayTracer*.png' \) -type f | sort)"
	if [ -z "${candidates}" ]; then
		warn "RayTracer not available or not configured in this Geant4/Meiga build."
		status_line "RayTracer" "no disponible"
		return 1
	fi

	local first_candidate
	first_candidate="$(printf '%s\n' "${candidates}" | head -n 1)"
	case "${first_candidate,,}" in
		*.png)
			cp "${first_candidate}" "${output_file}" >>"${log_file}" 2>&1 || true
			;;
		*.jpg|*.jpeg)
			convert_to_png "${first_candidate}" "${output_file}" || true
			;;
	esac

	if file_nonempty "${output_file}"; then
		echo "RayTracer rendered ${output_file}" | tee -a "${log_file}"
		status_line "RayTracer" "éxito"
		return 0
	fi

	warn "RayTracer produced ${first_candidate}, but conversion/copy to PNG failed."
	status_line "RayTracer" "fallo"
	return 1
}

write_readme() {
	cat >"${raw_dir}/README.md" <<'EOF'
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
EOF
}

: >"${log_file}"
{
	echo "# Render Status"
	echo
} >"${status_file}"
echo "Rendering Geant4/Meiga geometry assets" | tee -a "${log_file}"

write_environment
write_readme

rm -f "${raw_dir}"/g4_*.wrl "${raw_dir}/g4_trajectories.wrl" \
	"${figures_dir}/meiga_geometry_view.png" "${figures_dir}/meiga_trajectories_view.png"

if ! command -v cmake >/dev/null 2>&1; then
	warn "cmake is not available; cannot build or run Meiga geometry visualization."
else
	run_logged cmake -S "${repo_root}/src" -B "${repo_root}/build" || warn "cmake configure failed"
	run_logged cmake --build "${repo_root}/build" --target G4HodoscopeSimulator -j2 || warn "G4HodoscopeSimulator build failed"
fi

if [ -x "${app_dir}/G4HodoscopeSimulator" ]; then
	cp "${src_app_dir}/DetectorList.xml" "${raw_dir}/DetectorList.xml"
	write_geometry_config "${raw_dir}/G4HodoscopeSimulator_geometry.json" "VRML2FILE" "true" "false" "1" "./geometry_output.json"
	run_meiga_config "${raw_dir}/G4HodoscopeSimulator_geometry.json" "VRML geometry" || warn "Meiga VRML geometry run failed."
else
	warn "G4HodoscopeSimulator executable not found at ${app_dir}"
fi

if file_nonempty "${raw_dir}/g4_00.wrl"; then
	echo "Generated ${raw_dir}/g4_00.wrl" | tee -a "${log_file}"
	status_line "WRL generado" "sí"
else
	warn "No g4_00.wrl was generated."
	status_line "WRL generado" "no"
fi

geometry_manual_viewpoint_file="${raw_dir}/manual_viewpoint_geometry.wrl.txt"
trajectories_manual_viewpoint_file="${raw_dir}/manual_viewpoint_trajectories.wrl.txt"

geometry_front_wrl="${raw_dir}/g4_00_zoom_front.wrl"
geometry_oblique_wrl="${raw_dir}/g4_00_zoom_oblique.wrl"
geometry_side_wrl="${raw_dir}/g4_00_zoom_side.wrl"
trajectories_front_wrl="${raw_dir}/g4_trajectories_zoom_front.wrl"
trajectories_oblique_wrl="${raw_dir}/g4_trajectories_zoom_oblique.wrl"
trajectories_side_wrl="${raw_dir}/g4_trajectories_zoom_side.wrl"

geometry_front_png="${figures_dir}/meiga_geometry_view_front.png"
geometry_oblique_png="${figures_dir}/meiga_geometry_view_oblique.png"
geometry_side_png="${figures_dir}/meiga_geometry_view_side.png"
trajectories_front_png="${figures_dir}/meiga_trajectories_view_front.png"
trajectories_oblique_png="${figures_dir}/meiga_trajectories_view_oblique.png"
trajectories_side_png="${figures_dir}/meiga_trajectories_view_side.png"

final_geometry_viewpoint="HODOSCOPE_FRONT_CLOSE"
final_trajectories_viewpoint="HODOSCOPE_FRONT_CLOSE"

if [ -s "${geometry_manual_viewpoint_file}" ]; then
	geometry_viewpoint_name="$(grep -Eo 'DEF[[:space:]]+[A-Za-z0-9_]+[[:space:]]+Viewpoint' "${geometry_manual_viewpoint_file}" | head -n 1 | awk '{print $2}')"
	if [ -n "${geometry_viewpoint_name}" ]; then
		final_geometry_viewpoint="${geometry_viewpoint_name}"
	fi
fi

if [ -s "${trajectories_manual_viewpoint_file}" ]; then
	trajectories_viewpoint_name="$(grep -Eo 'DEF[[:space:]]+[A-Za-z0-9_]+[[:space:]]+Viewpoint' "${trajectories_manual_viewpoint_file}" | head -n 1 | awk '{print $2}')"
	if [ -n "${trajectories_viewpoint_name}" ]; then
		final_trajectories_viewpoint="${trajectories_viewpoint_name}"
	fi
fi

render_zoom_variant \
	"${raw_dir}/g4_00.wrl" \
	"${geometry_front_wrl}" \
	"${geometry_front_png}" \
	"${final_geometry_viewpoint}" \
	"0 0 220" \
	"0 0 1 0" \
	"0.18" \
	"${geometry_manual_viewpoint_file}" \
	"geometry / front close"

render_zoom_variant \
	"${raw_dir}/g4_00.wrl" \
	"${geometry_oblique_wrl}" \
	"${geometry_oblique_png}" \
	"HODOSCOPE_OBLIQUE_CLOSE" \
	"120 -160 180" \
	"0.65 0.35 0.67 0.9" \
	"0.22" \
	"" \
	"geometry / oblique close"

render_zoom_variant \
	"${raw_dir}/g4_00.wrl" \
	"${geometry_side_wrl}" \
	"${geometry_side_png}" \
	"HODOSCOPE_SIDE_Z" \
	"220 0 55" \
	"0 1 0 1.5708" \
	"0.20" \
	"" \
	"geometry / side z"

cp -f "${geometry_front_wrl}" "${raw_dir}/g4_00_zoom.wrl" || true
cp -f "${geometry_front_png}" "${figures_dir}/meiga_geometry_view.png" || true

if ! file_nonempty "${figures_dir}/meiga_geometry_view.png"; then
	try_raytracer || true
fi

if [ -x "${app_dir}/G4HodoscopeSimulator" ]; then
	write_geometry_config "${raw_dir}/G4HodoscopeSimulator_trajectories.json" "VRML2FILE" "true" "true" "10" "./trajectory_output.json"
	run_meiga_config "${raw_dir}/G4HodoscopeSimulator_trajectories.json" "VRML trajectories" || warn "Meiga trajectory VRML run failed."
	latest_wrl="$(find "${raw_dir}" -maxdepth 1 -name 'g4_*.wrl' -type f -printf '%T@ %p\n' | sort -n | tail -n 1 | cut -d' ' -f2-)"
	if [ -n "${latest_wrl}" ] && [ "${latest_wrl}" != "${raw_dir}/g4_00.wrl" ]; then
		cp "${latest_wrl}" "${raw_dir}/g4_trajectories.wrl" || true
	fi
fi

if file_nonempty "${raw_dir}/g4_trajectories.wrl"; then
	render_zoom_variant \
		"${raw_dir}/g4_trajectories.wrl" \
		"${trajectories_front_wrl}" \
		"${trajectories_front_png}" \
		"${final_trajectories_viewpoint}" \
		"0 0 220" \
		"0 0 1 0" \
		"0.18" \
		"${trajectories_manual_viewpoint_file}" \
		"trajectories / front close"

	render_zoom_variant \
		"${raw_dir}/g4_trajectories.wrl" \
		"${trajectories_oblique_wrl}" \
		"${trajectories_oblique_png}" \
		"HODOSCOPE_OBLIQUE_CLOSE" \
		"120 -160 180" \
		"0.65 0.35 0.67 0.9" \
		"0.22" \
		"" \
		"trajectories / oblique close"

	render_zoom_variant \
		"${raw_dir}/g4_trajectories.wrl" \
		"${trajectories_side_wrl}" \
		"${trajectories_side_png}" \
		"HODOSCOPE_SIDE_Z" \
		"220 0 55" \
		"0 1 0 1.5708" \
		"0.20" \
		"" \
		"trajectories / side z"

	cp -f "${trajectories_front_wrl}" "${raw_dir}/g4_trajectories_zoom.wrl" || true
	cp -f "${trajectories_front_png}" "${figures_dir}/meiga_trajectories_view.png" || true
else
	echo "No g4_trajectories.wrl found; skipping trajectory render." | tee -a "${log_file}"
fi

if file_nonempty "${figures_dir}/meiga_geometry_view.png"; then
	status_line "PNG generado" "sí"
else
	status_line "PNG generado" "no"
fi

if file_nonempty "${figures_dir}/meiga_trajectories_view.png"; then
	status_line "PNG trayectorias generado" "sí"
else
	status_line "PNG trayectorias generado" "no"
fi

if [ -n "${castle_model_viewer_exe}" ]; then
	status_line "Castle Model Viewer path usado" "${castle_model_viewer_exe}"
else
	status_line "Castle Model Viewer path usado" "no encontrado"
fi

status_line "Viewpoint final geometría" "${final_geometry_viewpoint}"
status_line "Viewpoint final trayectorias" "${final_trajectories_viewpoint}"
status_line "geometry front png" "$(if file_nonempty "${geometry_front_png}"; then echo sí; else echo no; fi)"
status_line "geometry oblique png" "$(if file_nonempty "${geometry_oblique_png}"; then echo sí; else echo no; fi)"
status_line "geometry side png" "$(if file_nonempty "${geometry_side_png}"; then echo sí; else echo no; fi)"
status_line "trajectories front png" "$(if file_nonempty "${trajectories_front_png}"; then echo sí; else echo no; fi)"
status_line "trajectories oblique png" "$(if file_nonempty "${trajectories_oblique_png}"; then echo sí; else echo no; fi)"
status_line "trajectories side png" "$(if file_nonempty "${trajectories_side_png}"; then echo sí; else echo no; fi)"

status_line "castle-model-viewer disponible" "$(if [ -n "${castle_model_viewer_exe}" ]; then echo sí; else echo no; fi)"
status_line "geometría WRL existe" "$(if file_nonempty "${raw_dir}/g4_00.wrl"; then echo sí; else echo no; fi)"
status_line "trayectorias WRL existe" "$(if file_nonempty "${raw_dir}/g4_trajectories.wrl"; then echo sí; else echo no; fi)"
status_line "meiga_geometry_view.png generado" "$(if file_nonempty "${figures_dir}/meiga_geometry_view.png"; then echo sí; else echo no; fi)"
status_line "meiga_trajectories_view.png generado" "$(if file_nonempty "${figures_dir}/meiga_trajectories_view.png"; then echo sí; else echo no; fi)"

if file_nonempty "${figures_dir}/meiga_geometry_view.png"; then
	status_line "método usado para geometría" "Castle Model Viewer con Viewpoint explícito"
else
	status_line "método usado para geometría" "ninguno"
fi

if file_nonempty "${figures_dir}/meiga_trajectories_view.png"; then
	status_line "método usado para trayectorias" "Castle Model Viewer con Viewpoint explícito"
else
	status_line "método usado para trayectorias" "ninguno"
fi

exit 0
