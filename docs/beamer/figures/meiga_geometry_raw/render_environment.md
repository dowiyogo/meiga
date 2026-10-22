# Render Environment

```text
$ /home/reriosto/castle-model-viewer/castle-model-viewer --help || true
Viewer for all 3D and 2D model formats supported by Castle Game Engine:
glTF, X3D, VRML, Spine, sprite sheets and more.
You can navigate in the (animated and interactive) scene,
with collision-checking, gravity, and a wealth of graphic effects.
You can also convert models in other formats to X3D.

Call as
  castle-model-viewer [OPTIONS]... [FILE-NAME-TO-OPEN]

You can provide FILE-NAME-TO-OPEN on the command-line.

Available options:
  -h / --help           Print this help message and exit.
  -v / --version        Print the version number and exit.
  -H / --hide-extras    Do not show anything extra (like status text or
                        toolbar or bounding box) when program starts. Show
                        only the loaded model.
  --hide-menu           Hide menu bar.
  --no-x3d-extensions   Do not use Castle Game Engine extensions to X3D.
                        Particularly useful when combined with --write, to
                        have X3D valid in all browsers (but less
                        functional).
  --screenshot TIME IMAGE-FILE-NAME
                        Take a screenshot of the loaded scene at given
                        TIME, and save it to IMAGE-FILE-NAME. You most
                        definitely want to pass 3D model file to load at
                        command-line too, otherwise we'll just make a
                        screenshot of the default black scene.
  --screenshot-range TIME-BEGIN TIME-STEP FRAMES-COUNT FILE-NAME
                        Take a FRAMES-COUNT number of screenshots from
                        TIME-BEGIN by step TIME-STEP. Save them to a single
                        movie file (like .avi) (ffmpeg must be installed
                        and available on $PATH for this) or to a sequence
                        of image files (FILE-NAME must then be specified
                        like image@counter(4).png).
  --screenshot-transparent
                        Screenshots background is transparent. Useful only
                        together with --screenshot-range or --screenshot
                        options.
  --viewpoint NAME      Use the viewpoint with given name or index as
                        initial. Especially useful to make a screenshot
                        from given viewpoint.
  --anti-alias AMOUNT   Use full-screen anti-aliasing. Argument AMOUNT is
                        an integer >= 0. Exact 0 means "no anti-aliasing",
                        this is the default. Each successive integer
                        generally makes method one step better. Especially
                        useful to make a screenshot with anti-aliasing
                        quality.
  --project DIR         Point to Castle Game Engine project directory (or
                        CastleEngineManifest.xml file) to resolve the
                        "castle-data:/" URLs in files.

Sound options:
  --audio-device DEVICE-NAME
                        Choose sound output device.
                        Available devices (2):
                          Default Device (internal name: "")
                          OpenAL Soft
  --no-sound            Turn off sound.

Window options:
  --fullscreen          Set window to full-screen (cover whole screen).
  --window              Set window to not be full-screen.
  --geometry WIDTHxHEIGHT<sign>XOFF<sign>YOFF
                        Set window to not be full-screen, and set initial
                        size and/or position.
  --fullscreen-custom WIDTHxHEIGHT
                        Change desktop resolution and set window to
                        full-screen.
  --pretend-touch-device
                        Pretend this is a device with a touch screen, for
                        debugging purposes.

Debug options:
  --debug-log-cache     Write log info, including cache.
  --debug-log-shaders   Write log info, including shader source and log.
  --debug-log-changes   Write log info, including VRML/X3D graph changes.
  --debug-log-videos    Write log info, including videos loading and cache.
  --debug-texture-memory
                        Profile GPU texture memory usage.
  --capabilities automatic|force-fixed-function|force-modern
                        Force OpenGL context to have specific capabilities,
                        to test rendering on modern or ancient GPUs.

Deprecated options:
  --write               DEPRECATED. Load the model and save it as X3D or
                        VRML to the standard output, exit. Use
                        --write-encoding to choose encoding. Use instead
                        castle-model-converter for command-line
                        conversions, it has more functionality with simpler
                        options.
  --write-encoding classic|xml
                        DEPRECATED. Choose X3D encoding to use with --write
                        option. Default is "classic". Use instead
                        castle-model-converter for command-line
                        conversions, it has more functionality with simpler
                        options.
  --write-force-x3d     DEPRECATED. Ignored now.

castle-model-viewer version 5.2.0.
Using Castle Game Engine ( https://castle-engine.io/ ) version 7.0-alpha.3.snapshot.
Compiled with FPC 3.2.2.
Platform: Desktop, OS: Linux, CPU: x86_64 (code compiled for 64-bit architecture).

$ /home/reriosto/castle-model-viewer/castle-model-converter --help || true
castle-model-converter: converter for all 3D and 2D model formats
supported by Castle Game Engine (like glTF, X3D, Collada, MD3...).

Call with these parameters on the command-line:

1. The input filename. It can be any filename, URL or - (stdin).
   Required.

2. The output filename. It can be any filename, URL or - (stdout).
   This parameter is optional. If not provided, the output is stdout,
   unless --validate is also used (then the output is not written anywhere).

Available options are:
  -h / --help           Print this help message and exit.
  -v / --version        Print the version number and exit.
  --no-x3d-extensions   Do not use Castle Game Engine extensions. This will
                        output file valid in all X3D browsers (but maybe
                        with some CGE-specific features missing).
  --enable-downloads    Enable (blocking) downloads from the net, e.g. to
                        download a texture or Inlined model referenced by
                        htt(s) protocol).
  --validate            Only validate the input, without any output (the
                        output URL is not allowed in this case). Moreover,
                        if there will be any warning or error, we will exit
                        with non-zero status (by default, only errors cause
                        non-zero status).
  --stdin-url           If input URL is "-", then we read file contents
                        from the standard input. In this case, you can use
                        this option to provide a "pretend" URL for the
                        input. We will use it to resolve relative URLs
                        inside the input (e.g. to glTF binary blobs) and to
                        guess the input file type. Default is "stdin.x3dv"
                        in current directory, so we assume it is X3D
                        (classic encoded), and resolve with respect to the
                        current directory.
  --stdout-url          If output URL is "-", then we write file contents
                        to the standard output. In this case, you can use
                        this option to provide a "pretend" URL for the
                        output. We will use it to determine the output file
                        type, e.g. "out.x3d" to output X3D XML encoding or
                        "out.x3dv" to output X3D classic encoding.
  --float-precision DIGITS
                        Number of digits after the decimal point when
                        writing floating-point numbers. Default is to write
                        all possibly relevant digits. Specify any value >=
                        0 to use this number of digits.
  --encoding classic|xml
                        DEPRECATED. Choose X3D encoding. Do not use -- the
                        2nd parameter should determine the output type,
                        ".x3d" extension says to make X3D XML, ".x3dv" says
                        to make X3D classic. Or use --stdout-url to provide
                        fake URL in case output is to stdout.

castle-model-converter version 5.2.0.
Using Castle Game Engine ( https://castle-engine.io/ ) version 7.0-alpha.3.snapshot.
Compiled with FPC 3.2.2.
Platform: Desktop, OS: Linux, CPU: x86_64 (code compiled for 64-bit architecture).

$ echo ${CASTLE_MODEL_VIEWER-}


$ echo ${CASTLE_MODEL_CONVERTER-}


$ CASTLE_MODEL_VIEWER_DEFAULT
/home/reriosto/castle-model-viewer/castle-model-viewer

$ CASTLE_MODEL_CONVERTER_DEFAULT
/home/reriosto/castle-model-viewer/castle-model-converter

$ CASTLE_MODEL_VIEWER_EXE
/home/reriosto/castle-model-viewer/castle-model-viewer

$ CASTLE_MODEL_CONVERTER_EXE
/home/reriosto/castle-model-viewer/castle-model-converter

$ uname -a
Linux MSI 6.6.114.1-microsoft-standard-WSL2 #1 SMP PREEMPT_DYNAMIC Mon Dec  1 20:46:23 UTC 2025 x86_64 x86_64 x86_64 GNU/Linux

$ cat /etc/os-release || true
NAME="AlmaLinux"
VERSION="9.7 (Moss Jungle Cat)"
ID="almalinux"
ID_LIKE="rhel centos fedora"
VERSION_ID="9.7"
PLATFORM_ID="platform:el9"
PRETTY_NAME="AlmaLinux 9.7 (Moss Jungle Cat)"
ANSI_COLOR="0;34"
LOGO="fedora-logo-icon"
CPE_NAME="cpe:/o:almalinux:almalinux:9::baseos"
HOME_URL="https://almalinux.org/"
DOCUMENTATION_URL="https://wiki.almalinux.org/"
BUG_REPORT_URL="https://bugs.almalinux.org/"

ALMALINUX_MANTISBT_PROJECT="AlmaLinux-9"
ALMALINUX_MANTISBT_PROJECT_VERSION="9.7"
REDHAT_SUPPORT_PRODUCT="AlmaLinux"
REDHAT_SUPPORT_PRODUCT_VERSION="9.7"
SUPPORT_END=2032-06-01

$ which blender || true
/usr/bin/blender

$ which wslpath || true
/usr/bin/wslpath

$ ls /mnt/c/Program\ Files/Blender\ Foundation || true

$ echo ${WINDOWS_BLENDER_EXE-}


$ WINDOWS_BLENDER_EXE=... --version
Windows Blender not detected

$ blender --version || true
Color management: using fallback mode for management
Color management: Error could not find role data role.
Blender 3.3.21
	build date: 2026-01-14
	build time: 00:00:00
	build commit date: 1970-01-01
	build commit time: 00:00
	build hash: unknown
	build platform: Linux
	build type: Release
	build c flags:  -Wall -Werror=implicit-function-declaration -Werror=return-type -Werror=vla -Wstrict-prototypes -Wmissing-prototypes -Wno-char-subscripts -Wno-unknown-pragmas -Wpointer-arith -Wunused-parameter -Wwrite-strings -Wlogical-op -Wundef -Winit-self -Wmissing-include-dirs -Wno-div-by-zero -Wtype-limits -Wformat-signedness -Wrestrict -Wnonnull -Wabsolute-value -Wuninitialized -Wredundant-decls -Wshadow -Wno-error=unused-but-set-variable -Wimplicit-fallthrough=5 -O2 -flto=auto -ffat-lto-objects -fexceptions -g -grecord-gcc-switches -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -Wp,-D_GLIBCXX_ASSERTIONS -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1 -fstack-protector-strong -specs=/usr/lib/rpm/redhat/redhat-annobin-cc1  -m64 -march=x86-64-v2 -mtune=generic -fasynchronous-unwind-tables -fstack-clash-protection -fcf-protection -Wl,--as-needed -fopenmp -std=gnu11 -msse -pipe -fPIC -funsigned-char -fno-strict-aliasing -ffp-contract=off -msse2  
	build c++ flags:  -Wuninitialized -Wredundant-decls -Wall -Wno-invalid-offsetof -Wno-sign-compare -Wlogical-op -Winit-self -Wmissing-include-dirs -Wno-div-by-zero -Wtype-limits -Werror=return-type -Wno-char-subscripts -Wno-unknown-pragmas -Wpointer-arith -Wunused-parameter -Wwrite-strings -Wundef -Wformat-signedness -Wrestrict -Wno-suggest-override -Wuninitialized -Wundef -Wmissing-declarations -Wimplicit-fallthrough=5 -O2 -flto=auto -ffat-lto-objects -fexceptions -g -grecord-gcc-switches -pipe -Wall -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -Wp,-D_GLIBCXX_ASSERTIONS -specs=/usr/lib/rpm/redhat/redhat-hardened-cc1 -fstack-protector-strong -specs=/usr/lib/rpm/redhat/redhat-annobin-cc1  -m64 -march=x86-64-v2 -mtune=generic -fasynchronous-unwind-tables -fstack-clash-protection -fcf-protection -Wl,--as-needed -fopenmp -msse -pipe -fPIC -funsigned-char -fno-strict-aliasing -ffp-contract=off -msse2  
	build link flags:  -Wl,--version-script='/builddir/build/BUILD/blender-3.3.21/source/creator/symbols_unix.map'
	build system: CMake

$ which xvfb-run || true
/usr/bin/xvfb-run

$ which glxinfo || true
/usr/bin/glxinfo

$ which glxgears || true
/usr/bin/glxgears

$ which python3 || true
/usr/bin/python3

$ which magick || true
/usr/bin/which: no magick in (/opt/hep/geant4/bin:/opt/hep/root/bin:/home/reriosto/.cache/cloud-code/m2c/bin:/home/reriosto/.vscode-server/data/User/globalStorage/github.copilot-chat/debugCommand:/home/reriosto/.vscode-server/data/User/globalStorage/github.copilot-chat/copilotCli:/home/reriosto/.vscode-server/bin/0958016b2af9f09bb4257e0df4a95e2f90590f9f/bin/remote-cli:/opt/hep/geant4/bin:/home/reriosto/.local/bin:/home/reriosto/bin:/usr/share/Modules/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/usr/lib/wsl/lib:/mnt/c/Users/renem/AppData/Local/Programs/Microsoft VS Code:/mnt/c/Program Files/ImageMagick-7.1.2-Q16-HDRI:/mnt/c/Program Files (x86)/Common Files/Oracle/Java/java8path:/mnt/c/Program Files (x86)/Common Files/Oracle/Java/javapath:/mnt/c/WINDOWS/system32:/mnt/c/WINDOWS:/mnt/c/WINDOWS/System32/Wbem:/mnt/c/WINDOWS/System32/WindowsPowerShell/v1.0/:/mnt/c/WINDOWS/System32/OpenSSH/:/mnt/c/Program Files/Wolfram Research/WolframScript/:/mnt/c/root_v6.24.02/bin:/mnt/c/Program Files/CMake/bin:/mnt/c/Program Files (x86)/QuickTime/QTSystem/:/mnt/c/Program Files (x86)/HIK-share/:/mnt/c/Program Files/Calibre2/:/mnt/d/ffmpeg/bin:/mnt/d/realesrgan:/mnt/c/Program Files/Git/cmd:/mnt/c/Program Files/PowerShell/7/:/mnt/c/Program Files/gs/gs10.07.0/bin:/mnt/c/Program Files (x86)/HP/HP OCR/DB_Lib/:/mnt/c/Program Files/HP/Common/HPDestPlgIn/:/mnt/c/Program Files (x86)/HP/Common/HPDestPlgIn/:/mnt/c/Users/renem/AppData/Local/Programs/Python/Python312/Scripts/:/mnt/c/Users/renem/AppData/Local/Programs/Python/Python312/:/mnt/c/Users/renem/AppData/Local/Programs/Python/Python39/Scripts/:/mnt/c/Users/renem/AppData/Local/Programs/Python/Python39/:/mnt/c/Users/renem/AppData/Local/Programs/Python/Python39-32/Scripts/:/mnt/c/Users/renem/AppData/Local/Programs/Python/Python39-32/:/mnt/c/Users/renem/AppData/Local/Microsoft/WindowsApps:/mnt/c/Users/renem/AppData/Local/Programs/Microsoft VS Code/bin:/mnt/c/Users/renem/AppData/Local/Programs/MiKTeX/miktex/bin/x64/:/mnt/c/Program Files (x86)/Skyworks/ClockBuilder Pro/Bin:/mnt/c/Users/renem/AppData/Local/Microsoft/WinGet/Links:/mnt/c/Program Files/HP/Common/HPDestPlgIn/:/mnt/c/Program Files (x86)/HP/Common/HPDestPlgIn/)

$ which convert || true
/usr/bin/convert

$ rpm -qa | grep -E "mesa|glvnd|blender|xorg-x11-server-Xvfb|ImageMagick" | sort || true
blender-3.3.21-3.el9.x86_64
ImageMagick-6.9.13.25-1.el9.x86_64
ImageMagick-libs-6.9.13.25-1.el9.x86_64
libglvnd-1.3.4-1.el9.x86_64
libglvnd-core-devel-1.3.4-1.el9.x86_64
libglvnd-devel-1.3.4-1.el9.x86_64
libglvnd-egl-1.3.4-1.el9.x86_64
libglvnd-gles-1.3.4-1.el9.x86_64
libglvnd-glx-1.3.4-1.el9.x86_64
libglvnd-opengl-1.3.4-1.el9.x86_64
mesa-demos-8.4.0-12.20210504git0f9e7d9.el9.x86_64
mesa-dri-drivers-25.0.7-5.el9_7.alma.1.x86_64
mesa-filesystem-25.0.7-5.el9_7.alma.1.x86_64
mesa-libEGL-25.0.7-5.el9_7.alma.1.x86_64
mesa-libEGL-devel-25.0.7-5.el9_7.alma.1.x86_64
mesa-libgbm-25.0.7-5.el9_7.alma.1.x86_64
mesa-libGL-25.0.7-5.el9_7.alma.1.x86_64
mesa-libGL-devel-25.0.7-5.el9_7.alma.1.x86_64
mesa-libGLU-9.0.1-6.el9.x86_64
mesa-libGLU-devel-9.0.1-6.el9.x86_64
mesa-vulkan-drivers-25.0.7-5.el9_7.alma.1.x86_64
xorg-x11-server-Xvfb-1.20.11-33.el9_7.x86_64

$ echo $DISPLAY
:0

$ echo $LIBGL_ALWAYS_SOFTWARE
0

$ echo $MESA_LOADER_DRIVER_OVERRIDE


$ LIBGL_ALWAYS_SOFTWARE=1 glxinfo -B || true
Dropped Escape call with ulEscapeCode : 0x03007703
name of display: :0
display: :0  screen: 0
direct rendering: Yes
Extended renderer info (GLX_MESA_query_renderer):
    Vendor: Microsoft Corporation (0xffffffff)
    Device: D3D12 (Radeon RX 5500M) (0xffffffff)
    Version: 24.0.5
    Accelerated: yes
    Video memory: 11940MB
    Unified memory: no
    Preferred profile: core (0x1)
    Max core profile version: 4.5
    Max compat profile version: 4.5
    Max GLES1 profile version: 1.1
    Max GLES[23] profile version: 3.1
OpenGL vendor string: Microsoft Corporation
OpenGL renderer string: D3D12 (Radeon RX 5500M)
OpenGL core profile version string: 4.5 (Core Profile) Mesa 24.0.5
OpenGL core profile shading language version string: 4.60
OpenGL core profile context flags: (none)
OpenGL core profile profile mask: core profile

OpenGL version string: 4.5 (Compatibility Profile) Mesa 24.0.5
OpenGL shading language version string: 4.60
OpenGL context flags: (none)
OpenGL profile mask: compatibility profile

OpenGL ES profile version string: OpenGL ES 3.1 Mesa 24.0.5
OpenGL ES profile shading language version string: OpenGL ES GLSL ES 3.10

```
