import platform
import subprocess
from pathlib import Path


_fmidhist = None
_output_counter = 0


def postgmsh(
    mshfile: str,
    analysis,
    nodes,
    post_flag: int,
    Dt: float = 0.0,
    step: int = 0,
    out_msh: str = "input/out.msh",
    opti_geo: str = "input/opti.geo",
    post_msh: str = "input/post.msh",
    launch_gmsh: bool = True,
):
    """
    Post-processing for transient structural analysis.

    Parameters
    ----------
    mshfile : str
        Original mesh file to merge/open in Gmsh.
    analysis : object
        Must provide analysis.NN.
    nodes : object
        Must provide nodes.U with shape (NN, 2) for planar displacement.
    post_flag : int
        0 -> initialize history file
        1 -> append one time frame
        2 -> finalize files and optionally launch Gmsh
    Dt : float
        Time step size.
    step : int
        Current time-step index.
    out_msh, opti_geo, post_msh : str
        Output filenames.
    launch_gmsh : bool
        Whether to launch Gmsh at finalization.
    """
    global _fmidhist, _output_counter

    NN = int(analysis.NN)

    out_path = Path(out_msh)
    opti_path = Path(opti_geo)
    post_path = Path(post_msh)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    opti_path.parent.mkdir(parents=True, exist_ok=True)
    post_path.parent.mkdir(parents=True, exist_ok=True)

    # ----------------------------
    # Initialize history file
    # ----------------------------
    if post_flag == 0:
        _output_counter = 0
        if _fmidhist is not None and not _fmidhist.closed:
            _fmidhist.close()
        _fmidhist = open(out_path, "w", encoding="utf-8")
        _fmidhist.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        return

    # ----------------------------
    # Append standard output frame
    # ----------------------------
    elif post_flag == 1:
        if _fmidhist is None or _fmidhist.closed:
            raise RuntimeError("History file is not open. Call postgmsh(..., post_flag=0) first.")

        _output_counter += 1
        _fmidhist.write("$NodeData\n")
        _fmidhist.write("1\n")
        _fmidhist.write('"U"\n')
        _fmidhist.write("1\n")
        _fmidhist.write(f"{Dt * step}\n")
        _fmidhist.write("3\n")
        _fmidhist.write(f"{step}\n")
        _fmidhist.write("3\n")
        _fmidhist.write(f"{NN}\n")

        for n in range(1, NN + 1):
            ux = float(nodes.U[n - 1, 0])
            uy = float(nodes.U[n - 1, 1])
            _fmidhist.write(f"{n} {ux:.15E} {uy:.15E} 0\n")

        _fmidhist.write("$EndNodeData\n")
        _fmidhist.flush()
        return

    # ----------------------------
    # Finalize files and launch Gmsh
    # ----------------------------
    elif post_flag == 2:
        if _fmidhist is not None and not _fmidhist.closed:
            _fmidhist.close()

        with open(opti_path, "w", encoding="utf-8") as f:
            f.write("General.BackgroundGradient=0;\n")
            f.write("Mesh.SurfaceEdges=0;\n")
            f.write("Mesh.SurfaceFaces=0;\n")
            f.write("View[0].Visible = 1;\n")
            f.write("View[0].VectorType = 5;\n")
            f.write("View[0].NbIso=21;\n")
            f.write("View[0].IntervalsType=3;\n")

        mshfile_rel = "../" + mshfile
        with open(post_path, "w", encoding="utf-8") as f:
            f.write(f'Merge "{mshfile_rel}";\n')
            f.write(f'Merge "{out_path.name}";\n')
            f.write(f'Merge "{opti_path.name}";\n')

        if not launch_gmsh:
            return

        sysname = platform.system().lower()
        post_path_str = str(post_path)

        if "windows" in sysname:
            subprocess.run(["gmsh.exe", post_path_str], check=False)

        elif "darwin" in sysname:
            subprocess.run(
                ["/Applications/gmsh.app/Contents/MacOS/gmsh", post_path_str],
                check=False,
            )

        elif "linux" in sysname:
            try:
                subprocess.run(["gmsh", post_path_str], check=False)
            except FileNotFoundError:
                gmsh_bins = sorted(Path(".").glob("gmsh-*/bin/gmsh"))
                if gmsh_bins:
                    subprocess.run([str(gmsh_bins[0]), post_path_str], check=False)
                else:
                    print("Gmsh not found (PATH or ./gmsh-*/bin/gmsh).")

        else:
            print("Platform not supported")

        return

    else:
        raise ValueError("post_flag must be 0, 1, or 2")
