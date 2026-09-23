import numpy as np
import platform
import subprocess
from pathlib import Path


def postgmsh(
    mshfile: str,
    analysis,
    nodes,
    Sn: np.ndarray,          # shape (NN, Sdim)
    ndof: int = 2,
    Sdim: int = 3,
    out_msh: str = "output/out.msh",
    opti_geo: str = "output/opti.geo",
    post_msh: str = "output/post.msh",
    launch_gmsh: bool = True,
):

    """
    Assumes:
      - analysis.NN = number of nodes
      - nodes.U is (NN, ndof) float array 
      - Sn is (NN, Sdim) float array
      - mshfile is the original mesh file to merge in Gmsh
    """

    NN = int(analysis.NN)

    # ----------------------------
    # Write out.msh (NodeData blocks)
    # ----------------------------
    with open(out_msh, "w", encoding="utf-8") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

        # displacement output
        for icomp in range(1, ndof + 1):
            # Same header lines as MATLAB fprintf with %s\n and pieces
            f.write("$NodeData\n")
            f.write("1\n")
            f.write(f"\"U{icomp}\"\n")
            f.write("1\n")
            f.write("0.0\n")
            f.write("3\n")
            f.write("0\n")
            f.write("1\n")
            f.write(f"{NN}\n")

            # MATLAB nodes are 1-based; your nodes.U array is 0-based
            for n in range(1, NN + 1):
                val = float(nodes.U[n - 1, icomp - 1])
                f.write(f"{n} {val:.15E}\n")

            f.write("$EndNodeData\n")

        # stress output
        for icomp in range(1, Sdim + 1):
            f.write("$NodeData\n")
            f.write("1\n")
            f.write(f"\"S{icomp}\"\n")
            f.write("1\n")
            f.write("0.0\n")
            f.write("3\n")
            f.write("0\n")
            f.write("1\n")
            f.write(f"{NN}\n")

            for n in range(1, NN + 1):
                val = float(Sn[n - 1, icomp - 1])
                f.write(f"{n} {val:.15E}\n")

            f.write("$EndNodeData\n")

    # ----------------------------
    # Write opti.geo
    # ----------------------------
    with open(opti_geo, "w", encoding="utf-8") as f:
        f.write("General.BackgroundGradient=0;\n")
        f.write("Mesh.SurfaceEdges=0;\n")
        f.write("Mesh.SurfaceFaces=0;\n")
        f.write("View[0].Visible = 1;\n")
        f.write("View[0].VectorType = 5;\n")
#        f.write("View[0].ShowElement = 1;\n")
        f.write("View[0].NbIso=21;\n")
        f.write("View[0].IntervalsType=3;\n")

        # MATLAB: for n=1:ndof+Sdim-1  -> Python: 1..(ndof+Sdim-1) inclusive
        for n in range(1, ndof + Sdim):
            f.write(f"View[{n}].Visible=0;\n")
            f.write(f"View[{n}].VectorType = 5;\n")
#            f.write(f"View[{n}].ShowElement = 1;\n")
            f.write(f"View[{n}].NbIso=21;\n")
            f.write(f"View[{n}].IntervalsType=3;\n")

    # ----------------------------
    # Write post.msh (merge script)
    # ----------------------------
    mshfile="../" + mshfile
    with open(post_msh, "w", encoding="utf-8") as f:
        f.write(f'Merge "{mshfile}";\n')
        f.write('Merge "out.msh";\n')
        f.write('Merge "opti.geo";\n')

    # ----------------------------
    # Launch Gmsh (optional)
    # ----------------------------
    if not launch_gmsh:
        return

    sysname = platform.system().lower()
    post_path = str(Path(post_msh))

    if "windows" in sysname:
        # assumes gmsh.exe is in PATH or current folder
        subprocess.run(["gmsh.exe", post_path], check=False)

    elif "darwin" in sysname:  # macOS
        subprocess.run(
            ["/Applications/gmsh.app/Contents/MacOS/gmsh", post_path],
            check=False,
        )

    elif "linux" in sysname:
        # Similar spirit to MATLAB: $(pwd)/gmsh-*/bin/gmsh post.msh
        # We try: gmsh in PATH first; if not found, try globbed local gmsh-*/bin/gmsh
        try:
            subprocess.run(["gmsh", post_path], check=False)
        except FileNotFoundError:
            gmsh_bins = sorted(Path(".").glob("gmsh-*/bin/gmsh"))
            if gmsh_bins:
                subprocess.run([str(gmsh_bins[0]), post_path], check=False)
            else:
                print("Gmsh not found (PATH or ./gmsh-*/bin/gmsh).")

    else:
        print("Platform not supported")