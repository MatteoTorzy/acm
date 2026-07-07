import numpy as np
import platform
import subprocess
from pathlib import Path


def postgmsh(
    mshfile: str,
    analysis,
    nodes,
    V: np.ndarray,
    ndof: int = 2,
    out_msh: str = "output/out.msh",
    opti_geo: str = "output/opti.geo",
    post_msh: str = "output/post.msh",
    launch_gmsh: bool = True,
):

    """
    Assumes:
      - analysis.NN = number of nodes
      - nodes.dof is (NN, ndof) with:
            -1  -> constrained dof
            >=0 -> global unknown number (1-based)
      - nodes.U is (NN, ndof) float array
      - V is (neq, neigs) array of eigenvectors
      - mshfile is the original mesh file to merge in Gmsh
    """

    NN = int(analysis.NN)

    V = np.asarray(V, dtype=float)
    if V.ndim == 1:
        V = V[:, None]

    neigs = int(V.shape[1])

    # ----------------------------
    # Write out.msh (NodeData blocks)
    # ----------------------------
    with open(out_msh, "w", encoding="utf-8") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

        for ieig in range(0, neigs):
            # Fill nodal displacements from eigenvector ieig
            for n in range(0, NN):
                for idir in range(0, ndof):
                    dof = nodes.dof[n, idir]
                    if dof >= 0:
                        nodes.U[n, idir] = V[dof,ieig]
                    else:
                        nodes.U[n,idir] = 0.0

            # displacement output: one vector field per eigenmode
            f.write("$NodeData\n")
            f.write("1\n")
            f.write(f'"U{ieig}"\n')
            f.write("1\n")
            f.write("0.0\n")
            f.write("3\n")
            f.write("0\n")
            f.write("3\n")
            f.write(f"{NN}\n")

            for n in range(0, NN):
                ux = float(nodes.U[n, 0])
                uy = float(nodes.U[n, 1]) 
                f.write(f"{n+1} {ux:.15E} {uy:.15E} {0.0:.15E}\n")

            f.write("$EndNodeData\n")

    # ----------------------------
    # Write opti.geo
    # ----------------------------
    with open(opti_geo, "w", encoding="utf-8") as f:
        f.write("General.BackgroundGradient=0;\n")

        for ieig in range(1, neigs + 1):
            if ieig == 1:
                f.write("View[0].Visible = 1;\n")
            else:
                f.write(f"View[{ieig - 1}].Visible = 0;\n")

            f.write(f"View[{ieig - 1}].VectorType = 5;\n")
            f.write(f"View[{ieig - 1}].NbIso=21;\n")
            f.write(f"View[{ieig - 1}].IntervalsType=3;\n")

    # ----------------------------
    # Write post.msh (merge script)
    # ----------------------------
    mshfile = "../" + mshfile
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
        subprocess.run(["gmsh.exe", post_path], check=False)

    elif "darwin" in sysname:
        subprocess.run(
            ["/Applications/gmsh.app/Contents/MacOS/gmsh", post_path],
            check=False,
        )

    elif "linux" in sysname:
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
