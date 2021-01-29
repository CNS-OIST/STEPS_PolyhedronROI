"""
This example imports a dendritic spine mesh reconstruction to STEPS,
and creates tetrahedral Regions of Interest (ROI)s using polyhedral surface boundaries.
The mesh with newly created ROIs is then exported to STEPS xml format,
and visualized using steps.visual module.

Glossary
* Dendritic spine: https://en.wikipedia.org/wiki/Dendritic_spine
* Cyto: Cytosol https://en.wikipedia.org/wiki/Cytosol
* ER: Endoplasmic reticulum https://en.wikipedia.org/wiki/Endoplasmic_reticulum
* PSD: Postsynaptic Density https://en.wikipedia.org/wiki/Postsynaptic_density
"""

import steps
import steps.utilities.meshio as meshio
import steps.visual

import pyqtgraph as pg

import polyhedronROI.steps as pRs

def main():
    print("Import mesh to STEPS")
    tetmesh = meshio.importAbaqus("meshes/tets.inp", 1e-6)[0]

    print("Create spatial index")
    spatial_index = pRs.gen_tet_spatial_index(tetmesh, 1e-6)

    print("Add ROIs according to boundaries and signatures")
    boundary_files = ["meshes/ER.stl", "meshes/PSD.stl"]
    roi_labels = {"ER": "-*",   # inside ER.stl, doesn't matter for PSD.stl
                  "PSD": "+-"}  # not inside ER.stl, inside PSD.stl

    pRs.add_tet_ROIs(tetmesh, 1e-6, spatial_index, boundary_files, roi_labels)

    print("Save to STEPS xml file")
    meshio.saveMesh("meshes/labeled_mesh", tetmesh)

    print("Load the xml file")
    loaded_mesh = meshio.loadMesh("meshes/labeled_mesh")[0]

    print("Visualize the ROIs")
    ER_tets = loaded_mesh.getROIData("ER")
    PSD_tets = loaded_mesh.getROIData("PSD")
    Cyto_tets = set(range(loaded_mesh.ntets)) - set(ER_tets)

    print(
        "PSD tets / Cyto tets / ER tets / All tets: %i / %i / %i / %i"
        % (len(PSD_tets), len(Cyto_tets), len(ER_tets), loaded_mesh.ntets)
    )

    tet_groups = {}
    for tet in Cyto_tets:
        if tet in PSD_tets:
            tet_groups[tet] = "PSD"
        else:
            tet_groups[tet] = "Cyto"
    for tet in ER_tets:
        tet_groups[tet] = "ER"

    color_map = dict(
        PSD=[0.0, 0.0, 1.0, 0.1],  # blue, alpha 0.1
        Cyto=[1.0, 1.0, 1.0, 0.1], # grey, alpha 0.1
        ER=[1.0, 0.0, 0.0, 0.1],  # red, alpha 0.1
    )

    app = pg.mkQApp()
    w = steps.visual.TetPartitionDisplay(
        loaded_mesh, tet_groups, color_map=color_map
    )
    app.exec_()


if __name__ == "__main__":
    main()
