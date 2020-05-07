import trimesh
import trimesh.proximity as proximity
from rtree import index
import numpy as np

try:
    import steps
    import steps.utilities.meshio as meshio
    import steps.geom as sgeom
except ImportError:
    raise ImportError("Unable to import STEPS. Please install it from http://steps.sourceforge.net/.")

def gen_mesh_tet_spatial_index(tetmesh, import_scale):
    """
    Generate a spatial index for all tetrahedrons in a STEPS Tetmesh object.

    Parameters:
        * tetmesh: A steps.geom.Tetmesh object
        * import_scale: The scale used when import the mesh data using steps.utilities.meshio

    Note: the import_scale is the scale value used when the mesh is imported
    to STEPS using steps.utilities.meshio. The value is usually 1e-6 (micrometer)
    or 1e-9 (nanometer).
    """
    p = index.Property()
    p.dimension = 3
    spatial_index = index.Index(properties=p)

    for tet in range(tetmesh.ntets):
        verts = tetmesh.getTet(tet)
        vert_cords = [tetmesh.getVertex(v) for v in verts]
        vert_max = np.maximum(vert_cords[0], vert_cords[1])
        vert_max = np.maximum(vert_cords[2], vert_max)
        vert_min = np.minimum(vert_cords[0], vert_cords[1])
        vert_min = np.minimum(vert_cords[2], vert_min)
        bounding_box = np.append(vert_min, vert_max) / import_scale
        spatial_index.insert(tet, bounding_box)
    return spatial_index

def add_tet_ROIs(tetmesh, import_scale, spatial_index, boundary_mesh_files, roi_labels):
    """
    Identify the tetrahedral elements in each ROI label, using a list of boundary mesh files.
    Create each ROI in the Tetmesh object.
    
    Parameters:
        * tetmesh: A steps.geom.Tetmesh object
        * import_scale: The scale used when import the mesh data using steps.utilities.meshio
        * spatial_index: A spatial index created from gen_mesh_tet_spatial_index
        * boundary_mesh_files: A list of file locations to surface meshes whose combinations define 
        the ROI boundaries.
        * roi_labels: A dict of [ROI_ID : ROI_SIGNATURES] pairs.
        ROI_ID is the name id of the ROI will be stored in tetmesh.
        ROI_SIGNITURES is a sign string consists of only "+", "-" or "*", with length of len(boundary_mesh_files)
        For a boundary mesh boundary_mesh_files[b], ROI_SIGNATURES[b] should be: 
            "+", if ROI elements are not inside the boundary
            "-", if ROI elements are inside the boundary
            "*", if ROI elements have no particular spatial relationship with the boundary 
    """
    boundary_meshes = [trimesh.load(f) for f in boundary_mesh_files]
    all_centers = [np.array(tetmesh.getTetBarycenter(t)) / import_scale for t in range(tetmesh.ntets)]
    tets_in_boundary_mesh = []
    for b_mesh in boundary_meshes:
        bounding_box = np.append(b_mesh.bounds[0], b_mesh.bounds[1])
        check_tets = list(set(spatial_index.intersection(bounding_box)))
        check_tet_centers = [all_centers[tet] for tet in check_tets]
        signed_distances = proximity.signed_distance(b_mesh, check_tet_centers)
        positive_idxs = np.where(signed_distances >= 0.0)[0]
        owned_tets = [check_tets[t] for t in positive_idxs]
        tets_in_boundary_mesh.append(owned_tets)

    for seg_label_key in roi_labels:
        label = roi_labels[seg_label_key]
        assert(len(label) == len(boundary_meshes))
        seg_tets = set(range(tetmesh.ntets))
        for b in range(len(label)):
            sign = label[b]
            if sign == '-':
                seg_tets = seg_tets.intersection(tets_in_boundary_mesh[b])
            elif sign == "+":
                seg_tets = seg_tets - set(tets_in_boundary_mesh[b])
            else:
                continue
        tetmesh.addROI(seg_label_key, sgeom.ELEM_TET, seg_tets)