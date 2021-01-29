import gmsh
import numpy as np
from rtree import index
import trimesh
import trimesh.proximity as proximity
from typing import NamedTuple
import sys

class _gmsh_ele_t(NamedTuple):
    etype: int       # gmsh type tag
    nnodes: int      # number of nodes
    dim: int         # dimensionality

# tetrahedron
tet_t = _gmsh_ele_t(4, 4, 3)

from contextlib import contextmanager

@contextmanager
def gmsh_fixture():
  try:
    yield gmsh.initialize()
  finally:
    gmsh.finalize()


def gen_tet_spatial_index(mesh, import_scale):
    """
    Generate a spatial index for all tetrahedrons in a gmsh model assuming only tets.

    Parameters:
        * mesh: a gmsh.model.mesh object
        * import_scale: The scale used when importing the mesh data

    """
    ele_ids, node_ids = mesh.getElementsByType(tet_t.etype)

    p = index.Property()
    p.dimension = 3
    spatial_index = index.Index(properties=p)

    # loop over tets
    for tet_id in range(len(ele_ids)):

        # get bounding box
        inf = 1.0e15
        vert_max = np.array([-inf, -inf, -inf], dtype=float)
        vert_min = np.array([inf, inf, inf], dtype=float)

        # loop over the node_ids
        for n in range(tet_t.nnodes):
            node = node_ids[tet_id*tet_t.nnodes+n]
            # Get coordinates (and parametric coor if available) from node id
            node_coo = gmsh.model.mesh.getNode(node)[0]
            vert_max = np.maximum(node_coo, vert_max)
            vert_min = np.minimum(node_coo, vert_min)

        bounding_box = np.append(vert_min, vert_max) / import_scale
        spatial_index.insert(tet_id, bounding_box)

    return spatial_index

def add_tet_ROIs(mesh, import_scale, spatial_index, boundary_mesh_files, roi_labels):
    """
    Identify the tetrahedral elements in each ROI label, using a list of boundary mesh files.
    Create each ROI in the gmsh.model.mesh object.

    Parameters:
        * mesh: a gmsh.model.mesh object
        * import_scale: The scale used when importing the mesh data using steps.utilities.meshio
        * spatial_index: A spatial index created from gen_mesh_tet_spatial_index
        * boundary_mesh_files: A list of file locations to surface meshes whose combinations define
        the ROI boundaries.
        * roi_labels: A dict of [ROI_ID : ROI_SIGNATURES] pairs.
        ROI_ID is the name id of the ROI that will be stored in tetmesh.
        ROI_SIGNATURES is a sign string consisting of only "+", "-" or "*", with length of len(boundary_mesh_files)
        For a boundary mesh boundary_mesh_files[b], ROI_SIGNATURES[b] should be:
            "+", if ROI elements are not inside the boundary
            "-", if ROI elements are inside the boundary
            "*", if ROI elements have no particular spatial relationship with the boundary
    """

    # Load and check stl files
    print("Loading stl files...")
    boundary_meshes = []
    for f in boundary_mesh_files:
        stl_mesh = trimesh.load(f)
        assert stl_mesh.is_watertight, f+" is not watertight"
        print(f)
        boundary_meshes.append(stl_mesh)

    # getBarycenters(elementType, tag, fast, primary)
    # tag = -1 --> all elements
    # fast = 0 --> normalize by number of coordinates
    # primary = 1 --> return only verteces
    all_centers = mesh.getBarycenters(tet_t.etype, -1, 0, 1) / import_scale
    all_centers = np.reshape(all_centers, (-1, 3))
    print("Done loading stl files.")

    print("Tagging ROIs...")
    tets_in_boundary_mesh = []
    for b_mesh in boundary_meshes:
        bounding_box = np.append(b_mesh.bounds[0], b_mesh.bounds[1])
        check_tets = list(set(spatial_index.intersection(bounding_box)))
        check_tet_centers = [all_centers[tet] for tet in check_tets]
        signed_distances = proximity.signed_distance(b_mesh, check_tet_centers)
        positive_idxs = np.where(signed_distances >= 0.0)[0]
        owned_tets = [check_tets[t] for t in positive_idxs]
        tets_in_boundary_mesh.append(owned_tets)

    ROIs = []
    for seg_label_key in roi_labels:
        label = roi_labels[seg_label_key]
        assert(len(label) == len(boundary_meshes))
        seg_tets = set(range(all_centers.shape[0]))
        for b in range(len(label)):
            sign = label[b]
            if sign == '-':
                seg_tets = seg_tets.intersection(tets_in_boundary_mesh[b])
            elif sign == "+":
                seg_tets = seg_tets - set(tets_in_boundary_mesh[b])
            else:
                continue
        ROIs.append(list(seg_tets))
    print("Done tagging ROIs.")

    return ROIs

def tag_mesh_entities(input_mesh, boundary_files, roi_labels, scale_ratio=1.0, \
    output_mesh=None, interactive=False):
    """Tag mesh entities
    Args:
        input_mesh: path to a GMSH mesh
        boundary_mesh_files: A list of file locations to surface meshes whose combinations define
        the ROI boundaries.
        roi_labels: A dict of [ROI_ID : ROI_SIGNATURES] pairs.
        scale_ratio: scale conversion factor between stl and msh
        interactive: to visualize the mesh using the graphical UI provide by gmsh API
        output_mesh: in-place edit if None, path to output mesh otherwise
    """

    gmsh.option.setNumber("General.Terminal", 1)

    # Open file
    gmsh.open(input_mesh)
    if gmsh.model.getDimension() < 0:
        print("Could not find " + input_mesh, file=sys.stderr)
        exit()

    # Print the model name and dimension:
    print('Model ' + gmsh.model.getCurrent() + ' (' +
          str(gmsh.model.getDimension()) + 'D)')

    # Create spatial index
    spatial_index = gen_tet_spatial_index(gmsh.model.mesh, scale_ratio)

    # Tag ROIs
    ROIs = add_tet_ROIs(gmsh.model.mesh, scale_ratio, spatial_index, \
        boundary_files, roi_labels)
    print("Found: ", len(ROIs), "ROIs")

    # Get element ids and node ids in mesh from element type
    ele_ids, node_ids = gmsh.model.mesh.getElementsByType(tet_t.etype)
    nTets = len(ele_ids)
    print("Number of tets in original mesh: ", nTets)

    # Select tets based on ROI list
    nROIs = len(ROIs)
    nEntities = nROIs + 1 # ROIs plus untagged/original volume
    tets_select = [[] for _ in range(nEntities)]     # selected elements
    nodes_select = [[] for _ in range(nEntities)]    # correspoding nodes
    nodes_coords = [[] for _ in range(nEntities)]    # nodes coordinates

    for i in range(nTets):
        is_in_ROIs = False
        # Loop over ROIs to select tag
        for id_roi in range(nROIs):
            for t in ROIs[id_roi]:
                if t == i:
                    is_in_ROIs = True
                    id_entity = id_roi + 1
                    tets_select[id_entity].extend([ele_ids[i]])
                    nodes_select[id_entity].extend(node_ids[tet_t.nnodes*i:\
                        tet_t.nnodes*i+tet_t.nnodes])
                    for n in range(tet_t.nnodes):
                        node = node_ids[i*tet_t.nnodes+n]
                        # Get coordinates (and parametric coor if available) from node id
                        node_coo = gmsh.model.mesh.getNode(node)[0]
                        nodes_coords[id_entity].extend(node_coo)

        if not is_in_ROIs:
            # Default tag, not in any ROI
            id_entity = 0
            tets_select[id_entity].extend([ele_ids[i]])
            nodes_select[id_entity].extend(node_ids[tet_t.nnodes*i:\
                tet_t.nnodes*i+tet_t.nnodes])
            for n in range(tet_t.nnodes):
                node = node_ids[i*tet_t.nnodes+n]
                # Get coordinates (and parametric coor if available) from node id
                node_coo = gmsh.model.mesh.getNode(node)[0]
                nodes_coords[id_entity].extend(node_coo)

    # Clear mesh and model
    if (output_mesh == None):
        output_mesh = gmsh.model.getCurrent()+".msh"
    gmsh.clear()

    # Add new discrete entities: non tagged tets + ROIs
    roi_names = list(roi_labels.keys())
    for i in range(nEntities):
        vnew = gmsh.model.addDiscreteEntity(tet_t.dim)
        gmsh.model.mesh.addNodes(tet_t.dim, vnew, nodes_select[i], nodes_coords[i])
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.mesh.addElementsByType(vnew, tet_t.etype, tets_select[i], nodes_select[i])
        volume_tag = i + 1000
        gmsh.model.addPhysicalGroup(dim=tet_t.dim, tags=[vnew], tag=volume_tag)
        if i == 0:
            gmsh.model.setPhysicalName(dim=tet_t.dim, tag=volume_tag, name="Initial domain")
        else:
            gmsh.model.setPhysicalName(dim=tet_t.dim, tag=volume_tag, name=roi_names[i-1])

    # To visualize the model we can run the graphical user interface from gmsh
    if interactive:
        gmsh.fltk.run()

    gmsh.write(output_mesh)