# Region of Interest (ROI) labeling for STEPS Tetmesh or gmsh msh using Polyhedral surface boundaries

This module reads a list of surface boundary meshes, and a dictionary of ROI labeling signatures,
then create the ROIs according to the boundaries and labeling signatures.

## Prerequisite
* Python3.6 or above
* [STEPS](http://steps.sourceforge.net/)
* numpy
* trimesh
* Rtree

## Download and Installation
```
git clone https://github.com/CNS-OIST/STEPS_PolyhedronROI.git
cd STEPS_PolyhedronROI
python -m pip install .
```

## Usage

### STEPS
in `Python` interface, import the module
```python
import polyhedronROI.steps as pRs
```

The first step is to create a spatial index of the mesh elements
```python
spatial_index = pRs.gen_tet_spatial_index(TETMESH, IMPORT_SCALE)
```
* `TETMESH`: A steps.geom.Tetmesh object
* `IMPORT_SCALE`: The scale used when the mesh is imported to STEPS.

The tetrahedral ROIs can be defined and added to `TETMESH` using
```python
pRs.add_tet_ROIs(TETMESH, IMPORT_SCALE, spatial_index, BOUNDARY_FILES, ROI_LABELS)
```

* `BOUNDARY_FILES`: A list of strings, each string stores the location of a boundary surface mesh file. for example
```python
boundary_files = ["meshes/ER.stl", "meshes/PSD.stl"]
```
* `ROI_LABELS`: A `dict` of [`ROI_ID` : `ROI_SIGNATURES`] pairs.
    * `ROI_ID`: Name id of the ROI
    * `ROI_SIGNATURES`: A sign string consisting of only `+`, `-` and `*`. The length of the
    string should be the same as `len(BOUNDARY_FILES)`. For a boundary file `BOUNDARY_FILES[b]`,
    the signature `ROI_SIGNATURES[b]` should be
        * `+`: if the ROI elements are not inside the boundary mesh
        * `-`: if the ROI elements are inside the boundary mesh
        * `*`: if ROI elements have no particular spatial relationship with the boundary


#### Example
An example is provided in [example/steps/spine_rois.py](example/steps/spine_rois.py).

This example imports a dendritic spine mesh reconstruction to STEPS,
and creates tetrahedral Regions of Interest (ROI)s using polyhedral surface boundaries.
The mesh with newly created ROIs is then exported to STEPS xml format,
and visualized using steps.visual module.

##### Prerequisite for the example
* [STEPS](http://steps.sourceforge.net/)
* pyqtgraph
* PyOpenGL
* PyQt5

##### Run the example
```
cd example/steps
python spine_rois.py
```
![visual](example/steps/visual.png)


### gmsh
A noteworthy difference between the `STEPS` and `gmsh` implementation is tets are duplicated
in the latter if they belong to multiple ROIs. This is a limitation related to the latest `gmsh` format

in `Python` interface, import the module
```python
import polyhedronROI.gmsh as pRg
```
either use the provided context manager or initialize gmsh API, then call the tagging method

```python
with pRg.gmsh_fixture():
        pRg.tag_mesh_entities(input_mesh=MESHIN, boundary_files=BOUNDARY_FILES, \
            roi_labels=ROI_LABELS, duplicate=DUPLICATE, output_mesh=MESHOUT, \
            interactive=INTERACTIVE)
# Alternatively
# gmsh.initialize()
# pRg.tag_mesh_entities(input_mesh=MESHIN, boundary_files=BOUNDARY_FILES, \
            # roi_labels=ROI_LABELS, duplicate=DUPLICATE, output_mesh=MESHOUT, \
            # interactive=INTERACTIVE)
# gmsh.finalize()

```
* `MESHIN`: A string containing the path to a gmsh file
* `DUPLICATE`: A optional boolean that allows for the duplication of tets (default `False`)
* `MESHOUT`: A optional string containing the path of the output gmsh file
* `INTERACTIVE`: A optional boolean that if set to `True` opens the gmsh GUI (default `False`)

##### Example
An example with few tests is provided in [example/gmsh/tests.py](example/gmsh/spine_rois.py).

##### Prerequisite for the example
* gmsh

##### Run the example
```
cd example/gmsh
python tests.py <cube|cube-overlap|overlap-ok|spine> [--interactive]
```
