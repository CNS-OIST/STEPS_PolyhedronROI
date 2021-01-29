import polyhedronROI.gmsh as pRg
import sys

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: " + sys.argv[0] + " <cube|cube-overlap|spine> [--interactive]")
        exit()

    test = sys.argv[1]

    if test == "cube":
        # cube test
        path = "./meshes/cube"
        mesh = path+"/cube.msh"
        boundary_files = [path+"/region_1.stl", path+"/region_2.stl"]
        roi_labels = {"roi1": "-*", "roi2": "*-"}

    elif test == "cube-overlap":
        # cube test with 2 ROIs overlapping
        path = "./meshes/cube"
        mesh = path+"/cube.msh"
        boundary_files = [path+"/region_2.stl", path+"/region_2.stl"]
        roi_labels = {"roi1": "-*", "roi2": "*-"}

    elif test == "spine":
        # spine ROIs
        path = "./meshes/spine"
        mesh = path+"/tets.msh"
        boundary_files = [path+"/ER.stl", path+"/PSD.stl"]
        roi_labels = {"ER": "-*",   # inside ER.stl, doesn't matter for PSD.stl
                      "PSD": "+-"}  # not inside ER.stl, inside PSD.stl
    else:
        print("Unknow test: "+test+ ". Options are <cube|cube-overlap|spine>")
        exit()

    interactive = '--interactive' in sys.argv

    with pRg.gmsh_fixture():
        pRg.tag_mesh_entities(mesh, boundary_files, roi_labels, \
            output_mesh="mod_"+test+".msh", interactive=interactive)
