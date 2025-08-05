import os
import re
import numpy as np


def import_elements_nodes(foldername_tests, foldername_Job,
                          foldername_ABAQUS,
                          jobname_head, abaqus_input_filename,
                          element_type, surface_element_set_name,
                          format_spec_elset, elset_col_num):
    """
    Parse ABAQUS .inp file to extract nodes, elements, and surface element sets.
    Returns a dict with:
      Nodes: (Num_Nodes, 3) array of nodal coords   
      Elements: (Num_Elements, N) array of element connectivity
      SurfaceElementSet: list of surface element indices
      SurfaceFlags: (Num_Elements,) array, 6 for surface elements, 0 otherwise
      Elements_Min, Elements_Max, Elements_Center: per-element bounding box and center
    Also writes OriginalCoordNodes.txt matching MATLAB's output format.
    """
    # Build input path
    inp_path = os.path.join('..','..',
                            foldername_tests, foldername_Job,
                            foldername_ABAQUS,
                            jobname_head + abaqus_input_filename + '.inp')

    # Read file lines without newline
    with open(inp_path, 'r') as f:
        ABAQUS_input_file = [line.rstrip('\n') for line in f]

    # 섹션 인덱스 추출
    idxNode = next(i for i, l in enumerate(ABAQUS_input_file) if '*Node' in l) + 1
    idxElement = next(i for i, l in enumerate(ABAQUS_input_file) if '*Element' in l) + 1
    idxNode_end = idxElement - 2
    idxElement_end = next(i for i, l in enumerate(ABAQUS_input_file[idxElement:], start=idxElement)
                          if l.strip().startswith('*')) - 1

    # SurfaceElementSet 인덱스 추출
    idxSurfaceElementSet = next(i for i, l in enumerate(ABAQUS_input_file)
                                 if surface_element_set_name in l) + 1
    idxSurfaceElementSet_end = next(i for i, l in enumerate(ABAQUS_input_file[idxSurfaceElementSet:], start=idxSurfaceElementSet)
                                     if l.strip().startswith('*')) - 1

    # Parse Nodes
    raw_nodes = []  # list of [id, x, y, z]
    for l in ABAQUS_input_file[idxNode:idxNode_end+1]:
        parts = re.split(r'[,\s]+', l.strip())
        nums = [float(p) for p in parts if p]
        raw_nodes.append(nums)
    raw_nodes = np.array(raw_nodes)
    raw_nodes = raw_nodes[raw_nodes[:, 0].argsort()]
    coords = raw_nodes[:, 1:4]  # Num_Nodes x 3
    Num_Nodes = coords.shape[0]

    # Parse Elements
    raw_elems = []
    for l in ABAQUS_input_file[idxElement:idxElement_end+1]:
        parts = re.split(r'[,\s]+', l.strip())
        nums = [int(p) for p in parts if p]
        raw_elems.append(nums)
    raw_elems = np.array(raw_elems)
    raw_elems = raw_elems[raw_elems[:, 0].argsort()]
    elems = raw_elems[:, 1:]   # Num_Elements x nodes_per_elem
    Num_Elements = elems.shape[0]
    Elements_col_num = raw_elems.shape[1]

    # Parse Surface Element Set
    surf_ids = []
    for l in ABAQUS_input_file[idxSurfaceElementSet:idxSurfaceElementSet_end+1]:
        parts = re.split(r'[,\s]+', l.strip())
        for p in parts:
            if p.isdigit():
                surf_ids.append(int(p))
    surf_ids = sorted(surf_ids)

    # Build surface flags
    surf_flag = np.zeros(Num_Elements, dtype=int)
    for sid in surf_ids:
        if 1 <= sid <= Num_Elements:
            surf_flag[sid-1] = 6

    # Write OriginalCoordNodes.txt in MATLAB format
    out_path = os.path.join('..', '..',
                            foldername_tests, foldername_Job,
                            foldername_ABAQUS, 'OriginalCoordNodes.txt')

    with open(out_path, 'w') as f:
        for j in range(Num_Elements):
            eid = int(raw_elems[j, 0])
            node_ids = raw_elems[j, 1:]
            coords_list = [coords[nid-1] for nid in node_ids]
            if Elements_col_num == 5:
                line_vals = [eid] + [c for coord in coords_list for c in coord] * 2
            else:
                line_vals = [eid] + [c for coord in coords_list for c in coord]
            fmt = '%d ' + ' '.join(['%.8e'] * (len(line_vals)-1)) + '\n'
            f.write(fmt % tuple(line_vals))

    # 추가: 요소 최소/최대/중심 계산
    Nodes_only = coords
    Elements_only = elems.astype(int)
    Elements_Min = np.zeros((Num_Elements, 3))
    Elements_Max = np.zeros((Num_Elements, 3))
    Elements_Center = np.zeros((Num_Elements, 3))
    for i in range(Num_Elements):
        # 요소 타입에 따라 노드 개수 구분
        if 'C3D4' in element_type:
            node_indices = Elements_only[i][:4]
        else:
            node_indices = Elements_only[i][:8]
        node_coords = Nodes_only[node_indices - 1]
        Elements_Min[i] = np.min(node_coords, axis=0)
        Elements_Max[i] = np.max(node_coords, axis=0)
        if 'C3D8' in element_type:
            Elements_Center[i] = np.mean(node_coords, axis=0)

    Num_SurfaceElement = len(surf_ids)

    return (
        Nodes_only,
        Elements_only,
        surf_ids,
        surf_flag,
        Elements_Min,
        Elements_Max,
        Elements_Center,
        Num_Nodes,
        Num_Elements,
        Num_SurfaceElement
    )
