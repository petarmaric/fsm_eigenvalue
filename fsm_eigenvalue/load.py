from beam_integrals.beam_types import BaseBeamType
import networkx as nx
import numpy as np
import yaml


from . import ASTIFF_BLOCK_SIZE


def parse_data_file(data_file):
    with open(data_file, 'r') as fp:
        return yaml.load(fp)

def get_beam_type_id(geometry):
    beam_type_id = geometry['beam_type_id']
    assert beam_type_id in BaseBeamType.plugins.valid_ids
    return beam_type_id

def linspace_with_step(start, stop, step, **kwargs):
    # Workaround until https://github.com/numpy/numpy/issues/630 is resolved
    num_samples = np.round((stop - start) / step + 1)
    return np.linspace(start, stop, num_samples, **kwargs)

def get_search_space_iterations(search_space):
    def dynamic_dtype(key):
        return int if key == 'm' else float # only 'mode' is int

    return {
        key: linspace_with_step(start, stop, step, dtype=dynamic_dtype(key))
        for key, (start, stop, step) in search_space.items()
    }

def get_transformation_matrix(dx, dz, b):
    sin_a = dz/b
    cos_a = dx/b

    # As per eq. 3.58,3.59 from [Milasinovic1997]
    R = np.asmatrix(np.eye(8))
    R[0, 0] =  cos_a
    R[2, 0] = -sin_a
    R[0, 2] =  sin_a
    R[2, 2] =  cos_a
    R[4, 4] =  cos_a
    R[6, 4] = -sin_a
    R[4, 6] =  sin_a
    R[6, 6] =  cos_a

    return R

def get_nodal_graph(geometry):
    # Mathematical graph of all nodal lines connected by finite strips
    nodal_graph = nx.DiGraph()

    # Add nodal lines
    for node_id, (x, z) in geometry['nodal_lines'].items():
        nodal_graph.add_node(node_id, x=x, z=z)

    # Add finite strips
    for strip_id, (node1_id, node2_id, material_id) in enumerate(geometry['finite_strips'], start=1):
        dx = nodal_graph.node[node2_id]['x'] - nodal_graph.node[node1_id]['x']
        dz = nodal_graph.node[node2_id]['z'] - nodal_graph.node[node1_id]['z']
        b = np.sqrt(dx**2 + dz**2) # [mm] strip width

        R = get_transformation_matrix(dx, dz, b) # global<->local coordinates transformation matrix

        # Deduced from Fortran block 83:94
        astiff_blocks = (node1_id-1, node2_id-1) # Python counts from 0
        astiff_fill_indices = []
        for row in xrange(2):
            for col in xrange(2):
                astiff_row_start = ASTIFF_BLOCK_SIZE * astiff_blocks[row]
                astiff_row_end   = ASTIFF_BLOCK_SIZE + astiff_row_start
                astiff_col_start = ASTIFF_BLOCK_SIZE * astiff_blocks[col]
                astiff_col_end   = ASTIFF_BLOCK_SIZE + astiff_col_start
                astiff_indices = (
                    slice(astiff_row_start, astiff_row_end),
                    slice(astiff_col_start, astiff_col_end)
                )

                segment_row_start = ASTIFF_BLOCK_SIZE * row
                segment_row_end   = ASTIFF_BLOCK_SIZE + segment_row_start
                segment_col_start = ASTIFF_BLOCK_SIZE * col
                segment_col_end   = ASTIFF_BLOCK_SIZE + segment_col_start
                segment_indices = (
                    slice(segment_row_start, segment_row_end),
                    slice(segment_col_start, segment_col_end)
                )

                astiff_fill_indices.append((astiff_indices, segment_indices))

        label = "(%d)" % strip_id

        edge_data_keys = 'material_id, b, R, astiff_fill_indices, label'.split(', ')
        nodal_graph.add_edge(node1_id, node2_id, dict(zip(
            edge_data_keys,
            (material_id, b, R, astiff_fill_indices, label)
        )))

    # Disallow further changes to the nodal_graph
    nx.freeze(nodal_graph)

    # Cache the traversal through strips for performance
    strip_data = nodal_graph.edges(data=True)

    return nodal_graph, strip_data

def precompute_material_properties(materials):
    for material in materials.values():
        material['ro'] /= 10**9 # convert mass density from [kg/m**3] to [kg/mm**3] before calc

        data_keys = 'E_x, E_y, mu_x, mu_y, G_xy'.split(', ')
        E_x , E_y, mu_x, mu_y, G_xy = (material[k] for k in data_keys)

        # As per eq. 2.20 from [Milasinovic1997]
        mu_xy = 1. - mu_x*mu_y
        material['K_x']  = E_x / mu_xy
        material['K_y']  = E_y / mu_xy
        material['K_1']  = mu_x * material['K_y']
        material['K_xy'] = G_xy

    return materials

def get_astiff_shape(nodal_graph):
    astiff_size = ASTIFF_BLOCK_SIZE * nodal_graph.number_of_nodes()
    astiff_shape = (astiff_size, astiff_size)
    return astiff_shape

def load_data_from(data_file):
    input_data = parse_data_file(data_file)

    beam_type_id = get_beam_type_id(input_data['geometry'])
    search_space = get_search_space_iterations(input_data['search_space'])
    nodal_graph, strip_data = get_nodal_graph(input_data['geometry'])
    materials = precompute_material_properties(input_data['materials'])
    astiff_shape = get_astiff_shape(nodal_graph)

    return (
        beam_type_id,
        search_space,
        nodal_graph,
        strip_data,
        materials,
        astiff_shape,
    )
