import numpy as np
import stim
from Surface import surface_data

# Function to map qubit IDs to numbers with specific prefixes
def map_qubit_ids(matrix):
    qubit_map = {}
    counter = 0
    rows, cols = matrix.shape
    for i in range(rows):
        for j in range(cols):
            if isinstance(matrix[i, j], str) and matrix[i, j] != '#':
                qubit_id = matrix[i, j]
                if qubit_id.startswith('A'):
                    mapped_id = '1' + qubit_id[1:]
                elif qubit_id.startswith('B'):
                    mapped_id = '2' + qubit_id[1:]
                elif qubit_id.startswith('T'):
                    mapped_id = '3' + qubit_id[1:]
                elif qubit_id.startswith('Q'):
                    mapped_id = '4' + qubit_id[1:]
                else:
                    mapped_id = qubit_id
                
                if mapped_id not in qubit_map:
                        qubit_map[qubit_id] = int(mapped_id)
    
    return qubit_map

# Function to initialize qubits with coordinates
def initialize_qubits(circuit, cut_matrix, qubit_map):
    rows, cols = cut_matrix.shape
    for i in range(rows):
        for j in range(cols):
            if isinstance(cut_matrix[i, j], str) and cut_matrix[i, j] != '#':
                qubit_id = qubit_map[cut_matrix[i, j]]
                circuit += stim.Circuit(f"QUBIT_COORDS({i},{j}) {qubit_id}")
    return list(qubit_map.values())

# Function to reset qubits
def reset_qubits(circuit, qubit_ids, noise,noise2,r):
    high=[]
    low=[]
    for i in qubit_ids:
        if str(i).startswith('3'):
            high.append(i)
        else:
            low.append(i)

    if r == 0:
        circuit.append('R', qubit_ids)
        circuit.append('TICK')
        circuit.append('DEPOLARIZE1', low, noise)
        circuit.append('DEPOLARIZE1', high, noise)


    else:    
        circuit.append('DEPOLARIZE1', low, noise)
        circuit.append('DEPOLARIZE1', high, noise)



# Function to create Bell pairs
def create_bell_pairs(circuit, t1_ancillas, t2_ancillas, t3_ancillas, qubit_map,noise2):
    t1_ancilla_ids = [qubit_map[q] for q in t1_ancillas]
    t2_ancilla_ids = [qubit_map[q] for q in t2_ancillas]
    t3_ancilla_ids = [qubit_map[q] for q in t3_ancillas]
    
    # Apply H gates to T1 ancillas
    circuit.append("H", t1_ancilla_ids)
    circuit.append('DEPOLARIZE1',t1_ancilla_ids, noise2)
    circuit.append('TICK')
    
    # Apply CNOT gates from T1 to T2 ancillas
    all_cx_pairs_t1_t2 = []
    for t1, t2 in zip(t1_ancilla_ids, t2_ancilla_ids):
        all_cx_pairs_t1_t2.append([t1, t2])
        
    circuit.append("CX", [q for pair in all_cx_pairs_t1_t2 for q in pair])
    circuit.append('DEPOLARIZE2', [q for pair in all_cx_pairs_t1_t2 for q in pair], noise2)
    circuit.append('TICK')
    
    # Apply CNOT gates from T2 to T3 ancillas
    all_cx_pairs_t2_t3 = []
    for t2, t3 in zip(t2_ancilla_ids, t3_ancilla_ids):
        all_cx_pairs_t2_t3.append([t2, t3])
        
    circuit.append("CX", [q for pair in all_cx_pairs_t2_t3 for q in pair])
    circuit.append('DEPOLARIZE2', [q for pair in all_cx_pairs_t2_t3 for q in pair], noise2)
    circuit.append('TICK')

# Function to create Bell pairs, CAT STATE
def create_cat_state(circuit, t1_ancillas, t2_ancillas, t3_ancillas, qubit_map, noise2):
    t1_ancilla_ids = [qubit_map[q] for q in t1_ancillas]
    t2_ancilla_ids = [qubit_map[q] for q in t2_ancillas]
    t3_ancilla_ids = [qubit_map[q] for q in t3_ancillas]
    
    # Apply H gates to T3 ancillas
    circuit.append("H", t3_ancilla_ids)
    circuit.append('DEPOLARIZE1',t3_ancilla_ids, noise2)

    circuit.append('TICK')
    
    # Apply CNOT gates from T1 to T2 ancillas   (3 to 1)
    all_cx_pairs_t3_t1 = []
    for t3, t1 in zip(t3_ancilla_ids, t1_ancilla_ids):
        all_cx_pairs_t3_t1.append([t3, t1])
        
    circuit.append("CX", [q for pair in all_cx_pairs_t3_t1 for q in pair])
    circuit.append('DEPOLARIZE2', [q for pair in all_cx_pairs_t3_t1 for q in pair], noise2)
    circuit.append('TICK')

# Function to apply Hadamard gates to all X ancillas
def apply_hadamard_x_ancillas(circuit, x_ancillas, t1_ancillas, qubit_map,noise,noise2):
    t1_ids = set(qubit_map[q] for q in t1_ancillas)
    x_ancilla_ids = [qubit_map[q] for q in x_ancillas]
    circuit.append("H", x_ancilla_ids)
    high=[]
    low=[]
    for i in x_ancillas:
        if i.startswith('T'):
            high.append(qubit_map[i])
        else:
            low.append(qubit_map[i]) 
    circuit.append('DEPOLARIZE1',high, noise)
    circuit.append('DEPOLARIZE1',low, noise)

    circuit.append('TICK')

## Function to create stabilizers for a given direction
def create_stabilizers_direction(circuit, data, qubit_map, noise,noise2, direction):
    # Get diagonal neighbors and extract the given direction neighbors
    diagonal_neighbors = data['diagonal_neighbors']
    index = {'nw': 2, 'ne': 3, 'se': 4, 'sw': 5}[direction]
    
    # Initialize a string to collect all qubits involved in CX gates
    depolarize_string1 = f'DEPOLARIZE2({noise}) '
    depolarize_string2 = f'DEPOLARIZE2({noise2}) '

    # For Z ancillas, apply CX between the neighbor and the ancilla
    for ancilla_type, ancilla, *neighbors in diagonal_neighbors:
        neighbor = neighbors[index - 2]
        if ancilla_type == 'Z' and neighbor is not None and neighbor != '#':
            neighbor_id = qubit_map[neighbor]
            ancilla_id = qubit_map[ancilla]
            circuit.append("CX", [neighbor_id, ancilla_id])
            if str(ancilla_id).startswith('3'):
                depolarize_string2 += f'{neighbor_id} {ancilla_id} '
            else:
                depolarize_string1 += f'{neighbor_id} {ancilla_id} '

    
    # For X ancillas, apply CX between the ancilla and its neighbor
    for ancilla_type, ancilla, *neighbors in diagonal_neighbors:
        neighbor = neighbors[index - 2]
        if ancilla_type == 'X' and neighbor is not None and neighbor != '#':
            neighbor_id = qubit_map[neighbor]
            ancilla_id = qubit_map[ancilla]
            circuit.append("CX", [ancilla_id, neighbor_id])
            if str(ancilla_id).startswith('3'):
                depolarize_string2 += f'{ancilla_id} {neighbor_id} '
            else:
                depolarize_string1 += f'{ancilla_id} {neighbor_id} '
    # Apply DEPOLARIZE2 to all involved qubits
    circuit += stim.Circuit(f'{depolarize_string1.strip()}')
    circuit += stim.Circuit(f'{depolarize_string2.strip()}')

    circuit.append('TICK')



# Function to create all stabilizers
def create_stabilizers(circuit, data, qubit_map, noise,noise2):
    create_stabilizers_direction(circuit, data, qubit_map, noise,noise, 'se')
    create_stabilizers_direction(circuit, data, qubit_map, noise,noise, 'ne')
    create_stabilizers_direction(circuit, data, qubit_map, noise,noise2, 'sw')
    create_stabilizers_direction(circuit, data, qubit_map, noise,noise2, 'nw')


def apply_hadamard_t2(circuit, t2_ancillas, qubit_map,noise2):
    t2_ancilla_ids = [qubit_map[q] for q in t2_ancillas]
    circuit.append("H", t2_ancilla_ids)
    circuit.append('DEPOLARIZE1',t2_ancilla_ids, noise2)
    circuit.append('TICK')


# Function to apply Measure and Reset (MR) to all elements except Qs
def apply_measure_and_reset(circuit, data, qubit_map,noise,noise2):
    elements = data['T1_Ancillas'] + data['T2_Ancillas'] + data['T3_Ancillas'] + data['A_Ancillas'] + data['B_Ancillas']
    element_ids = [qubit_map[q] for q in elements]
    high_elements = data['T1_Ancillas'] + data['T2_Ancillas'] + data['T3_Ancillas']
    high_element_ids = [qubit_map[q] for q in high_elements]
    low_elements = data['A_Ancillas'] + data['B_Ancillas']
    low_element_ids = [qubit_map[q] for q in low_elements]

    circuit.append('DEPOLARIZE1',high_element_ids, noise)
    circuit.append('DEPOLARIZE1',low_element_ids, noise)
    circuit.append("MR", element_ids)
    circuit.append('DEPOLARIZE1',high_element_ids, noise)
    circuit.append('DEPOLARIZE1',low_element_ids, noise)

    circuit.append('TICK')
    
    return element_ids  # Return the order of measurements

# Function to add Z detectors for all Z ancillas (excluding T1)
def add_z_detectors(circuit, data, qubit_map, measurement_order, round_index):
    z_ancillas = [qubit_map[q] for q in data['Z_Ancillas'] if q not in data['T1_Ancillas']]
    t1_ancillas = {qubit_map[q]: q for q in data['T1_Ancillas']}
    t3_ancillas = {qubit_map[q]: q for q in data['T3_Ancillas']}
    
    total_mr = len(measurement_order)
    
    for z_ancilla in z_ancillas:
        # Initialize the recs string
        recs = " "
        # Find the index of the measurement in the measurement order
        index = measurement_order.index(z_ancilla)
        # Calculate the backwards reference
        backwards_reference = -(len(measurement_order) - index)
        
        if round_index == 0:
            recs += f" rec[{backwards_reference}]"
        else:
            recs += f" rec[{backwards_reference}] rec[{backwards_reference - total_mr}]"
        
            """         if z_ancilla in t3_ancillas:
            print(f"t3_ancillas are: {t3_ancillas}")
            # Replace the last digit of the T3 ancilla with '1' to find the corresponding T1 ancilla
            corresponding_t1_name = t3_ancillas[z_ancilla][:-1] + '1'
            corresponding_t1 = qubit_map[corresponding_t1_name]
            # Find the index of the corresponding T1 in the measurement order
            t1_index = measurement_order.index(corresponding_t1)
            # Calculate the backwards reference for the corresponding T1
            t1_backwards_reference = -(len(measurement_order) - t1_index)
            
            if round_index == 0:
                recs += f" rec[{t1_backwards_reference}] "
            else:
                recs += f" rec[{t1_backwards_reference}] rec[{t1_backwards_reference - total_mr}] "
            """            
            """             # Additional correction logic
            row_index = int(t3_ancillas[z_ancilla][1:-1]) - 1
            up_correction_row = row_index - 2
            down_correction_row = row_index + 2

            up_correction_name = f'T{up_correction_row + 1}2'
            down_correction_name = f'T{down_correction_row + 1}2'
            if up_correction_name in qubit_map:
                up_correction_ancilla = qubit_map[up_correction_name]
                up_correction_index = measurement_order.index(up_correction_ancilla)
                up_correction_backwards_reference = -(len(measurement_order) - up_correction_index)
                if round_index == 0:
                    recs += f" rec[{up_correction_backwards_reference}] "
                else:
                    recs += f" rec[{up_correction_backwards_reference}] "
            if down_correction_name in qubit_map:
                down_correction_ancilla = qubit_map[down_correction_name]
                down_correction_index = measurement_order.index(down_correction_ancilla)
                down_correction_backwards_reference = -(len(measurement_order) - down_correction_index)
                if round_index > 0:
                    recs += f"  rec[{down_correction_backwards_reference - total_mr}]" """

        # Append DETECTOR instruction
        circuit += stim.Circuit(f"DETECTOR{recs}")
    if round_index ==0:
        circuit.append("TICK")

# Function to add X detectors for all X ancillas (excluding T1)
def add_x_detectors(circuit, data, qubit_map, measurement_order, round_index):
    x_ancillas = [qubit_map[q] for q in data['X_Ancillas'] if q not in data['T1_Ancillas']]
    t1_ancillas = {qubit_map[q]: q for q in data['T1_Ancillas']}
    t3_ancillas = {qubit_map[q]: q for q in data['T3_Ancillas']}
    
    total_mr = len(measurement_order)
    
    for x_ancilla in x_ancillas:
        # Initialize the recs string
        recs = " "
        
        # Find the index of the measurement in the measurement order
        index = measurement_order.index(x_ancilla)
        # Calculate the backwards reference
        backwards_reference = -(len(measurement_order) - index)
        
        if round_index == 0:
            recs += f" rec[{backwards_reference}]"
        else:
            recs += f" rec[{backwards_reference}] rec[{backwards_reference - total_mr}]"
        
            """         if x_ancilla in t3_ancillas:
            # Replace the last digit of the T3 ancilla with '1' to find the corresponding T1 ancilla
            corresponding_t1_name = t3_ancillas[x_ancilla][:-1] + '1'
            corresponding_t1 = qubit_map[corresponding_t1_name]
            # Find the index of the corresponding T1 in the measurement order
            t1_index = measurement_order.index(corresponding_t1)
            # Calculate the backwards reference for the corresponding T1
            t1_backwards_reference = -(len(measurement_order) - t1_index)
            
            if round_index == 0:
                recs += f" rec[{t1_backwards_reference}] "
            else:
                recs += f" rec[{t1_backwards_reference}] rec[{t1_backwards_reference - total_mr}] "
              """            
            """             # Additional correction logic
            row_index = int(t3_ancillas[x_ancilla][1:-1]) - 1
            up_correction_row = row_index - 2
            down_correction_row = row_index + 2

            up_correction_name = f'T{up_correction_row + 1}2'
            down_correction_name = f'T{down_correction_row + 1}2'
            if up_correction_name in qubit_map:
                up_correction_ancilla = qubit_map[up_correction_name]
                up_correction_index = measurement_order.index(up_correction_ancilla)
                up_correction_backwards_reference = -(len(measurement_order) - up_correction_index)
                if round_index == 0:
                    recs += f" rec[{up_correction_backwards_reference}] "
                else:
                    recs += f" rec[{up_correction_backwards_reference}] "
            if down_correction_name in qubit_map:
                down_correction_ancilla = qubit_map[down_correction_name]
                down_correction_index = measurement_order.index(down_correction_ancilla)
                down_correction_backwards_reference = -(len(measurement_order) - down_correction_index)
                if round_index > 0:
                    recs += f"  rec[{down_correction_backwards_reference - total_mr}]" """

        # Append DETECTOR instruction
        circuit += stim.Circuit(f"DETECTOR{recs}")
    circuit.append("TICK")

# Function to apply "M" gate to all Q qubits and keep track of the measurement order
def measure_q_qubits(circuit, qubit_map, data, measurement_order,noise):
    q_qubits = [qubit_map[q] for q in data['Qubits']]
    circuit.append('DEPOLARIZE1',q_qubits, noise)

    circuit.append("M", q_qubits)
    circuit.append('TICK')
    measurement_order.extend(q_qubits)  # Track measurement order

# Function to add final detectors for all Z ancillas
def add_final_detectors(circuit, data, qubit_map, measurement_order):
    z_ancillas = [qubit_map[q] for q in data['Z_Ancillas']]
    t1_ancillas = {qubit_map[q]: q for q in data['T1_Ancillas']}
    t3_ancillas = {qubit_map[q]: q for q in data['T3_Ancillas']}
    diagonal_neighbors = data['diagonal_neighbors']

    for ancilla_type, ancilla, nw, ne, se, sw in diagonal_neighbors:
        if qubit_map[ancilla] in z_ancillas:
            recs = ""
            neighbors = [nw, ne, se, sw]
            # Include the ancilla itself
            ancilla_id = qubit_map[ancilla]
            index = measurement_order.index(ancilla_id)
            backwards_reference = -(len(measurement_order) - index)
            if qubit_map[ancilla] not in t1_ancillas:
                recs += f" rec[{backwards_reference}]"
            
                # Include only valid neighbors
                valid_neighbors = [n for n in neighbors if n is not None and n != '#']
                for neighbor in valid_neighbors:
                    neighbor_id = qubit_map[neighbor]
                    index = measurement_order.index(neighbor_id)
                    backwards_reference = -(len(measurement_order) - index)
                    recs += f" rec[{backwards_reference}]"
                
                    """                 if qubit_map[ancilla] in t3_ancillas:
                    # Handle T3 ancillas separately
                    corresponding_t1_name = t3_ancillas[qubit_map[ancilla]][:-1] + '1'
                    corresponding_t1 = qubit_map[corresponding_t1_name]
                    
                    # Include the corresponding T1 ancilla itself
                    corresponding_t1_id = qubit_map[corresponding_t1_name]
                    index = measurement_order.index(corresponding_t1_id)
                    backwards_reference = -(len(measurement_order) - index)
                    recs += f" rec[{backwards_reference}]"
                    
                    # Find the neighbors of the corresponding T1 ancilla
                    corresponding_t1_neighbors = next((nw, ne, se, sw) for ancilla_type, anc, nw, ne, se, sw in diagonal_neighbors if anc == corresponding_t1_name)
                    valid_t1_neighbors = [n for n in corresponding_t1_neighbors if n is not None and n != '#']
                    
                    for neighbor in valid_t1_neighbors:
                        neighbor_id = qubit_map[neighbor]
                        index = measurement_order.index(neighbor_id)
                        backwards_reference = -(len(measurement_order) - index)
                        recs += f" rec[{backwards_reference}]"
                     """
                    """                     # Additional correction logic
                    row_index = int(t3_ancillas[qubit_map[ancilla]][1:-1]) - 1
                    down_correction_row = row_index + 2

                    down_correction_name = f'T{down_correction_row + 1}2'
                    if down_correction_name in qubit_map:
                        down_correction_ancilla = qubit_map[down_correction_name]
                        down_correction_index = measurement_order.index(down_correction_ancilla)
                        down_correction_backwards_reference = -(len(measurement_order) - down_correction_index)
                        recs += f" rec[{down_correction_backwards_reference}]" """

                # Append DETECTOR instruction
                circuit += stim.Circuit(f"DETECTOR{recs}")
    circuit.append("TICK")


def include_observable(circuit, data, qubit_map, measurement_order):
    # Get the penultimate row
    penultimate_row = data['cut_matrix'][-2]
    
    # Extract Q qubits in the penultimate row
    q_qubits = [qubit_map[q] for q in penultimate_row if isinstance(q, str) and q.startswith('Q')]
    
    # Create the list of references
    references = []
    for q in q_qubits:
        index = measurement_order.index(q)
        backwards_reference = -(len(measurement_order) - index)
        references.append(backwards_reference)
    
    # Create the OBSERVABLE_INCLUDE line
    recs = " ".join([f"rec[{ref}]" for ref in references])
    circuit += stim.Circuit(f"OBSERVABLE_INCLUDE(0) {recs}")
    circuit.append("TICK")


# Main function to generate the circuit
def generate_circuit(data, noise,noise2, rounds):
    circuit = stim.Circuit()
    
    # Map qubit IDs to numbers
    qubit_map = map_qubit_ids(data['cut_matrix'])
    
    # Initialize qubits and get their IDs
    qubit_ids = initialize_qubits(circuit, data['cut_matrix'], qubit_map)
    
    # Reset qubits

    for i in range(rounds):
        reset_qubits(circuit, qubit_ids, noise,noise2,i)

        #CAT STATE PREP
        #create_cat_state(circuit, data['T1_Ancillas'], data['T2_Ancillas'], data['T3_Ancillas'], qubit_map, noise2)
        # Apply Hadamard gates to all X ancillas
        apply_hadamard_x_ancillas(circuit, data['X_Ancillas'], data['T1_Ancillas'], qubit_map,noise,noise2)

        # Create stabilizers
        # Function to create all stabilizers
        create_stabilizers_direction(circuit, data, qubit_map, noise,noise, 'se')
        create_stabilizers_direction(circuit, data, qubit_map, noise,noise, 'ne')
        create_stabilizers_direction(circuit, data, qubit_map, noise,noise2, 'sw')
        create_stabilizers_direction(circuit, data, qubit_map, noise,noise2, 'nw')

        apply_hadamard_x_ancillas(circuit, data['X_Ancillas'], data['T1_Ancillas'], qubit_map,noise,noise2)

        # Apply Hadamard gates to all T2 ancillas
        #apply_hadamard_t2(circuit, data['T2_Ancillas'], qubit_map,noise2)

        # Apply Measure and Reset to all elements except Qs
        measurement_order = apply_measure_and_reset(circuit, data, qubit_map,noise,noise2)

        if i== 0:
            # Add Z detectors
            add_z_detectors(circuit, data, qubit_map, measurement_order,i)
        else:
            add_z_detectors(circuit, data, qubit_map, measurement_order,i)
            # Add X detectors
            add_x_detectors(circuit, data, qubit_map, measurement_order, i)
   
    measure_q_qubits(circuit, qubit_map, data, measurement_order,noise)

    # Add final detectors
    add_final_detectors(circuit, data, qubit_map, measurement_order)
    
    # Include observable
    include_observable(circuit, data, qubit_map, measurement_order)

    return circuit


