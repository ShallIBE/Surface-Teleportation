import numpy as np
import svgwrite
from IPython.display import SVG, display


def generate_surface_code_matrix_with_shifts(n, m):
    rows = 2 * n + 1
    cols = 2 * m + 1
    
    # Initialize the matrix with '#' (representing empty spaces)
    matrix = np.full((rows, cols), '#', dtype=object)
    
    q_counter = 1  # Counter for qubits
    b_counter = 1  # Counter for bulk ancillas
    a_counter = 1  # Counter for auxiliary ancillas
    
    # Place physical qubits (Q)
    qubit_positions = {}
    for i in range(1, rows, 2):
        for j in range(1, cols, 2):
            matrix[i, j] = f'Q{q_counter}'
            qubit_positions[f'Q{q_counter}'] = (i, j)
            q_counter += 1
    
    # Place bulk ancillas (B)
    for i in range(2, rows-2, 2):
        for j in range(2, cols-2, 2):
            matrix[i, j] = f'B{b_counter}'
            b_counter += 1

    # Create the boundary qubit list
    boundary_qubits = []

    # Top row (left to right)
    for j in range(1, cols, 2):
        boundary_qubits.append(f'Q{matrix[1, j][1:]}')

    # Right column (top to bottom)
    for i in range(3, rows, 2):
        boundary_qubits.append(f'Q{matrix[i, cols - 2][1:]}')

    # Bottom row (right to left)
    for j in range(cols - 2, 0, -2):
        boundary_qubits.append(f'Q{matrix[rows - 2, j][1:]}')

    # Left column (bottom to top)
    for i in range(rows - 4, 0, -2):
        boundary_qubits.append(f'Q{matrix[i, 1][1:]}')

    # Remove duplicates from the boundary list
    boundary_qubits = list(dict.fromkeys(boundary_qubits))

    
    # Place auxiliary ancillas (A) between pairs of boundary qubits with shifts
    for k in range(1, len(boundary_qubits), 2):
        q1 = boundary_qubits[k - 1]
        q2 = boundary_qubits[k]

        # Get positions of q1 and q2
        q1_pos = qubit_positions[q1]
        q2_pos = qubit_positions[q2]
        
        if q1_pos[0] == q2_pos[0]:  # Same row
            ancilla_position = (q1_pos[0], (q1_pos[1] + q2_pos[1]) // 2)
            if ancilla_position[0] == 1:
                ancilla_position = (0, ancilla_position[1])  # Shift to top row
            elif ancilla_position[0] == rows - 2:
                ancilla_position = (rows - 1, ancilla_position[1])  # Shift to bottom row
        else:  # Same column
            ancilla_position = ((q1_pos[0] + q2_pos[0]) // 2, q1_pos[1])
            if ancilla_position[1] == 1:
                ancilla_position = (ancilla_position[0], 0)  # Shift to left column
            elif ancilla_position[1] == cols - 2:
                ancilla_position = (ancilla_position[0], cols - 1)  # Shift to right column

        # Place auxiliary ancilla (A)
        matrix[ancilla_position] = f'A{a_counter}'
        a_counter += 1
    
    return matrix

def make_vertical_cut(matrix, c, spacing):
    if c == 0:
        return matrix
    rows, cols = matrix.shape
    keep_columns = 2 * c
    if keep_columns >= cols:
        raise ValueError("Not enough columns in the matrix to make the cut.")

    # Create a new matrix with the additional columns
    new_matrix = np.full((rows, cols + spacing), '#', dtype=object)

    # Copy the initial unchanged part of the matrix
    for i in range(rows):
        new_matrix[i, :keep_columns] = matrix[i, :keep_columns]

    # Determine the column to modify and copy
    modify_column = keep_columns

    # Change IDs in the modify column to begin with 'T'
    for i in range(rows):
        if isinstance(matrix[i, modify_column], str) and matrix[i, modify_column][0] in {'A', 'B'}:
            new_matrix[i, modify_column] = 'T' + f'{i+1}' 
        else:
            new_matrix[i, modify_column] = matrix[i, modify_column]

    # Make 3 copies of the modify column
    for i in range(rows):
        if isinstance(matrix[i, modify_column], str) and matrix[i, modify_column][0] in {'A', 'B'}:
            og= new_matrix[i, modify_column]
            for k in range(2):
                new_matrix[i, modify_column + spacing -k] = og + f'{3-k}'
            new_matrix[i, modify_column] = new_matrix[i, modify_column] +'1'
    # Copy the rest of the matrix, adjusting for the added columns
    for i in range(rows):
        new_matrix[i, modify_column + spacing+1 :cols + spacing] = matrix[i, modify_column + 1:cols]

    return new_matrix

def classify_ancillas(matrix):
    rows, cols = matrix.shape
    B_Ancillas = []
    A_Ancillas = []
    T3_Ancillas = []
    T2_Ancillas = []
    T1_Ancillas = []
    Z_Ancillas = []
    X_Ancillas = []
    Qubits = []

    for i in range(rows):
        for j in range(cols):
            if isinstance(matrix[i, j], str) and matrix[i, j][0] in {'A'}:
                A_Ancillas.append(matrix[i, j])
            if isinstance(matrix[i, j], str) and matrix[i, j][0] in {'B'}:
                B_Ancillas.append(matrix[i, j])
            if isinstance(matrix[i, j], str) and matrix[i, j][0] in {'T'} and matrix[i, j][-1] == '3':
                T3_Ancillas.append(matrix[i, j])
            if isinstance(matrix[i, j], str) and matrix[i, j][0] in {'T'} and matrix[i, j][-1] == '2':
                T2_Ancillas.append(matrix[i, j])
            if isinstance(matrix[i, j], str) and matrix[i, j][0] in {'T'} and matrix[i, j][-1] == '1':
                T1_Ancillas.append(matrix[i, j])
            if isinstance(matrix[i, j], str) and matrix[i, j][0] in {'Q'}:
                Qubits.append(matrix[i, j])

    X_Ancillas = []
    Z_Ancillas = []
    for i in range(2, rows - 1):
        counter = ((i / 2) - 1) % 2
        for j in range(cols):
            if isinstance(matrix[i, j], str) and matrix[i, j][0] in {'B'}:
                if counter % 2 == 0:
                    X_Ancillas.append(matrix[i, j])
                    counter += 1
                else:
                    Z_Ancillas.append(matrix[i, j])
                    counter += 1
            elif isinstance(matrix[i, j], str) and matrix[i, j][0] in {'T'} and matrix[i, j][-1] == '1':
                if counter % 2 == 0:
                    X_Ancillas.append(matrix[i, j])
                    counter += 1
                else:
                    Z_Ancillas.append(matrix[i, j])
                    counter += 1

    for i in (0, rows - 1):
        for j in range(cols):
            if isinstance(matrix[i, j], str) and matrix[i, j][0] in {'A'}:
                Z_Ancillas.append(matrix[i, j])
            elif isinstance(matrix[i, j], str) and matrix[i, j][0] in {'T'} and matrix[i, j][-1] == '1':
                Z_Ancillas.append(matrix[i, j])

    for i in range(0, rows):
        for j in (0, cols - 1):
            if isinstance(matrix[i, j], str) and matrix[i, j][0] in {'A'}:
                X_Ancillas.append(matrix[i, j])
            elif isinstance(matrix[i, j], str) and matrix[i, j][0] in {'T'} and matrix[i, j][-1] == '1':
                X_Ancillas.append(matrix[i, j])

    for element in Z_Ancillas.copy():  # Use a copy to avoid modifying the list while iterating
        if element.startswith('T') and element.endswith('1'):
            # Create a new element with the last digit changed to '3'
            new_element = element[:-1] + '3'
            # Append the new element to Z_Ancillas
            Z_Ancillas.append(new_element)

    for element in X_Ancillas.copy():  # Use a copy to avoid modifying the list while iterating
        if element.startswith('T') and element.endswith('1'):
            # Create a new element with the last digit changed to '3'
            new_element = element[:-1] + '3'
            # Append the new element to X_Ancillas
            X_Ancillas.append(new_element)

    return X_Ancillas, Z_Ancillas, B_Ancillas, T3_Ancillas, T2_Ancillas, T1_Ancillas, A_Ancillas, Qubits

def find_diagonal_neighbors_with_type(matrix, X_Ancillas, Z_Ancillas):
    rows, cols = matrix.shape
    diagonal_neighbors = []

    for i in range(rows):
        for j in range(cols):
            if isinstance(matrix[i, j], str) and matrix[i, j][0] in {'A', 'B', 'T'}:
                neighbors = []
                # North-West
                if i > 0 and j > 0:
                    neighbors.append(matrix[i-1, j-1])
                else:
                    neighbors.append(None)
                
                # North-East
                if i > 0 and j < cols - 1:
                    neighbors.append(matrix[i-1, j+1])
                else:
                    neighbors.append(None)
                
                # South-East
                if i < rows - 1 and j < cols - 1:
                    neighbors.append(matrix[i+1, j+1])
                else:
                    neighbors.append(None)
                
                # South-West
                if i < rows - 1 and j > 0:
                    neighbors.append(matrix[i+1, j-1])
                else:
                    neighbors.append(None)
                
                # Determine type and create a tuple with the element and its neighbors
                ancilla_type = 'X' if matrix[i, j] in X_Ancillas else 'Z' if matrix[i, j] in Z_Ancillas else ''
                diagonal_neighbors.append((ancilla_type, matrix[i, j], *neighbors))
    
    return diagonal_neighbors

def draw_surface_code_svg(matrix, x_ancillas, z_ancillas, n, m):
    rows, cols = matrix.shape
    cell_size = 40
    legend_height = 80  # Increased space for the legend at the bottom
    legend_cell_size = 10
    title_height = 40   # Space for the title at the top
    dwg_height = rows * cell_size + legend_height + title_height
    dwg = svgwrite.Drawing(size=(cols * cell_size, dwg_height))

    # Add title
    title = f"Modular Surface Code {n}x{m}"
    dwg.add(dwg.text(title, insert=(cols * cell_size / 2, title_height / 2),
                     text_anchor="middle", font_size=20, dominant_baseline="middle"))

    # Draw the matrix
    for i in range(rows):
        for j in range(cols):
            if matrix[i, j] != '#':
                # Determine the fill color
                if matrix[i, j].startswith('Q'):
                    fill_color = 'grey'
                elif matrix[i, j].startswith('T') and matrix[i, j].endswith('1'):
                    fill_color = 'green'
                elif matrix[i, j] in x_ancillas:
                    fill_color = 'red'
                elif matrix[i, j] in z_ancillas:
                    fill_color = 'blue'
                else:
                    fill_color = 'green'
                # Draw the rectangle
                dwg.add(dwg.rect((j * cell_size, i * cell_size + title_height), (cell_size, cell_size), fill=fill_color, stroke='black'))
                # Add the text
                dwg.add(dwg.text(matrix[i, j], insert=((j + 0.5) * cell_size, (i + 0.5) * cell_size + title_height),
                                 text_anchor="middle", dominant_baseline="middle", font_size=12, fill='white'))

    # Add legend
    legend_y = rows * cell_size + title_height + 10
    legend_items = [
        ('grey', 'Physical Qubits'),
        ('red', 'X Ancilla'),
        ('blue', 'Z Ancilla'),
        ('green', 'Bell Pair')
    ]
    for idx, (color, label) in enumerate(legend_items):
        x_start = 10
        y_start = legend_y + idx * (legend_cell_size + 5)
        dwg.add(dwg.rect((x_start, y_start), (legend_cell_size, legend_cell_size), fill=color, stroke='black'))
        dwg.add(dwg.text(label, insert=(x_start + legend_cell_size + 10, y_start + legend_cell_size / 2),
                         text_anchor="start", dominant_baseline="middle", font_size=12, fill='black'))

    return dwg.tostring()




def surface_data(n, m, c=0, spacing=3, display_svg=False):
    initial_matrix = generate_surface_code_matrix_with_shifts(n, m)
    cut_matrix = make_vertical_cut(initial_matrix, c, spacing)
    X_Ancillas, Z_Ancillas, B_Ancillas, T3_Ancillas, T2_Ancillas, T1_Ancillas, A_Ancillas, Qubits = classify_ancillas(cut_matrix)
    diagonal_neighbors = find_diagonal_neighbors_with_type(cut_matrix, X_Ancillas, Z_Ancillas)

    result = {
        "initial_matrix": initial_matrix,
        "cut_matrix": cut_matrix,
        "X_Ancillas": X_Ancillas,
        "Z_Ancillas": Z_Ancillas,
        "B_Ancillas": B_Ancillas,
        "T3_Ancillas": T3_Ancillas,
        "T2_Ancillas": T2_Ancillas,
        "T1_Ancillas": T1_Ancillas,
        "A_Ancillas": A_Ancillas,
        "Qubits": Qubits,
        "diagonal_neighbors": diagonal_neighbors
    }

    if display_svg:
        initial_svg = draw_surface_code_svg(initial_matrix, X_Ancillas, Z_Ancillas, n, m)
        cut_svg = draw_surface_code_svg(cut_matrix, X_Ancillas, Z_Ancillas,n,m)
        display(SVG(initial_svg))
        display(SVG(cut_svg))

    return result

