def generate_input_ABAQUS_stressControl(foldername_restart, foldername_ABAQUS, ep_rate, paradis_timeNow, N_DD, critical_eps_rate):
    import numpy as np
    import os
    import re

    # Read restart.data file and find keyword line
    restart_data_path = os.path.join(foldername_restart, 'restart.data')
    with open(restart_data_path, 'r') as file:
        lines = file.readlines()

    keyword = 'Secondary lines:'
    num_data_start = -1
    for i, line in enumerate(lines[:50]):
        if keyword in line:
            num_data_start = i
            break

    # Parse node data from restart.data
    rn = []

    # Simulate fscanf() style line-by-line parsing
    for line in lines[num_data_start+1:]:
        if not line.strip():
            break  # end of data

        # Try splitting using comma + whitespace combination
        parts = re.split(r'[,\s]+', line.strip())
        parts = [p for p in parts if p != '']  # remove empty strings

        # If line has exactly 7 numbers, treat as primary line
        if len(parts) == 7:
            data = [float(x) for x in parts]
            rn.append([data[2], data[3], data[4], data[6], int(data[1]), int(data[5])])
        # If line has more than 7 numbers, treat as secondary line
        elif len(parts) >= 13:
            data = [float(x) for x in parts]
            row = [
                data[-13],  # pos x
                data[-12],  # pos y
                data[-11],  # pos z
                data[-9],   # flag
                int(data[-14]),  # node #
                int(data[-10])   # domain
            ]
            rn.append(row)

    rn = np.array(rn)

    # Read restart.cn for current time
    restart_cn_path = os.path.join(foldername_restart, 'restart.cn')
    paradis_timeUpdate = None
    with open(restart_cn_path, 'r') as file:
        for i in range(50):
            line = file.readline()
            if 'maxstep' in line:
                break
        for line in file:
            match = re.search(r'timeNow\s*=\s*([-\d.eE+]+)', line)
            if match:
                paradis_timeUpdate = float(match.group(1))
                break

    # Calculate Dt
    paradis_Dt = paradis_timeUpdate - paradis_timeNow

    # Apply load
    loadMode = 2
    stress_inc = 0.05  # Normalized torque
    FEM_Du = stress_inc * N_DD if ep_rate <= critical_eps_rate else 0.0

    # Write outputs
    os.makedirs(foldername_ABAQUS, exist_ok=True)

    with open(os.path.join(foldername_ABAQUS, 'nodes_num.txt'), 'w') as f:
        f.write(f'{rn.shape[0]}')

    with open(os.path.join(foldername_ABAQUS, 'FEM_deltaU.txt'), 'w') as f:
        f.write(f'{loadMode} {FEM_Du:.6e}')

    with open(os.path.join(foldername_ABAQUS, 'rn_matrix.txt'), 'w') as f:
        for row in rn:
            f.write(
                f"{int(row[4])} {int(row[3])} {int(row[5])} "
                f"{row[0]:>12.10g} {row[1]:>12.10g} {row[2]:>12.10g}\n"
            )

    return paradis_Dt, FEM_Du
