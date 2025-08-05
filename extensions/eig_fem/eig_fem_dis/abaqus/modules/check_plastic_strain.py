def check_plastic_strain(filename, n_step, N_DD):
    """
    This function calculates the average plastic strain rate from a ParaDiS result file.
    
    Parameters:
        filename (str): Path to the input file.
        n_step (int): Current simulation step.
        N_DD (int): Number of dislocation dynamics steps to average over.
    
    Returns:
        float: Averaged plastic strain rate.
    """
    rate_plastic_strain = 0.0

    with open(filename, 'r') as file:
        # Skip lines before the target data
        for _ in range(n_step - N_DD):
            # next(file)
            line = file.readline()

        # Read and sum up the plastic strain rate values
        for _ in range(N_DD):
            line = file.readline()
            if not line:
                break  # Avoid error if file ends early
            parts = line.strip().split()
            if len(parts) >= 2:
                rate_plastic_strain += float(parts[1])

    rate_plastic_strain /= N_DD
    return rate_plastic_strain