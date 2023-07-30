# Importing the required modules
import matplotlib.pyplot as plt
import numpy as np

# Parameters for the problem
# All the data are in S.I units
K = 10      # Thermal conductivity (constant)
L = 0.5     # Length of the bar
A = 0.01    # Area of the bar (uniform)
Tb = 500    # Dirichlet boundary condition
Qa = 1000   # Von-Neumann boundary condition
Sc = 500    # Source term-1
Sp = -30    # Source term-2
N = 15      # No. of cells
dX = L/N    # 1-D Grid size
dV = A * dX # Differential volume

# Constants for the calculation
a1 = K * A / dX
a2 = K * A / dX
a3 = a1 + a2 - Sp * dV
a4 = 2 * a1 + a2 - Sp * dV
a5 = 2 * a1
a6 = 2 * a2
a7 = a1 - Sp * dV
b = Sc * dV
conv_error = 0.0001
value = 0

# Declaration of the variable required for the calculation
T = np.zeros(N)
T_old = np.zeros(N)
error = np.zeros(N)
iteration_no = 1
converge = False

# Creating the log file
# Change the directory as per the configuration
file = open(file = "/Users/nirupam/Desktop/Python/CFD (FVM)/log_file", mode='w')

# Using the jacobi method to solve the simultaneous equations
while converge == False:
    for i in range(N):
        if i == 0:
            T[i] = 1/a7*(a1*T[i+1] + Qa*A + b)
        elif i == N-1:
            T[i] = 1/a4*(a2*T[i-1] + a5*Tb + b)
        else:
            T[i] = 1/a3*(a1*T[i+1] + a2*T[i-1] + b)
    print(f"Iteration-{iteration_no} : {T}")
    file.write(f"Iteration-{iteration_no} : {T}\n\n")
    iteration_no += 1
    # Checking for convergence
    for j in range(N):
        error[j] = T[j] - T_old[j]
        T_old[j] = T[j]
        if (error[j] < conv_error):
            value += 1
    if value == N:
        converge = True
        print("\nThe solution has converged!\n")
        file.write("The solution has converged!\n\n")

# Interplolating the solution
Ta = (Qa / K) * (dX / 2) + T[0]
T = np.append(T, Tb)
T = np.insert(T, 0, Ta)

# Print the solution
print(T)

# Creating the spatial domain
# Plotting the solution
x = np.linspace(dX/2, L-dX/2, N)
x = np.insert(x, 0, 0)
x = np.append(x, L)
plt.xlim(0, L)
plt.plot(x, T, color='green', linestyle='--')
plt.show()

# Writing the final solution to the log file
file.write(f"X = {x}\n\n")
file.write(f"T = {T}")
file.close()
