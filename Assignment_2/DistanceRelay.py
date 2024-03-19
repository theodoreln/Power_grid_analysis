import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Wedge

# Line Impedances
Z_12 = 7 + 70j
Z_23 = 6 + 60j
Z_24 = 8 + 80j

#Zone impedances 
Z_z1 = 0.8*Z_12
Z_z2 = 1.2 *Z_12 
Z_z3 = Z_12 + 1.2*Z_24 
# Mho relay settings 
VTR = 4500
CTR = 500

Z_r1 = CTR/VTR * Z_z1
Z_r2 = CTR/VTR * Z_z2
Z_r3 = CTR/VTR * Z_z3

# Radius
Z1_abs = np.absolute(Z_r1)
Z2_abs = np.absolute(Z_r2)
Z3_abs = np.absolute(Z_r3)

# R and X values of Z
R1 = np.real(Z_r1)
R2 = np.real(Z_r2)
R3 = np.real(Z_r3)

X1 = np.imag(Z_r1)
X2 = np.imag(Z_r2)
X3 = np.imag(Z_r3)

# Plot
R_values = range(-20, 20)
X_values = range(-20, 20)

# Create a new figure
fig, axes = plt.subplots(dpi = 400)

# Plot the x-axis
plt.plot([-20, 20], [0, 0], 'k-', linewidth=1)

# Plot the y-axis
plt.plot([0, 0], [-20, 20], 'k-', linewidth=1)


# Set x and y limits
plt.xlim(-20, 20)
plt.ylim(-20, 20)

# Set x and y labels
plt.xlabel('Re(Z)')
plt.ylabel('Im(Z)')

# Grey out area that#s behind the relay
theta=-np.arccos(X3/Z3_abs)*180/np.pi
grey_area = Wedge((0,0), Z3_abs, 180+theta, theta, fc='lightgrey', edgecolor='k')
axes.add_artist(grey_area)
# Plot vectors
plt.quiver(0, 0, R3, X3, angles='xy', scale_units='xy', scale=1, color='b', label='Z_r3')
plt.quiver(0, 0, R2, X2, angles='xy', scale_units='xy', scale=1, color='m', label='Z_r2')
plt.quiver(0, 0, R1, X1, angles='xy', scale_units='xy', scale=1, color='r', label='Z_r1')

# Plot circles
circle1 = plt.Circle((0,0), Z1_abs, color ='r', fill = False)
axes.add_artist(circle1)
circle2 = plt.Circle((0,0), Z2_abs, color ='m', fill = False)
axes.add_artist(circle2)
circle3 = plt.Circle((0,0), Z3_abs, color ='b', fill = False)
axes.add_artist(circle3)

# Plot orthogonals to add the directional restraint
plt.plot([-X3, X3], [R3,-R3], 'k--', linewidth=3)

# Set plot title
plt.title('Relay characteristics')
plt.legend()

# Set aspect ratio to be equal
plt.gca().set_aspect('equal', adjustable='box')

# Show plot
plt.show()

#%%
# mho Relay
# new center
# Create a new figure
fig2, axes2 = plt.subplots(dpi = 400)

# Plot the x-axis
plt.plot([-10, 30], [0, 0], 'k-', linewidth=1)

# Plot the y-axis
plt.plot([0, 0], [-10, 30], 'k-', linewidth=1)


# Set x and y limits
plt.xlim(-10, 30)
plt.ylim(-10, 30)

# Set x and y labels
plt.xlabel('Re(Z)')
plt.ylabel('Im(Z)')

# Plot vectors
plt.quiver(0, 0, R3, X3, angles='xy', scale_units='xy', scale=1, color='b', label='Z_r3')

CenterR = R3/2
CenterX = X3/2
rad_mho = np.absolute(CenterR+CenterX*1j)

circle_mho = plt.Circle((CenterR,CenterX), rad_mho, color ='k', fill = False)
center = plt.Circle((CenterR, CenterX), 0.3,color = 'b', fill = True)
axes2.add_artist(circle_mho)
axes2.add_artist(center)
# Set aspect ratio to be equal
plt.gca().set_aspect('equal', adjustable='box')


#calculate the impedance vector for the relay with 1800A emergency loading
pf =0.95
Im = 1800 * np.sqrt(3)
Theta = np.arccos(pf)
Voltage = 500000
I = Im * np.exp(-1j*Theta) 
Z_emergency = 500000/I
Z_emergency_relay = CTR/VTR* Z_emergency
R4 = np.real(Z_emergency_relay)
X4 = np.imag(Z_emergency_relay)

plt.quiver(0,0, R4, X4, angles='xy', scale_units='xy', scale=1, color='g', label='Emergency load')
plt.legend()
plt.show()

#%%

# Create a new figure
fig, axes = plt.subplots(dpi = 400)

# Plot the x-axis
plt.plot([-20, 20], [0, 0], 'k-', linewidth=1)

# Plot the y-axis
plt.plot([0, 0], [-20, 20], 'k-', linewidth=1)


# Set x and y limits
plt.xlim(-20, 20)
plt.ylim(-20, 20)

# Set x and y labels
plt.xlabel('Re(Z)')
plt.ylabel('Im(Z)')

# Grey out area that#s behind the relay
theta=-np.arccos(X3/Z3_abs)*180/np.pi
grey_area = Wedge((0,0), Z3_abs, 180+theta, theta, fc='lightgrey', edgecolor='k')
axes.add_artist(grey_area)
# Plot vectors
plt.quiver(0, 0, R3, X3, angles='xy', scale_units='xy', scale=1, color='b', label='Z_r3')
plt.quiver(0, 0, R2, X2, angles='xy', scale_units='xy', scale=1, color='m', label='Z_r2')
plt.quiver(0, 0, R1, X1, angles='xy', scale_units='xy', scale=1, color='r', label='Z_r1')
plt.quiver(0,0, R4, X4, angles='xy', scale_units='xy', scale=1, color='g', label='Emergency load')

# Plot circles
circle1 = plt.Circle((0,0), Z1_abs, color ='r', fill = False)
axes.add_artist(circle1)
circle2 = plt.Circle((0,0), Z2_abs, color ='m', fill = False)
axes.add_artist(circle2)
circle3 = plt.Circle((0,0), Z3_abs, color ='b', fill = False)
axes.add_artist(circle3)

# Plot orthogonals to add the directional restraint
plt.plot([-X3, X3], [R3,-R3], 'k--', linewidth=3)

# Set plot title
plt.title('Relay characteristics')
plt.legend()

# Set aspect ratio to be equal
plt.gca().set_aspect('equal', adjustable='box')

# Show plot
plt.show()

