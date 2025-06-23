

import sys
sys.path.append('/mnt/data/dg765/cryoCOMP/PEDS/src')

import matplotlib.pyplot as plt

from peds.distributions_fibres import FibreRadiusDistribution, FibreDistribution2d



# Create a fibre radius distribution instance
fibre_dist = FibreRadiusDistribution(
    r_avg=3.5e-3,
    r_min=3e-3,
    r_max=4e-3,
    sigma=0.5e-3,
    gaussian=True,
    seed=42,
)


# Draw 5 sample radii
samples = fibre_dist.draw(5)
print("Sample fibre radii:", samples)

# Create a 2D fibre distribution instance
fibres_2d = FibreDistribution2d(
    n=255,
    domain_size=50e-3,
    r_fibre_dist=fibre_dist,
    volume_fraction=0.45,
    seed=42,
    fast_code=True    #Set to False if you don't have the C++ code or compiler
)

# Generate one fibre distribution alpha matrix from the iterator
alpha = next(iter(fibres_2d))
print("Alpha matrix shape:", alpha.shape)


# Plot the alpha matrix
plt.figure(figsize=(8, 6))
plt.imshow(alpha, cmap='viridis', origin='lower')
plt.colorbar(label='Alpha Values')
plt.title('2D Fibre Distribution Alpha Matrix')
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.show()

