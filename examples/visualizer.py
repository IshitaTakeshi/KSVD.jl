from matplotlib import pyplot as plt
from sklearn.datasets import load_digits
import numpy as np

digits = load_digits()

D = np.loadtxt("digits256.dlm").T

fig1 = plt.figure(1)
fig2 = plt.figure(2)

fig1.suptitle("Digits in the dataset", fontsize=18)
for i, image in enumerate(digits.images[:256]):
    ax = fig1.add_subplot(16, 16, i+1)
    ax.imshow(image, cmap='Greys_r')
    ax.set_axis_off()
fig1.savefig("digits256.png", bbox_inches="tight")

fig2.suptitle("Obtained dictionary", fontsize=18)
for i, atom in enumerate(D):
    ax = fig2.add_subplot(16, 16, i+1)
    ax.imshow(atom.reshape(8, 8), cmap='Greys_r')
    ax.set_axis_off()
fig2.savefig("digit_images.png", bbox_inches="tight")

plt.show()
