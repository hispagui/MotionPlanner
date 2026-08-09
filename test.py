import se3
import so3
import numpy as np

T1 = se3.make(so3.exp(np.array([0.3, 1.0, 0.5]) / np.linalg.norm([0.3, 1.0, 0.5]) * 2.4), [3.0, 1.5, 2.0]) # kind of a random pose

print(T1)