"""
kinameatics.py : 6-DOF manipulator
    class Manipulator keeps track of screws, whereeach screw is a pose in SE(3)
    Forward Kinematics (from angle of each joint, determines the end effector pose)
    Inverse Kinematics (from end effector pose, determines angle of each joint)
"""


import numpy as np
import so3
import se3
import math

def screw_axis(omega : list, q : list) -> np.asarray:
    # omega = unit rotation axis and q = any point on that axis 
    # like (rho, phi) in SE(3)
    omega = np.asarray(omega, float)
    q = np.asarray(q, float)
    return np.concatenate([-np.cross(omega, q), omega])


class Manipulator :
    def __init__(self, screws : np.ndarray, home : np.ndarray, name = "arm", joint_limits = None, draw_points = None):
        self.screws = np.asarray(screws, float).reshape(-1, 6) 
        self.M = np.asarray(home,float)
        if self.M.shape != (4,4):
            raise ValueError("hom pose must be a 4x4 matrix (in se3)")
        self.dof = len(self.screws)
        self.name = name
        self.draw_points = draw_points

    @classmethod
    def from_revolute(cls, axes_and_points : list, home, **kwargs) -> "Manipulator":
        # builds objects from a list of (omega, q), one per revolution joint
        screws = [screw_axis(w, q) for (w, q) in axes_and_points]
        return cls(screws, home, **kwargs)

    def _theta(self, thetas : np.ndarray) -> np.ndarray:
        # angles for each frame (also a sanity check)
        thetas = np.asarray(thetas, float).reshape(-1)
        if thetas.size != self.dof:
            raise ValueError(f"{self.name}: expected {self.dof} joint angles, " f"got {thetas.size}")
        return thetas

    # ---------------------------------------------------------------
    # kinematics (thetas => EE pose)
    # ---------------------------------------------------------------
    def fk(self, thetas : np.ndarray) -> np.ndarray: # forward kinematics
        # compute pose of end-effector (basically a chain of se3.exp calls)
        """ ADD LENGTHS, OFFSETS AND ANGLED JOINTS"""
        T = np.eye(4)
        for S, th in zip(self.screws, self._theta(thetas)):
            T = T @ se3.exp(th * S)  # forward kinematics equation
        final_eq = T @ self.M
        return final_eq
 
    def ee_position(self, thetas : np.ndarray) -> np.ndarray:
        # end-effector position, p in (R,p) for element of se3
        ee_pos = se3.split(self.fk(thetas))[1]
        return ee_pos

    def skeleton(self, thetas : np.ndarray) -> np.ndarray:
        # nb of screws = nb of vertices in skeleton
        """  check again  """
        thetas = self._theta(thetas)
        G = [np.eye(4)]                
        T = np.eye(4)
        # G[k] = exp(th_1 * S_1) * ... * exp(th_k * S_k)
        for S, th in zip(self.screws, thetas):
            T = T @ se3.exp(th * S)
            G.append(T.copy())
        dp = self.draw_points or [(0, (0, 0, 0)), (self.dof, se3.split(self.M)[1])] 
        pts = []
        # k is pt index + coord and p is deg of freedom + orientation
        for k, p in dp: 
            hom = np.array([p[0], p[1], p[2], 1.0])
            pts.append((G[k] @ hom)[:3])
        return np.array(pts)


    # ---------------------------------------------------------------
    # Jacobians
    # ---------------------------------------------------------------
    def space_jacobian(self, thetas : np.ndarray) -> np.ndarray:
        # V = J*theta' is the linear+angular velocity (theta' joint velocity vector)
        # i-th column in J is the i-th screw axis transported through the joints ahead of it 
        # J_i = Ad_{exp([S1]t1) ... exp([S_{i-1}]t_{i-1})} S_i
        thetas = self._theta(thetas)
        J = np.zeros((6, self.dof)) # Jacobian
        T = np.eye(4) 
        for i, (S, th) in enumerate(zip(self.screws, thetas)):
            J[:, i] = se3.adjoint(T) @ S
            T = T @ se3.exp(th * S)
        return J
 
    def body_jacobian(self, thetas : np.ndarray) -> np.ndarray:
        # Ad_{T^{-1}} * space jacobian 
        T = self.fk(thetas) 
        return se3.adjoint(se3.inverse(T)) @ self.space_jacobian(thetas)

    # (self evaluation)
    def finite_diff_space_jacobian(self, thetas : np.ndarray, eps=1e-6) -> np.ndarray: 
        # second way of computing Jacobian
        # log( Tn * T0^{-1} ) / eps  +  small wiggle
        thetas = self._theta(thetas)
        T = self.fk(thetas)  
        J = np.zeros((6, self.dof))
        for i in range(self.dof):
            dth = np.zeros(self.dof)
            dth[i] = eps
            J[:, i] = se3.log(self.fk(thetas + dth) @ se3.inverse(T)) / eps    # infinitessimal movement 
        return J

    def validate(self, trials=200, seed=0, tol=1e-6) -> float:
        # analytic vs finite difference Jacobian, returns max_error (validating se3.exp, se3.log and se3.adjoint)
        rng = np.random.default_rng(seed)
        max_err = 0.0
        for _ in range(trials):
            th = rng.uniform(-np.pi, np.pi, self.dof)
            err = float(np.max(np.abs(self.space_jacobian(th) - self.finite_diff_space_jacobian(th)))) # comparison
            max_err = max(max_err, err)
        assert max_err < tol, f"{self.name}: Jacobian mismatch {max_err:.2e} > {tol:.0e}"
        return max_err


    # ---------------------------------------------------------------
    # inverse kinematics (EE pose => thetas)
    # ---------------------------------------------------------------
    # analytic inverse kinematics
    def puma6r(self, T : np.ndarray, lengths : np.ndarray, offsets : np.ndarray) -> np.ndarray: # with an offset length[0]
        # recall T is the ee_pose in se3
        # see drawing of 6R PUMA-type arm (alwyas has 4 solutions, except when p_x = p_y = 0)
        # returns each of the 4 solutions
        R, p = se3.split(T) 
        a_1, a_2, a_3 = lengths # in 6RPUMA with offset as in ModernRobotics, a_1 supposed to be 0
        d_1, d_2, d_3 = offsets # in 6RPUMA with offset as in ModernRobotics, d_1 and d_2 supposed to be 0
        p_x, p_y, p_z = p

        # inverse position problem (2D kinematics first)
        D = (np.sum(p**2) - a_2**2 - a_3**2) / (2 * a_2 * a_3) # law of cosines (for cos theta_3)
        theta_3L = np.atan2( math.sqrt(1 - D**2), D )
        theta_3R = np.atan2( -math.sqrt(1 - D**2), D )
        theta_2L = np.atan2(p_z, math.sqrt(p_x**2 + p_y**2 - d_1**2)) - np.atan2(a_3*np.sin(theta_3L), a_2 + a_3*np.cos(theta_3L))
        theta_2R = np.atan2(p_z, math.sqrt(p_x**2 + p_y**2 - d_1**2)) - np.atan2(a_3*np.sin(theta_3R), a_2 + a_3*np.cos(theta_3R))

        if not ((-1e-9 < p_x < 1e-9) and (-1e-9 < p_y < 1e-9)): # otherwise infinitly many solutions for theta_1
            theta_1L = np.atan2(p_y,p_x) - np.atan2(d_1, math.sqrt(p_x**2 + p_y**2 - d_1**2))
            theta_1R = np.pi + np.atan2(p_y,p_x) + np.atan2( -math.sqrt(p_x**2 + p_y**2 - d_1**2), d_1)
            thetas = np.asarray([[theta_1L, theta_2L, theta_3L],
                                 [theta_1L, theta_2R, theta_3R],
                                 [theta_1R, theta_2L, theta_3L],
                                 [theta_1R, theta_2R, theta_3R]])
            for j in range(len(thetas)): # inverse orientation problem
                R03 = self.denavit_hartenberg(thetas[j], offsets, lengths, (np.pi/2, 0, 0)) 
                R36 = R03.T @ R # 
                thetas_next = self.ZYX_euler_angles(R36)
                thetas[j].append([thetas_next])
            return thetas
        
        thetas = np.asarray([[0.0, theta_2L, theta_3L],
                             [0.0, theta_2R, theta_3R],
                             [0.0, theta_2L, theta_3L],
                             [0.0, theta_2R, theta_3R]])
        for j in range(len(thetas)): # inverse orientation problem
            R03 = self.denavit_hartenberg(thetas[j], offsets, lengths, (np.pi/2, 0, 0)) 
            R36 = R03.T @ R
            thetas_next = self.ZYX_euler_angles(R36)
            thetas[j].append([thetas_next])

        return thetas

    def denavit_hartenberg(self, thetas : np.ndarray, offsets : np.ndarray, lengths : np.ndarray, alpha : np.ndarray) -> np.ndarray: 
        # DH parameters algorithm (works for any number of joints)
        T = np.eye(4)
        for i in range(len(thetas)):
            T_next = [[np.cos(thetas[i]), -np.sin(thetas[i])*np.cos(alpha[i]), np.sin(thetas[i])*np.sin(alpha[i]), lengths[i]*np.cos(thetas[i])],
                      [np.sin(thetas[i]), np.cos(thetas[i])*np.cos(alpha[i]), -np.cos(thetas[i])*np.sin(alpha[i]), lengths[i]*np.sin(thetas[i])],
                      [0,                 np.sin(alpha[i]),                    np.cos(alpha[i]),                   offsets[i]],
                      [0,                 0,                                   0,                                  1]]
            T = T @ T_next
        return T
    
    def ZYX_euler_angles(self, R : np.ndarray) -> np.ndarray:
        # algo for finding alpha, beta, gamma knowing R (algo in B.1.1 of Modern Robotics)
        # 
        # for any R in SO(3) such alpha, beta, gamma exists
        if R[2][0] != 1 and R[2][0] != -1:
            beta = np.atan2(-R[2][0], math.sqrt(R[0][0]**2 + R[1][0]**2))
            alpha = np.atan2(R[1][0], R[0][0])
            gamma = np.atan2(R[2][1], R[2][2])
        if R[2][0] == -1: # here a one parameter family of solutions for alpha and gamma producing the same rotation matrix
            beta = np.pi/2
            alpha = 0 # one solution is the following
            gamma = np.atan2(R[0][1], R[1][1])
        if R[2][0] == 1: # similar
            beta = - np.pi/2
            alpha = 0
            gamma = - np.atan2(R[0][1], R[1][1])
        return np.array([alpha, beta, gamma])



    
    # numerical inverse kinematics
    






    # ---------------------------------------------------------------
    # interpolation 
    # ---------------------------------------------------------------
    def geodesic_distance(self, theta_a : float, theta_b : float, length_scale=1.0) -> float:
        # computes geodesic distance (see se3.py)
        geod_dist = se3.geodesic_distance(self.fk(theta_a), self.fk(theta_b), w_rot=length_scale ** 2, w_trans=1.0)
        return geod_dist
    
    def cartesian_trajectory(self, theta_a : float, theta_b : float, n = 20) -> np.ndarray:
        # computes screw interpolate (continuous evolution of pose, see se3)
        Ta , Tb = self.fk(theta_a), self.fk(theta_b)
        traj_pts = [se3.screw_interpolate(Ta, Tb, t) for t in np.linspace(0, 1, n)]
        return traj_pts






 
    def __repr__(self):
        return f"Manipulator(name={self.name!r}, dof={self.dof})"