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

def screw_axis(omega : np.ndarray, q : np.ndarray) -> np.ndarray:
    # omega = unit rotation axis and q = any point on that axis 
    # basically an element in SE(3)
    omega = np.asarray(omega, float)
    q = np.asarray(q, float)
    return np.concatenate([-np.cross(omega, q), omega])


class Manipulator :
    def __init__(self, screws : np.ndarray, home : np.ndarray, name = "arm", joint_limits = None, draw_points = None):
        self.screws = np.asarray(screws, float).reshape(-1, 6) # revolute joints
        self.M = np.asarray(home,float) # pose of ee when every joint angle is 0
        if self.M.shape != (4,4):
            raise ValueError("home pose must be a 4x4 matrix (in se3)")
        self.dof = len(self.screws)
        self.name = name
        self.draw_points = draw_points

    @classmethod
    def from_revolute(cls, axes_and_points : np.ndarray, home, **kwargs) -> "Manipulator":
        # builds objects from a list of (omega, q), one per revolution joint
        screws = [screw_axis(w, q) for (w, q) in axes_and_points]
        return cls(screws, home, **kwargs)

    def _theta(self, thetas : np.ndarray) -> np.ndarray:
        # angles for each frame (also a sanity check)
        thetas = np.asarray(thetas, float).reshape(-1)
        if thetas.size != self.dof:
            raise ValueError(f"{self.name}: expected {self.dof} joint angles, " f"got {thetas.size}")
        return thetas

    def __repr__(self):
        return f"Manipulator(name={self.name!r}, dof={self.dof})"


    # ---------------------------------------------------------------
    # forward kinematics (thetas => EE pose)
    # ---------------------------------------------------------------
    def fk(self, thetas : np.ndarray) -> np.ndarray: # forward kinematics
        # PoE product for given angles and home pose (basically a chain of se3.exp calls)
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
    # Jacobians (velocities)
    # ---------------------------------------------------------------
    def space_jacobian(self, thetas : np.ndarray) -> np.ndarray:
        # V = J*theta' is the linear+angular velocity (theta' joint velocity vector)
        # i-th column in J is the i-th screw axis transported through the joints ahead of it 
        # J_i = Ad_{exp([S1]t1) ... exp([S_{i-1}]t_{i-1})} S_i
        thetas = self._theta(thetas)
        J = np.zeros((6, self.dof))
        T = np.eye(4) 
        for i, (S, th) in enumerate(zip(self.screws, thetas)):  # "analytic" construction
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
    def _wrap(self, a : float) -> float:
        return (a + np.pi) % (2 * np.pi) - np.pi

    def _wrist_zyx(self, Rw :np.ndarray) -> list:
            # computes (alpha, beta, gamma) from Rw = R_z(alpa) R_y(beta) R_x(gamma)
            # see matrix R in B.1 (Modern Robotics)
            if abs(Rw[2,0]) < 1- 1e-9: # check that sin(beta) != ±1
                beta = np.arctan2(-Rw[2,0], np.sqrt(Rw[0,0]**2 + Rw[1,0]**2))
                out = []
                for betas in (beta, -beta):
                    out.append([(np.arctan2(Rw[1,0], Rw[0,0]),
                                betas,
                                np.arctan2(Rw[2,1], Rw[2,2]))])
                    return out
            if Rw[2,0] > 0: # gimbal lock 
                return [(0.0, -np.pi/2, -np.arctan2(Rw[0,1], Rw[1,1]))]
            return [(0.0, np.pi/2, np.arctan2(Rw[0,1], Rw[1,1]))]
    
    def _wrist_zyz(self, Rw : np.ndarray) -> list:
        #
        # computes (alpha, beta, gamma) from Rw = R_z(alpha) R_y(beta) R_z(gamma)
        if abs(Rw[2,2]) < 1 - 1e-9: # if cos(beta) != ±1 then sin(beta) != 0     (recall numpy.array)
            beta = np.arctan2(np.sqrt(Rw[0,2]**2 + Rw[1,2]**2), Rw[2,2])
            out = []
            for betas in (beta, -beta):
                c = np.sign(np.sin(beta))
                out.append((np.arctan2(c * Rw[1,2], c * Rw[0,2]),
                          betas,
                          np.arctan2(c * Rw[2,1], -c * Rw[2,0]))) 
                return out
        if  Rw[2,2] > 0 : # gimbal lock (beta = 0)
            return [(0.0, 0.0, np.arctan2(Rw[1,0], Rw[0,0]))]
        return [(0,0, np.pi, np.arctan2(-Rw[1,0], -Rw[0,0]))]

    def puma_6r(self, T : np.ndarray, L1 : float, a2 : float, a3 : float, tool_offset=0.0) -> list: # with no shoulder offset
        # T is ee pose (in se3) from that we got a set of solutions for angles of joints
        # at most 2^3 solutions 
        # with a zyz_wrist (pinch and rotate)
        R, p = se3.split(T)
        R_M = se3.split(self.M)[0]
        p_wc = p - tool_offset * R[:, 2] # position wrist center
        px, py, pz = p_wc
        r_xy = np.hypot(px,py) # horizontal length
        s = pz - L1 # height above shoulder
        D = (r_xy**2 + s**2 - a2**2 - a3**2) / (2*a2*a3)   # law of cosines (for cos theta_3)
        if abs(D) > 1: # target out of reach
            return []
        base_angle = np.arctan2(px,py)
        solutions = []
        # 2 base solutions (left and right)
        for theta_1, r in [(base_angle, r_xy), (self._wrap(base_angle + np.pi), -r_xy)]:
            # 2 elbow solutions (up and down)
            for elbow in (+1.0, -1.0):
                theta_3 = np.arctan2(elbow * np.sqrt(1 - D**2), D)
                theta_2 = np.arctan2(r,s) - np.arctan2(a3 * np.sin(theta_3), a2 + a3 * np.cos(theta_3))
                # PoE equations
                R03 = se3.split(se3.exp(theta_1 * self.screws[0]) @ 
                                se3.exp(theta_2 * self.screws[1]) @
                                se3.exp(theta_3 * self.screws[2]))[0] 
                R36 = R03.T @ R @ R_M.T
                for theta_4, theta_5, theta_6 in self._wrist_zyz(R36) : # 2 wrist solutions (except gimbal lock)
                    solutions.append(np.array([theta_1, theta_2, theta_3, 
                                               theta_4, theta_5, theta_6]))
        return solutions






    
    # numerical inverse kinematics
    def numerical_ik(self, T : np.ndarray) -> list:
        
        pass




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
    





def practice_6dof():
    # test encoding of 6r arm (no offset)
    L1, L2, L3, L4 = 0.40, 0.40, 0.20, 0.10 # lengts
    joints = [
        ([0, 0, 1], [0, 0, 0]),                    # base yaw
        ([0, 1, 0], [0, 0, L1]),                   # shoulder pitch
        ([0, 1, 0], [0, 0, L1 + L2]),              # elbow pitch
        ([0, 0, 1], [0, 0, L1 + L2 + L3]),         # wrist roll
        ([0, 1, 0], [0, 0, L1 + L2 + L3]),         # wrist pitch
        ([0, 0, 1], [0, 0, L1 + L2 + L3]),         # wrist roll
    ]
    M = se3.make(np.eye(3), [0, 0, L1 + L2 + L3 + L4]) # zero-position 
    # skeleton for drawing: base, shoulder, elbow, wrist, tool tip.
    draw_points = [
        (0, [0, 0, 0]),                        # base (fixed)
        (1, [0, 0, L1]),                       # after joint 1
        (2, [0, 0, L1 + L2]),                  # after joint 2
        (3, [0, 0, L1 + L2 + L3]),             # wrist (after joint 3)
        (6, [0, 0, L1 + L2 + L3 + L4]),        # tool tip (after all joints)
    ]
    return Manipulator.from_revolute(joints, M, name="didactic_6dof",
                                     draw_points=draw_points)





if __name__ == "__main__":
    arm = practice_6dof()
    print(arm)
 
    # fk sanity: at the zero config the tool sits at the top of the chain.
    assert np.allclose(arm.ee_position(np.zeros(6)), [0, 0, 1.10], atol=1e-12)
 
    # The evaluation: analytic vs finite-difference Jacobian.
    err = arm.validate()
    print(f"space Jacobian analytic vs finite-difference: max error {err:.2e}")
 
    # Task-space use.
    rng = np.random.default_rng(1)
    a, b = rng.uniform(-1, 1, 6), rng.uniform(-1, 1, 6)
    print(f"EE geodesic distance (length_scale=1.0): "
          f"{arm.geodesic_distance(a, b):.4f}")
    traj = arm.cartesian_trajectory(a, b, n=20)
    print(f"generated a {len(traj)}-pose Cartesian trajectory")
    assert arm.body_jacobian(a).shape == (6, 6)
 
    # Inverse kinematics: recover joint angles that reach random reachable poses.
    rng2 = np.random.default_rng(7)
    reached = 0
    trials = 50
    '''
    for _ in range(trials):
        T_goal = arm.fk(rng2.uniform(-np.pi, np.pi, 6))
        for attempt in range(6):                       # a few random restarts
            th_sol, ok = arm.ik(T_goal, seed=attempt)
            if ok:
                break
        if ok and se3.geodesic_distance(arm.fk(th_sol), T_goal) < 1e-4:
            reached += 1
    print(f"inverse kinematics reached {reached}/{trials} random poses")
    '''
    print("\nkinematics.py: all checks passed "
          "(so3/se3 validated on the 6-DOF chain).")