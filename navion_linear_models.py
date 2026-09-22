"""Hai mô hình Navion 4×4 độc lập; đơn vị ft, s, slug và rad."""

from dataclasses import dataclass
from pathlib import Path
import numpy as np


LONGITUDINAL_DATA = dict(CL0=0.41, CD0=0.05, CLalpha=4.44, CDalpha=0.33,
                        Cmalpha=-0.683, Cmalpha_dot=-4.36, Cmq=-9.96,
                        dCD_dMach=0.0, dCm_dMach=0.0)
LATERAL_DATA = dict(CL0=0.41, Cybeta=-0.564, Clbeta=-0.074, Cnbeta=0.0701,
                   Clp=-0.410, Cnp=-0.0575, Clr=0.107, Cnr=-0.125)


@dataclass(frozen=True)
class LongitudinalParameters:
    S: float = 184.0
    b: float = 33.4
    c: float = 5.7
    S_t: float = 43.0
    AR_t: float = 4.0
    W: float = 2750.0
    g: float = 32.2
    rho: float = 0.002378
    u0: float = 176.0
    speed_of_sound: float = 1117.0
    I_y: float = 3000.0
    x_cg_c: float = 0.295
    x_ac_c: float = 0.25
    x_tail_ac: float = 16.0 + 0.295 * 5.7
    gamma0_deg: float = 0.0
    a0_w_per_deg: float = 0.097
    a0_t_per_deg: float = 0.1
    incompressible: bool = True

    fuselage_dx: tuple = (1.5,1.5,1.5,1.5,2.9,2.9,2.9,2.9,2.9)
    fuselage_width: tuple = (3.0,3.4,3.8,4.2,3.8,3.1,2.3,1.5,0.8)
    fuselage_flow_gradient: tuple = (1.2,1.3,1.4,3.2,0.06,0.18,0.31,0.43,0.55)


LONGITUDINAL_REFERENCE = LongitudinalParameters()


def finite_wing_slope(a0_per_deg, aspect_ratio):

    a0 = a0_per_deg * 180.0 / np.pi
    return a0 / (1.0 + a0 / (np.pi * aspect_ratio))


def longitudinal_geometry(p):
    if min(p.S, p.b, p.c, p.S_t, p.AR_t, p.W, p.g, p.rho, p.u0, p.I_y) <= 0:
        raise ValueError("Các kích thước, khối lượng và điều kiện bay phải dương.")
    lt = p.x_tail_ac - p.x_cg_c * p.c
    if lt <= 0:
        raise ValueError("Đuôi ngang phải nằm sau CG trong mô hình này.")
    AR = p.b**2 / p.S
    aw = finite_wing_slope(p.a0_w_per_deg, AR)
    at = finite_wing_slope(p.a0_t_per_deg, p.AR_t)
    if not (len(p.fuselage_dx) == len(p.fuselage_width) == len(p.fuselage_flow_gradient)):
        raise ValueError("Các mảng Multhopp phải có cùng độ dài.")
    integral = sum(dx*w*w*de for dx,w,de in
                   zip(p.fuselage_dx, p.fuselage_width, p.fuselage_flow_gradient))
    Q = 0.5 * p.rho * p.u0**2
    theta = np.radians(p.gamma0_deg)
    return dict(m=p.W/p.g, Q=Q, QS=Q*p.S, QSc=Q*p.S*p.c, AR=AR,
                b_t=np.sqrt(p.AR_t*p.S_t), l_t=lt, V_H=p.S_t*lt/(p.S*p.c),
                a_w=aw, a_t=at, eps_geometry=2*aw/(np.pi*AR),
                fuselage_integral=integral,
                Cmf_geometry=(180/np.pi)*integral/(36.5*p.S*p.c),
                theta0=theta, Mach=p.u0/p.speed_of_sound,
                CL0=p.W*np.cos(theta)/(Q*p.S))


def longitudinal_calibration():


    r = LONGITUDINAL_REFERENCE
    g = longitudinal_geometry(r)
    d = LONGITUDINAL_DATA
    eta = -d['Cmq']/(2*g['a_t']*g['V_H']*g['l_t']/r.c)
    eps = d['Cmalpha_dot']/d['Cmq']
    cmf = d['Cmalpha'] - g['a_w']*(r.x_cg_c-r.x_ac_c) + eta*g['V_H']*g['a_t']*(1-eps)
    cla_raw = g['a_w'] + eta*(r.S_t/r.S)*g['a_t']*(1-eps)

    e = 2*d['CL0']*d['CLalpha']/(np.pi*g['AR']*d['CDalpha'])
    return dict(eta_effective=eta, eps_factor=eps/g['eps_geometry'],
                Cmf_factor=cmf/g['Cmf_geometry'],
                CLalpha_residual=d['CLalpha']-cla_raw,
                e_effective=e, CD_profile=d['CD0']-d['CL0']**2/(np.pi*e*g['AR']))


def longitudinal_aerodynamics(p):
    g = longitudinal_geometry(p)
    k = longitudinal_calibration()
    r = LONGITUDINAL_REFERENCE
    eta = k['eta_effective']
    eps = k['eps_factor']*g['eps_geometry']
    cmf = k['Cmf_factor']*g['Cmf_geometry']
    cla_raw = g['a_w']+eta*(p.S_t/p.S)*g['a_t']*(1-eps)
    cla = cla_raw+k['CLalpha_residual']*(r.S/p.S)
    cl = g['CL0']
    cd = k['CD_profile'] + cl**2/(np.pi*k['e_effective']*g['AR'])
    cda = 2*cl*cla/(np.pi*k['e_effective']*g['AR'])
    cma = g['a_w']*(p.x_cg_c-p.x_ac_c)+cmf-eta*g['V_H']*g['a_t']*(1-eps)
    cmq = -2*g['a_t']*eta*g['V_H']*g['l_t']/p.c
    cmadot = cmq*eps


    tail_slope = eta*(p.S_t/p.S)*g['a_t']*(1-eps)
    h_np = (g['a_w']*p.x_ac_c + tail_slope*p.x_tail_ac/p.c - cmf)/(g['a_w']+tail_slope)
    mach2 = g['Mach']**2
    if mach2 >= 1:
        raise ValueError("Công thức hiệu chỉnh Mach này chỉ dùng dưới âm.")
    clu_raw = mach2/(1-mach2)*cl
    return {**g, **k, 'eps_effective':eps, 'Cmf_effective':cmf,
            'CLalpha_geometry':cla_raw, 'C_Lalpha':cla, 'C_D0':cd, 'C_Dalpha':cda,
            'C_malpha':cma, 'C_mq':cmq, 'C_malpha_dot':cmadot,
            'C_Du':g['Mach']*LONGITUDINAL_DATA['dCD_dMach'],
            'C_mu':g['Mach']*LONGITUDINAL_DATA['dCm_dMach'],
            'C_Lu_compressible':clu_raw, 'C_Lu':0.0 if p.incompressible else clu_raw,
            'x_np_c':h_np, 'SM':h_np-p.x_cg_c,
            'C_malpha_geometry':g['a_w']*(p.x_cg_c-p.x_ac_c)+g['Cmf_geometry']
               -g['V_H']*g['a_t']*(1-g['eps_geometry']),
            'C_mq_geometry':-2*g['a_t']*g['V_H']*g['l_t']/p.c,
            'C_malpha_dot_geometry':-2*g['a_t']*g['V_H']*g['l_t']/p.c*g['eps_geometry']}


def build_longitudinal_matrix(p):
    c = longitudinal_aerodynamics(p)
    Xu=-(c['C_Du']+2*c['C_D0'])*c['QS']/(p.u0*c['m'])
    Xw=-(c['C_Dalpha']-c['CL0'])*c['QS']/(p.u0*c['m'])
    Zu=-(c['C_Lu']+2*c['CL0'])*c['QS']/(p.u0*c['m'])
    Zw=-(c['C_Lalpha']+c['C_D0'])*c['QS']/(p.u0*c['m'])
    Mu=c['C_mu']*c['QSc']/(p.u0*p.I_y)
    Mw=c['C_malpha']*c['QSc']/(p.u0*p.I_y)
    Mwdot=c['C_malpha_dot']*(p.c/(2*p.u0))*c['QSc']/(p.u0*p.I_y)
    Mq=c['C_mq']*(p.c/(2*p.u0))*c['QSc']/p.I_y
    gu=-p.g*np.cos(c['theta0']); gw=-p.g*np.sin(c['theta0'])
    A=np.array([[Xu,Xw,0,gu],[Zu,Zw,p.u0,gw],
                [Mu+Mwdot*Zu,Mw+Mwdot*Zw,Mq+Mwdot*p.u0,Mwdot*gw],
                [0,0,1,0]],dtype=float)
    return A,{**c,'X_u':Xu,'X_w':Xw,'Z_u':Zu,'Z_w':Zw,'M_u':Mu,'M_w':Mw,
              'M_w_dot':Mwdot,'M_q':Mq}


@dataclass(frozen=True)
class LateralParameters:
    S: float = 184.0
    b: float = 33.4
    c: float = 5.7
    S_v: float = 12.5
    W: float = 2750.0
    g: float = 32.2
    rho: float = 0.002378
    u0: float = 176.0
    Ix: float = 1048.0
    Iz: float = 3530.0
    x_cg_c: float = 0.295
    x_fin_ac: float = 16.0 + 0.295 * 5.7
    wing_sweep_deg: float = 2 + 59/60 + 46/3600
    wing_dihedral_deg: float = 7.5
    wing_taper: float = 0.54
    a0_w_per_deg: float = 0.097
    gamma0_deg: float = 0.0
    neglect_Yp_Yr: bool = True


LATERAL_REFERENCE = LateralParameters()


def lateral_geometry(p):
    if min(p.S,p.b,p.c,p.S_v,p.W,p.g,p.rho,p.u0,p.Ix,p.Iz) <= 0:
        raise ValueError("Tham số vật lý ngang-hướng phải dương.")
    lv=p.x_fin_ac-p.x_cg_c*p.c
    if lv<=0:
        raise ValueError("Đuôi đứng phải nằm sau CG.")
    Q=0.5*p.rho*p.u0**2; AR=p.b**2/p.S
    theta=np.radians(p.gamma0_deg)
    return dict(m=p.W/p.g,Q=Q,QS=Q*p.S,QS_b=Q*p.S*p.b,AR=AR,
                l_v=lv,V_v=p.S_v*lv/(p.S*p.b),theta0=theta,
                sweep=np.radians(p.wing_sweep_deg),dihedral=np.radians(p.wing_dihedral_deg),
                a_w=finite_wing_slope(p.a0_w_per_deg,AR),CL0=p.W*np.cos(theta)/(Q*p.S))


def lateral_calibration():


    p=LATERAL_REFERENCE; g=lateral_geometry(p); d=LATERAL_DATA
    Kbeta=-d['Cybeta']*p.S/p.S_v
    Kr=-d['Cnr']/(2*g['V_v']*g['l_v']/p.b)
    Cyr=2*Kr*(p.S_v/p.S)*(g['l_v']/p.b)
    z=(d['Clr']-d['CL0']/4)*p.b/Cyr
    Cnb_wf=d['Cnbeta']-(-d['Cybeta']*g['l_v']/p.b)
    Clb_tail=(z/p.b)*d['Cybeta']
    wing_dihedral_gain=(d['Clbeta']-Clb_tail)/g['dihedral']
    Clp_geometry=-g['a_w']/12*(1+3*p.wing_taper)/(1+p.wing_taper)
    return dict(K_beta=Kbeta,K_r=Kr,z_effective=z,Cnb_wf_reference=Cnb_wf,
                dihedral_gain=wing_dihedral_gain,Clp_factor=d['Clp']/Clp_geometry,
                Cnp_factor=d['Cnp']/(-d['CL0']/8))


def build_lateral_matrix(p):
    g=lateral_geometry(p); k=lateral_calibration(); r=LATERAL_REFERENCE
    Cyb=-k['K_beta']*p.S_v/p.S
    Cnb_wf=k['Cnb_wf_reference']*(r.S*r.b)/(p.S*p.b)
    Cnb=Cnb_wf-Cyb*g['l_v']/p.b
    Cy_p=g['CL0']*(g['AR']+np.cos(g['sweep']))/(g['AR']+4*np.cos(g['sweep']))*np.tan(g['sweep'])
    Cy_r=2*k['K_r']*(p.S_v/p.S)*(g['l_v']/p.b)
    Cnr=-2*k['K_r']*g['V_v']*(g['l_v']/p.b)
    Clr=g['CL0']/4+(k['z_effective']/p.b)*Cy_r
    Clb=k['dihedral_gain']*g['dihedral']+(k['z_effective']/p.b)*Cyb
    Clp_raw=-g['a_w']/12*(1+3*p.wing_taper)/(1+p.wing_taper)
    Clp=k['Clp_factor']*Clp_raw
    Cnp=k['Cnp_factor']*(-g['CL0']/8)
    Yb=g['QS']*Cyb/g['m']; Yp_full=g['QS_b']*Cy_p/(2*g['m']*p.u0)
    Yr_full=g['QS_b']*Cy_r/(2*g['m']*p.u0)
    Yp=0.0 if p.neglect_Yp_Yr else Yp_full
    Yr=0.0 if p.neglect_Yp_Yr else Yr_full
    Lb=g['QS_b']*Clb/p.Ix; Nb=g['QS_b']*Cnb/p.Iz
    Lp=g['QS']*p.b**2*Clp/(2*p.Ix*p.u0)
    Lr=g['QS']*p.b**2*Clr/(2*p.Ix*p.u0)
    Np=g['QS']*p.b**2*Cnp/(2*p.Iz*p.u0)
    Nr=g['QS']*p.b**2*Cnr/(2*p.Iz*p.u0)
    A=np.array([[Yb/p.u0,Yp/p.u0,Yr/p.u0-1,p.g*np.cos(g['theta0'])/p.u0],
                [Lb,Lp,Lr,0],[Nb,Np,Nr,0],[0,1,np.tan(g['theta0']),0]],dtype=float)
    return A,{**g,**k,'Cy_beta':Cyb,'Cn_beta':Cnb,'Cn_beta_wf':Cnb_wf,
              'Cy_p':Cy_p,'Cy_r':Cy_r,'Cl_beta':Clb,'Cl_p':Clp,'Clp_geometry':Clp_raw,
              'Cn_p':Cnp,'Cl_r':Clr,'Cn_r':Cnr,
              'Y_beta':Yb,'Y_p_full':Yp_full,'Y_r_full':Yr_full,'Y_p':Yp,'Y_r':Yr,
              'L_beta':Lb,'L_p':Lp,'L_r':Lr,'N_beta':Nb,'N_p':Np,'N_r':Nr}


DT=0.01
T_END=60.0
X0_LONGITUDINAL=np.array([0.,0.,0.1,0.])
X0_LATERAL=np.array([0.,0.,0.1,0.])
OUTPUT_DIR=Path(__file__).resolve().parent/'results'


def integrate_rk4(A, x0, dt=DT, t_end=T_END):
    if dt <= 0 or t_end <= 0:
        raise ValueError("DT và T_END phải dương.")
    n_steps = int(round(t_end / dt))
    if not np.isclose(n_steps * dt, t_end):
        raise ValueError("T_END phải là bội số dương của DT.")
    if A.shape != (4, 4) or np.shape(x0) != (4,):
        raise ValueError("Mỗi mô hình phải có A 4x4 và x0 gồm 4 phần tử.")
    t = np.arange(n_steps + 1, dtype=float) * dt
    x = np.empty((n_steps + 1, 4), dtype=float)
    x[0] = x0
    for i in range(n_steps):
        k1 = A @ x[i]
        k2 = A @ (x[i] + 0.5 * dt * k1)
        k3 = A @ (x[i] + 0.5 * dt * k2)
        k4 = A @ (x[i] + dt * k3)
        x[i + 1] = x[i] + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    if not np.all(np.isfinite(x)):
        raise FloatingPointError("Đáp ứng chứa NaN/Inf.")
    return t, x


def spectrum_groups(A, model):
    eigs=np.linalg.eigvals(A)
    pairs=sorted([z for z in eigs if z.imag>1e-8],key=abs)
    real=sorted([z for z in eigs if abs(z.imag)<=1e-8],key=lambda z:z.real)
    if model=='longitudinal' and len(pairs)==2:
        return dict(zip(('Phugoid','Short-period'),pairs))
    if model=='lateral_directional' and len(pairs)==1 and len(real)==2:
        return {'Dutch roll':pairs[0],'Roll':real[0],'Spiral':real[1]}
    return {f'root_{i+1}':z for i,z in enumerate(eigs)}


def simulate_models():
    results = {}
    for model, p, builder, x0 in [
        ("longitudinal", LongitudinalParameters(), build_longitudinal_matrix, X0_LONGITUDINAL),
        ("lateral_directional", LateralParameters(), build_lateral_matrix, X0_LATERAL),
    ]:
        A, coefficients = builder(p)
        t, x = integrate_rk4(A, x0)
        results[model] = dict(A=A, parameters=p, coefficients=coefficients,
                              t=t, x=x, eigenvalues=np.linalg.eigvals(A))
    return results


def main():
    results = simulate_models()
    for model, result in results.items():
        print(model)
        print("A =", result["A"], sep="\n")
        print("Eigenvalues:", result["eigenvalues"])
    return results


if __name__ == "__main__":
    main()
