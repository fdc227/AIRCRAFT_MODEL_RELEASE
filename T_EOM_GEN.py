import sys
from sympy import *
from shape_gen import shape_gen
from torsion_shape_gen import torsion_shape_gen
import numpy as np
from pathos.multiprocessing import ProcessingPool as Pool
from multiprocessing.pool import Pool
from sympylist_to_txt import sympylist_to_txt
from iseven import iseven
import pickle

print(sys.argv[1])
print(sys.argv[2])
num_of_elements, num_of_processes = int(sys.argv[1]), int(sys.argv[2])
# num_of_elements = total number of finite element sections on two beams
# num_of_processes = number of processes of CPU processes for the function, multiprocessing is assumed
if not iseven(num_of_elements):
    raise Exception ("Number of finite beam elements must be even")
else:
    np = num_of_elements // 2
    nq = num_of_elements // 2 + 1

#########################################
######### Define Variables ##############
#########################################

t = symbols('t')
x, w, L, theta_0 = symbols('x, w, L, theta_0')
M, m, x, y, z, g, h, E, I, G, J, x_f, c, s, K, A = symbols('M, m, x, y, z, g, h, E, I, G, J, x_f, c, s, K, A')
rho, V, a_w, gamma, M_thetadot, e = symbols('rho, V, a_w, gamma, M_thetadot, e')
beta, P, Q, R = symbols('beta, P, Q, R')
W_x, W_y, W_z = symbols('W_x, W_y, W_z')
P_s, gamma_alpha = symbols('P_s, gamma_alpha')
MOI = symbols('MOI')

######## Define FUnctions ################
theta = Function('theta')(t)
phi = Function('phi')(t)
psi = Function('psi')(t)
X = Function('X')(t)
Y = Function('Y')(t)
Z = Function('Z')(t)
short_var_list = [theta, phi, psi, X, Y, Z]

short_func_to_sym = {}
for func in short_var_list:
    short_func_to_sym[func] = symbols(str(func)[0:-3])

######### Define Symbols ###################

theta_dt = symbols('theta_dt')
phi_dt = symbols('phi_dt')
psi_dt = symbols('psi_dt')
X_dt = symbols('X_dt')
Y_dt = symbols('Y_dt')
Z_dt = symbols('Z_dt')
short_var_list_dt = [theta_dt, phi_dt, psi_dt, X_dt, Y_dt, Z_dt]
short_var_list_dt_raw = [diff(i, t) for i in short_var_list]

theta_dt_dt = symbols('theta_dt_dt')
phi_dt_dt = symbols('phi_dt_dt')
psi_dt_dt = symbols('psi_dt_dt')
X_dt_dt = symbols('X_dt_dt')
Y_dt_dt = symbols('Y_dt_dt')
Z_dt_dt = symbols('Z_dt_dt')
short_var_list_dt_dt = [theta_dt_dt, phi_dt_dt, psi_dt_dt, X_dt_dt, Y_dt_dt, Z_dt_dt]

short_subs_dict = {}
for i in range(len(short_var_list)):
    short_subs_dict[diff(short_var_list[i], t)] = short_var_list_dt[i]
    short_subs_dict[diff(short_var_list_dt_raw[i], t)] = short_var_list_dt_dt[i]

#########################################

shape_func = shape_gen(4)
shape_func = [i.subs({x:y}) for i in shape_func]
S1, S2, S3, S4 = shape_func[0], shape_func[1], shape_func[2], shape_func[3]
S5_tilde, S6_tilde = torsion_shape_gen()
torsion_shape = [S5_tilde, S6_tilde]
S5 = S5_tilde * (x - x_f)
S6 = S6_tilde * (x - x_f)

# print(shape_func, S5, S6)

##########################################
## Injecting q and q_dot into globals() ##
##########################################

var_q_bending = []
for i in range(np, 0, -1):
    globals()[f'p{i}_b'] = Function(f'p{i}_b')(t)
    var_q_bending.append(globals()[f'p{i}_b'])
for i in range(1, nq):
    globals()[f'q{i}_b'] = Function(f'q{i}_b')(t)
    var_q_bending.append(globals()[f'q{i}_b'])

var_q_bending_dot = []
for i in range(np, 0, -1):
    globals()[f'p{i}_b_dot'] = Function(f'p{i}_b_dot')(t)
    var_q_bending_dot.append(globals()[f'p{i}_b_dot'])
for i in range(1, nq):
    globals()[f'q{i}_b_dot'] = Function(f'q{i}_b_dot')(t)
    var_q_bending_dot.append(globals()[f'q{i}_b_dot'])

var_q_torsion = []
for i in range(np, 0, -1):
    globals()[f'p{i}_t'] = Function(f'p{i}_t')(t)
    var_q_torsion.append(globals()[f'p{i}_t'])
for i in range(1, nq):
    globals()[f'q{i}_t'] = Function(f'q{i}_t')(t)
    var_q_torsion.append(globals()[f'q{i}_t'])

var_q_inplane = []
for i in range(np, 0, -1):
    globals()[f'p{i}_i'] = Function(f'p{i}_i')(t)
    var_q_inplane.append(globals()[f'p{i}_i'])
for i in range(1, nq):
    globals()[f'q{i}_i'] = Function(f'q{i}_i')(t)
    var_q_inplane.append(globals()[f'q{i}_i'])

var_q_inplane_dot = []
for i in range(np, 0, -1):
    globals()[f'p{i}_i_dot'] = Function(f'p{i}_i_dot')(t)
    var_q_inplane_dot.append(globals()[f'p{i}_i_dot'])
for i in range(1, nq):
    globals()[f'q{i}_i_dot'] = Function(f'q{i}_i_dot')(t)
    var_q_inplane_dot.append(globals()[f'q{i}_i_dot'])

var_q_list = [*var_q_bending, *var_q_bending_dot, *var_q_torsion, *var_q_inplane, *var_q_inplane_dot]
var_q_list_to_sym = {}
for term in var_q_list:
    var_q_list_to_sym[term] = symbols(str(term)[0:-3])

###########################################

var_q_bending_dt = []
for i in range(np, 0, -1):
    globals()[f'p{i}_b_dt'] = symbols(f'p{i}_b_dt')
    var_q_bending_dt.append(globals()[f'p{i}_b_dt'])
for i in range(1, nq):
    globals()[f'q{i}_b_dt'] = symbols(f'q{i}_b_dt')
    var_q_bending_dt.append(globals()[f'q{i}_b_dt'])

var_q_bending_dot_dt = []
for i in range(np, 0, -1):
    globals()[f'p{i}_b_dot_dt'] = symbols(f'p{i}_b_dot_dt')
    var_q_bending_dot_dt.append(globals()[f'p{i}_b_dot_dt'])
for i in range(1, nq):
    globals()[f'q{i}_b_dot_dt'] = symbols(f'q{i}_b_dot_dt')
    var_q_bending_dot_dt.append(globals()[f'q{i}_b_dot_dt'])

var_q_torsion_dt = []
for i in range(np, 0, -1):
    globals()[f'p{i}_t_dt'] = symbols(f'p{i}_t_dt')
    var_q_torsion_dt.append(globals()[f'p{i}_t_dt'])
for i in range(1, nq):
    globals()[f'q{i}_t_dt'] = symbols(f'q{i}_t_dt')
    var_q_torsion_dt.append(globals()[f'q{i}_t_dt'])

var_q_inplane_dt = []
for i in range(np, 0, -1):
    globals()[f'p{i}_i_dt'] = symbols(f'p{i}_i_dt')
    var_q_inplane_dt.append(globals()[f'p{i}_i_dt'])
for i in range(1, nq):
    globals()[f'q{i}_i_dt'] = symbols(f'q{i}_i_dt')
    var_q_inplane_dt.append(globals()[f'q{i}_i_dt'])

var_q_inplane_dot_dt = []
for i in range(np, 0, -1):
    globals()[f'p{i}_i_dot_dt'] = symbols(f'p{i}_i_dot_dt')
    var_q_inplane_dot_dt.append(globals()[f'p{i}_i_dot_dt'])
for i in range(1, nq):
    globals()[f'q{i}_i_dot_dt'] = symbols(f'q{i}_i_dot_dt')
    var_q_inplane_dot_dt.append(globals()[f'q{i}_i_dot_dt'])

var_q_list_dt = [*var_q_bending_dt, *var_q_bending_dot_dt, *var_q_torsion_dt, *var_q_inplane_dt, *var_q_inplane_dot_dt]
var_q_list_dt_raw = [diff(i, t) for i in var_q_list]

###########################################

var_q_bending_dt_dt = []
for i in range(np, 0, -1):
    globals()[f'p{i}_b_dt_dt'] = symbols(f'p{i}_b_dt_dt')
    var_q_bending_dt_dt.append(globals()[f'p{i}_b_dt_dt'])
for i in range(1, nq):
    globals()[f'q{i}_b_dt_dt'] = symbols(f'q{i}_b_dt_dt')
    var_q_bending_dt_dt.append(globals()[f'q{i}_b_dt_dt'])

var_q_bending_dot_dt_dt = []
for i in range(np, 0, -1):
    globals()[f'p{i}_b_dot_dt_dt'] = symbols(f'p{i}_b_dot_dt_dt')
    var_q_bending_dot_dt_dt.append(globals()[f'p{i}_b_dot_dt_dt'])
for i in range(1, nq):
    globals()[f'q{i}_b_dot_dt_dt'] = symbols(f'q{i}_b_dot_dt_dt')
    var_q_bending_dot_dt_dt.append(globals()[f'q{i}_b_dot_dt_dt'])

var_q_torsion_dt_dt = []
for i in range(np, 0, -1):
    globals()[f'p{i}_t_dt_dt'] = symbols(f'p{i}_t_dt_dt')
    var_q_torsion_dt_dt.append(globals()[f'p{i}_t_dt_dt'])
for i in range(1, nq):
    globals()[f'q{i}_t_dt_dt'] = symbols(f'q{i}_t_dt_dt')
    var_q_torsion_dt_dt.append(globals()[f'q{i}_t_dt_dt'])

var_q_inplane_dt_dt = []
for i in range(np, 0, -1):
    globals()[f'p{i}_i_dt_dt'] = symbols(f'p{i}_i_dt_dt')
    var_q_inplane_dt_dt.append(globals()[f'p{i}_i_dt_dt'])
for i in range(1, nq):
    globals()[f'q{i}_i_dt_dt'] = symbols(f'q{i}_i_dt_dt')
    var_q_inplane_dt_dt.append(globals()[f'q{i}_i_dt_dt'])

var_q_inplane_dot_dt_dt = []
for i in range(np, 0, -1):
    globals()[f'p{i}_i_dot_dt_dt'] = symbols(f'p{i}_i_dot_dt_dt')
    var_q_inplane_dot_dt_dt.append(globals()[f'p{i}_i_dot_dt_dt'])
for i in range(1, nq):
    globals()[f'q{i}_i_dot_dt_dt'] = symbols(f'q{i}_i_dot_dt_dt')
    var_q_inplane_dot_dt_dt.append(globals()[f'q{i}_i_dot_dt_dt'])

var_q_list_dt_dt = [*var_q_bending_dt_dt, *var_q_bending_dot_dt_dt, *var_q_torsion_dt_dt, *var_q_inplane_dt_dt, *var_q_inplane_dot_dt_dt]

q_sub_dict = {}
for i in range(len(var_q_list)):
    q_sub_dict[diff(var_q_list[i], t)] = var_q_list_dt[i]
    q_sub_dict[diff(diff(var_q_list[i], t), t)] = var_q_list_dt_dt[i]


###########################################
###### Combine var with shape_func ########
###########################################

var_q_bending.insert(np, 0)
var_q_bending_dot.insert(np, 0)
var_q_torsion.insert(np, 0)
var_q_inplane.insert(np, 0)
var_q_inplane_dot.insert(np, 0)

bending_shape_func = []
for i in range(num_of_elements):
    output = S1 * var_q_bending[i] + S2 * var_q_bending_dot[i] + S3 * var_q_bending[i+1] + S4 * var_q_bending_dot[i+1]
    bending_shape_func.append(output)
# print(f'bening shape func is {bending_shape_func}')

torsion_shape_func = []
for i in range(num_of_elements):
    output = S5 * var_q_torsion[i] + S6 * var_q_torsion[i+1]
    torsion_shape_func.append(output)
# print(f'torsion_shape_func is {torsion_shape_func}')

torsion_shape_tilde_func = []
for i in range(num_of_elements):
    output = S5_tilde * var_q_torsion[i] + S6_tilde * var_q_torsion[i+1]
    torsion_shape_tilde_func.append(output)

inplane_shape_func = []
for i in range(num_of_elements):
    output = S1 * var_q_inplane[i] + S2 * var_q_inplane_dot[i] + S3 * var_q_inplane[i+1] + S4 * var_q_inplane_dot[i+1]
    inplane_shape_func.append(output)

#############################################

bending_shape_func_dt = []
for i in range(num_of_elements):
    output = S1 * diff(var_q_bending[i], t) + S2 * diff(var_q_bending_dot[i], t) + S3 * diff(var_q_bending[i+1], t) + S4 * diff(var_q_bending_dot[i+1], t)
    bending_shape_func_dt.append(output)

torsion_shape_func_dt = []
for i in range(num_of_elements):
    output = S5 * diff(var_q_torsion[i], t) + S6 * diff(var_q_torsion[i+1], t)
    torsion_shape_func_dt.append(output)

torsion_shape_tilde_func_dt = []
for i in range(num_of_elements):
    output = S5_tilde * diff(var_q_torsion[i], t) + S6_tilde * diff(var_q_torsion[i+1], t)
    torsion_shape_tilde_func_dt.append(output)

inplane_shape_func_dt = []
for i in range(num_of_elements):
    output = S1 * diff(var_q_inplane[i], t) + S2 * diff(var_q_inplane_dot[i], t) + S3 * diff(var_q_inplane[i+1], t) + S4 * diff(var_q_inplane_dot[i+1], t)
    inplane_shape_func_dt.append(output)
# print(torsion_shape_func_dt)
# print(torsion_shape_tilde_func_dt)

##############################################

T_variables = [*short_var_list_dt_raw, *var_q_list_dt_raw]
W_variables = [*short_var_list, *var_q_list]

subs_dict = {**short_subs_dict, **q_sub_dict, **short_func_to_sym, **var_q_list_to_sym}

###############################################

A = Matrix([[cos(theta)*cos(psi), cos(psi)*sin(theta)*sin(phi)-cos(phi)*sin(psi),cos(phi)*cos(psi)*sin(theta)+sin(phi)*sin(psi)],
            [cos(theta)*sin(psi), cos(phi)*cos(psi)+sin(theta)*sin(phi)*sin(psi), -cos(psi)*sin(phi)+cos(phi)*sin(theta)*sin(psi)],
            [-sin(theta), cos(theta)*sin(phi), cos(theta)*cos(phi)]])

A_omega = Matrix([[-sin(theta), 0, 1],
                [cos(theta)*sin(phi), cos(phi), 0],
                [cos(theta)*cos(phi), -sin(phi), 0]])

Omega = A_omega * Matrix([[diff(psi, t)],[diff(theta, t)],[diff(phi, t)]])

eye_mat = Matrix([[1,0,0],[0,1,0],[0,0,1]])

def diff_gen(n):
    count = 0
    T_local_differentiate = []
    for j in T_variables:
        print(f'Now differentiating {count+1}/{len(T_variables)} of {i+1}/{num_of_elements}th T_term')
        diffed = diff(diff(n, j), t)
        T_local_differentiate.append(diffed.xreplace(subs_dict))
        count += 1
    return T_local_differentiate

def diff_others(n):
    T_local_differentiate = []
    for j in T_variables:
        diffed = diff(diff(n, j), t)
        T_local_differentiate.append(diffed.xreplace(subs_dict))
    return T_local_differentiate

def integral_gen(i):
    bs = bending_shape_func[i]
    bs_dt = bending_shape_func_dt[i]
    ts = torsion_shape_func[i]
    ts_dt = torsion_shape_func_dt[i]
    tst = torsion_shape_tilde_func[i]
    tst_dt = torsion_shape_func_dt[i]
    ips = inplane_shape_func[i]
    ips_dt = inplane_shape_func_dt[i]

    P_b = Matrix([[ips],[y + (i - np) * L],[bs]]) # Note that np here deniotes np pieces of wings
    V_M = Matrix([[diff(X, t)],[diff(Y, t)],[diff(Z, t)]])

    P1 = Matrix([[0], [0], [bs_dt]])
    P2 = Matrix([[0],[0],[tst_dt]])
    P3 = Matrix([[ips_dt],[0],[0]])
    P4 = eye_mat * V_M
    P5 = Omega.cross(P_b)

    S1 = P1.dot(P1)
    S2 = 2 * P1.dot(P2)
    S3 = 2 * P1.dot(P3)
    S4 = 2 * P1.dot(P4)
    S5 = 2 * P1.dot(P5)
    S6 = P2.dot(P2)
    S7 = 2 * P2.dot(P3)
    S8 = 2 * P2.dot(P4)
    S9 = 2 * P2.dot(P5)
    S10 = P3.dot(P3)
    S11 = 2 * P3.dot(P4)
    S12 = 2 * P3.dot(P5)
    S13 = P4.dot(P4)
    S14= 2 * P4.dot(P5)
    S15= P5.dot(P5)

    S_list = [S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15]

    S_integrated = []
    for k in range(len(S_list)):
        integrand = Rational(1/2) * m * integrate(S_list[k], (x, 0, c), (y, 0, L))
        print(f'now integrating {k+1}/{len(S_list)} on {i+1}th wing')
        S_integrated.append(integrand)
    
    n = nsimplify(sum(S_integrated))
    T_local_differentiate = diff_gen(n)
    return T_local_differentiate

R = [i for i in range(num_of_elements)]
p = Pool(num_of_processes)
T_terms_differentiated = p.map(integral_gen, R)

T_linear = Rational(1/2) * M * (diff(X, t)**2 + diff(Y, t)**2 + diff(Z, t)**2)
T_angular = Rational(1/2) * MOI * Omega.dot(Omega)

T_linear_differentiated = diff_others(T_linear)
T_angular_differentiated = diff_others(T_angular)

T_terms_differentiated.append(T_linear_differentiated)
T_terms_differentiated.append(T_angular_differentiated)

print('T_terms_differentiated generated')

T_M = Matrix(T_terms_differentiated)
T_m = T_M.T.tolist()

T_final = []
for i in T_m:
    T_final.append(sum(i))

T_raw = open('T_final.pkl', 'wb')
pickle.dump(T_final, T_raw)
