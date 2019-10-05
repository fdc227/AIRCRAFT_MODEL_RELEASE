import sys
from sympy import *
from shape_gen import shape_gen
from torsion_shape_gen import torsion_shape_gen
from cross_product import Cross
from dot_product import Dot
import numpy as np
from pathos.multiprocessing import ProcessingPool as Pool
from sympylist_to_txt import sympylist_to_txt
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

t = symbols('t')
x, w, L, theta_0 = symbols('x, w, L, theta_0')
M, m, x, y, z, g, h, E, I, G, J, x_f, c, s, K = symbols('M, m, x, y, z, g, h, E, I, G, J, x_f, c, s, K')
rho, V, a_w, gamma, M_thetadot, e = symbols('rho, V, a_w, gamma, M_thetadot, e')
beta, P, Q, R = symbols('beta, P, Q, R')
W_x, W_y, W_z = symbols('W_x, W_y, W_z')
P_s, gamma_alpha = symbols('P_s, gamma_alpha')

theta = symbols('theta')
phi = symbols('phi')
psi = symbols('psi')
X = symbols('X')
Y = symbols('Y')
Z = symbols('Z')
short_var_list = [theta, phi, psi, X, Y, Z]

theta_dt = symbols('theta_dt')
phi_dt = symbols('phi_dt')
psi_dt = symbols('psi_dt')
X_dt = symbols('X_dt')
Y_dt = symbols('Y_dt')
Z_dt = symbols('Z_dt')
short_var_list_dt = [theta_dt, phi_dt, psi_dt, X_dt, Y_dt, Z_dt]

theta_dt_dt = symbols('theta_dt_dt')
phi_dt_dt = symbols('phi_dt_dt')
psi_dt_dt = symbols('psi_dt_dt')
X_dt_dt = symbols('X_dt_dt')
Y_dt_dt = symbols('Y_dt_dt')
Z_dt_dt = symbols('Z_dt_dt')
short_var_list_dt_dt = [theta_dt_dt, phi_dt_dt, psi_dt_dt, X_dt_dt, Y_dt_dt, Z_dt_dt]

var_q_bending = []
for i in range(np, 0, -1):
    globals()[f'p{i}_b'] = symbols(f'p{i}_b')
    var_q_bending.append(globals()[f'p{i}_b'])
for i in range(1, nq):
    globals()[f'q{i}_b'] = symbols(f'q{i}_b')
    var_q_bending.append(globals()[f'q{i}_b'])

var_q_bending_dot = []
for i in range(np, 0, -1):
    globals()[f'p{i}_b_dot'] = symbols(f'p{i}_b_dot')
    var_q_bending_dot.append(globals()[f'p{i}_b_dot'])
for i in range(1, nq):
    globals()[f'q{i}_b_dot'] = symbols(f'q{i}_b_dot')
    var_q_bending_dot.append(globals()[f'q{i}_b_dot'])

var_q_torsion = []
for i in range(np, 0, -1):
    globals()[f'p{i}_t'] = symbols(f'p{i}_t')
    var_q_torsion.append(globals()[f'p{i}_t'])
for i in range(1, nq):
    globals()[f'q{i}_t'] = symbols(f'q{i}_t')
    var_q_torsion.append(globals()[f'q{i}_t'])

var_q_inplane = []
for i in range(np, 0, -1):
    globals()[f'p{i}_i'] = symbols(f'p{i}_i')
    var_q_inplane.append(globals()[f'p{i}_i'])
for i in range(1, nq):
    globals()[f'q{i}_i'] = symbols(f'q{i}_i')
    var_q_inplane.append(globals()[f'q{i}_i'])

var_q_inplane_dot = []
for i in range(np, 0, -1):
    globals()[f'p{i}_i_dot'] = symbols(f'p{i}_i_dot')
    var_q_inplane_dot.append(globals()[f'p{i}_i_dot'])
for i in range(1, nq):
    globals()[f'q{i}_i_dot'] = symbols(f'q{i}_i_dot')
    var_q_inplane_dot.append(globals()[f'q{i}_i_dot'])


var_q_list = [*var_q_bending, *var_q_bending_dot, *var_q_torsion, *var_q_inplane, *var_q_inplane_dot]


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


q_list = [*short_var_list, *var_q_list]
q_list_dt = [*short_var_list_dt, *var_q_list_dt]
q_list_dt_dt = [*short_var_list_dt_dt, *var_q_list_dt_dt]

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

torsion_shape_func = []
for i in range(num_of_elements):
    output = S5 * var_q_torsion[i] + S6 * var_q_torsion[i+1]
    torsion_shape_func.append(output)

torsion_shape_tilde_func = []
for i in range(num_of_elements):
    output = S5_tilde * var_q_torsion[i] + S6_tilde * var_q_torsion[i+1]
    torsion_shape_tilde_func.append(output)

inplane_shape_func = []
for i in range(num_of_elements):
    output = S1 * var_q_inplane[i] + S2 * var_q_inplane_dot[i] + S3 * var_q_inplane[i+1] + S4 * var_q_inplane_dot[i+1]
    inplane_shape_func.append(output)


##############################################

W_variables = [*short_var_list, *var_q_list]

###############################################


def strain_terms(i):
    bs = bending_shape_func[i]
    tst = torsion_shape_tilde_func[i]
    ips = inplane_shape_func[i]

    V_bs = integrate(Rational(1/2) * E * A * (diff(diff(bs, y), y)) ** 2, (y, 0, L))
    V_tst = integrate(Rational(1/2) * G * J * (diff(tst, y)) ** 2, (y, 0, L))
    V_ips = integrate(Rational(1/2) * E * A * (diff(diff(ips, y), y)) ** 2, (y, 0, L))

    strain_var_each = []

    for j in W_variables:
        strain_var_each.append(diff(V_bs, j) + diff(V_tst, j) + diff(V_ips, j))
    
    return strain_var_each

R = [i for i in range(6)]
strain_arrays = []
for i in R:
    strain_arrays.append(strain_terms(i))

strain_terms_final = []
for j in range(36):
    local_sum = 0
    for i in range(6):
        local_sum += strain_arrays[i][j]
    strain_terms_final.append(local_sum)

try:
    strain_raw = open('strain_terms.pkl', 'wb')
    pickle.dump(strain_terms_final, strain_raw)
except:
    sympylist_to_txt(strain_terms_final, 'strain_terms.txt')

