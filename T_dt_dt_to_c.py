import pickle
from sympy import *
from sympy.printing.ccode import C99CodePrinter
from sympy.printing.codeprinter import Assignment
from iseven import iseven
import sys

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
A = symbols('A')

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

y_sym = [*q_list, *q_list_dt]
str1_list = []
for i in range(len(y_sym)):
    str1_list.append('double ' + str(y_sym[i]) + '=' + f'state_var[{i}];')
str1 = '\n'.join(str1_list)
print('str1 generated')

T_raw = open('T_dt_dt.pkl', 'rb')
T_dt_dt = Matrix(pickle.load(T_raw))

class CMatrixPrinter(C99CodePrinter):
    def _print_ImmutableDenseMatrix(self, expr):
        sub_exprs, simplified = cse(expr)
        lines = []
        for var, sub_expr in sub_exprs:
            lines.append('double ' + self._print(Assignment(var, sub_expr)))
        M = MatrixSymbol('T_dt_dt', len(q_list_dt_dt), len(q_list_dt_dt))
        return '\n'.join(lines) + '\n' + self._print(Assignment(M, simplified[0]))



p = CMatrixPrinter()
str2 = p.doprint(T_dt_dt)
print('str2 generated')

str3 = str1 + '\n' + str2
str0 = '#include <iostream>'+'\n'+'#include <cmath>'+'\n'+'#include "parameters.h"'+'\n'+'\n'+'void T_dt_dt_f(double* state_var, double* T_dt_dt)'+'\n'+'{\n'
str_end = '\n}'

str_final = str0 + str3 + str_end


T_c = open('T_dt_dt_c.cpp', 'w')
T_c.write(str_final)
T_c.close()





