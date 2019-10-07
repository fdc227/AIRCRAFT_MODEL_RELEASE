#include <iostream>
#include <mkl.h>
#include <boost/numeric/odeint.hpp>

void strain_terms_f(double* state_var, double* strain_terms);
void T_dt_dt_f(double* state_var, double* T_dt_dt);
void T_others_f(double* state_var, double* T_others);

void ode_gen(const double* Ic, double* dydt, double t);
void write_ode(const double* y, const double t);

int main(void)
{
    double Initial_condition[212];
    Initial_condition[0] = 0;
    Initial_condition[1] =  0;
    Initial_condition[2] =  0;
    Initial_condition[3] =  0;
    Initial_condition[4] =  0;
    Initial_condition[5] =  0;
    Initial_condition[6] =  2.0;
    Initial_condition[7] =  1.8;
    Initial_condition[8] =  1.6;
    Initial_condition[9] =  1.4;
    Initial_condition[10] =  1.2;
    Initial_condition[11] =  1.0;
    Initial_condition[12] =  0.8;
    Initial_condition[13] =  0.6;
    Initial_condition[14] =  0.4;
    Initial_condition[15] =  0.2;
    Initial_condition[16] =  0.2;
    Initial_condition[17] =  0.4;
    Initial_condition[18] =  0.6;
    Initial_condition[19] =  0.8;
    Initial_condition[20] =  1.0;
    Initial_condition[21] =  1.2;
    Initial_condition[22] =  1.4;
    Initial_condition[23] =  1.6;
    Initial_condition[24] =  1.8;
    Initial_condition[25] =  2.0;
    Initial_condition[26] =  -0.2;
    Initial_condition[27] =  -0.2;
    Initial_condition[28] =  -0.2;
    Initial_condition[29] =  -0.2;
    Initial_condition[30] =  -0.2;
    Initial_condition[31] =  -0.2;
    Initial_condition[32] =  -0.2;
    Initial_condition[33] =  -0.2;
    Initial_condition[34] =  -0.2;
    Initial_condition[35] =  -0.2;
    Initial_condition[36] =  0.2;
    Initial_condition[37] =  0.2;
    Initial_condition[38] =  0.2;
    Initial_condition[39] =  0.2;
    Initial_condition[40] =  0.2;
    Initial_condition[41] =  0.2;
    Initial_condition[42] =  0.2;
    Initial_condition[43] =  0.2;
    Initial_condition[44] =  0.2;
    Initial_condition[45] =  0.2;
    Initial_condition[46] =  0;
    Initial_condition[47] =  0;
    Initial_condition[48] =  0;
    Initial_condition[49] =  0;
    Initial_condition[50] =  0;
    Initial_condition[51] =  0;
    Initial_condition[52] =  0;
    Initial_condition[53] =  0;
    Initial_condition[54] =  0;
    Initial_condition[55] =  0;
    Initial_condition[56] =  0;
    Initial_condition[57] =  0;
    Initial_condition[58] =  0;
    Initial_condition[59] =  0;
    Initial_condition[60] =  0;
    Initial_condition[61] =  0;
    Initial_condition[62] =  0;
    Initial_condition[63] =  0;
    Initial_condition[64] =  0;
    Initial_condition[65] =  0;
    Initial_condition[66] =  0;
    Initial_condition[67] =  0;
    Initial_condition[68] =  0;
    Initial_condition[69] =  0;
    Initial_condition[70] =  0;
    Initial_condition[71] =  0;
    Initial_condition[72] =  0;
    Initial_condition[73] =  0;
    Initial_condition[74] =  0;
    Initial_condition[75] =  0;
    Initial_condition[76] =  0;
    Initial_condition[77] =  0;
    Initial_condition[78] =  0;
    Initial_condition[79] =  0;
    Initial_condition[80] =  0;
    Initial_condition[81] =  0;
    Initial_condition[82] =  0;
    Initial_condition[83] =  0;
    Initial_condition[84] =  0;
    Initial_condition[85] =  0;
    Initial_condition[86] =  0;
    Initial_condition[87] =  0;
    Initial_condition[88] =  0;
    Initial_condition[89] =  0;
    Initial_condition[90] =  0;
    Initial_condition[91] =  0;
    Initial_condition[92] =  0;
    Initial_condition[93] =  0;
    Initial_condition[94] =  0;
    Initial_condition[95] =  0;
    Initial_condition[96] =  0;
    Initial_condition[97] =  0;
    Initial_condition[98] =  0;
    Initial_condition[99] =  0;
    Initial_condition[100] =  0;
    Initial_condition[101] =  0;
    Initial_condition[102] =  0;
    Initial_condition[103] =  0;
    Initial_condition[104] =  0;
    Initial_condition[105] =  0;
    Initial_condition[106] =  0;
    Initial_condition[107] =  0;
    Initial_condition[108] =  0;
    Initial_condition[109] =  0;
    Initial_condition[110] =  0;
    Initial_condition[111] =  0;
    Initial_condition[112] =  0;
    Initial_condition[113] =  0;
    Initial_condition[114] =  0;
    Initial_condition[115] =  0;
    Initial_condition[116] =  0;
    Initial_condition[117] =  0;
    Initial_condition[118] =  0;
    Initial_condition[119] =  0;
    Initial_condition[120] =  0;
    Initial_condition[121] =  0;
    Initial_condition[122] =  0;
    Initial_condition[123] =  0;
    Initial_condition[124] =  0;
    Initial_condition[125] =  0;
    Initial_condition[126] =  0;
    Initial_condition[127] =  0;
    Initial_condition[128] =  0;
    Initial_condition[129] =  0;
    Initial_condition[130] =  0;
    Initial_condition[131] =  0;
    Initial_condition[132] =  0;
    Initial_condition[133] =  0;
    Initial_condition[134] =  0;
    Initial_condition[135] =  0;
    Initial_condition[136] =  0;
    Initial_condition[137] =  0;
    Initial_condition[138] =  0;
    Initial_condition[139] =  0;
    Initial_condition[140] =  0;
    Initial_condition[141] =  0;
    Initial_condition[142] =  0;
    Initial_condition[143] =  0;
    Initial_condition[144] =  0;
    Initial_condition[145] =  0;
    Initial_condition[146] =  0;
    Initial_condition[147] =  0;
    Initial_condition[148] =  0;
    Initial_condition[149] =  0;
    Initial_condition[150] =  0;
    Initial_condition[151] =  0;
    Initial_condition[152] =  0;
    Initial_condition[153] =  0;
    Initial_condition[154] =  0;
    Initial_condition[155] =  0;
    Initial_condition[156] =  0;
    Initial_condition[157] =  0;
    Initial_condition[158] =  0;
    Initial_condition[159] =  0;
    Initial_condition[160] =  0;
    Initial_condition[161] =  0;
    Initial_condition[162] =  0;
    Initial_condition[163] =  0;
    Initial_condition[164] =  0;
    Initial_condition[165] =  0;
    Initial_condition[166] =  0;
    Initial_condition[167] =  0;
    Initial_condition[168] =  0;
    Initial_condition[169] =  0;
    Initial_condition[170] =  0;
    Initial_condition[171] =  0;
    Initial_condition[172] =  0;
    Initial_condition[173] =  0;
    Initial_condition[174] =  0;
    Initial_condition[175] =  0;
    Initial_condition[176] =  0;
    Initial_condition[177] =  0;
    Initial_condition[178] =  0;
    Initial_condition[179] =  0;
    Initial_condition[180] =  0;
    Initial_condition[181] =  0;
    Initial_condition[182] =  0;
    Initial_condition[183] =  0;
    Initial_condition[184] =  0;
    Initial_condition[185] =  0;
    Initial_condition[186] =  0;
    Initial_condition[187] =  0;
    Initial_condition[188] =  0;
    Initial_condition[189] =  0;
    Initial_condition[190] =  0;
    Initial_condition[191] =  0;
    Initial_condition[192] =  0;
    Initial_condition[193] =  0;
    Initial_condition[194] =  0;
    Initial_condition[195] =  0;
    Initial_condition[196] =  0;
    Initial_condition[197] =  0;
    Initial_condition[198] =  0;
    Initial_condition[199] =  0;
    Initial_condition[200] =  0;
    Initial_condition[201] =  0;
    Initial_condition[202] =  0;
    Initial_condition[203] =  0;
    Initial_condition[204] =  0;
    Initial_condition[205] =  0;
    Initial_condition[206] =  0;
    Initial_condition[207] =  0;
    Initial_condition[208] =  0;
    Initial_condition[209] =  0;
    Initial_condition[210] =  0;
    Initial_condition[211] =  0;

    boost::numeric::odeint::integrate(ode_gen, Initial_condition, 0.0, 1.0, 0.1, write_ode);

    return 0;
}

void ode_gen(const double& init_condition, double& dxdt, const double t)
{
    double state_var[212];
    for(int i = 0; i < 212; i++)
    {
        state_var[i] = init_condition[i];
    }

    double T_others_out [106];
    double strain_terms_out[106];
    double T_dt_dt_out [11236];

    T_dt_dt_f(state_var, T_dt_dt_out);
    T_others_f(state_var, T_others_out);
    strain_terms_f(state_var, strain_terms_out);

    int N = 106, NRHS=1, LDA=106, LDB=1;
    int ipiv[N];
    int info;
    double* a = T_dt_dt_out;
    double b[106];
    for(int i = 0; i < 106; i++)
    {
        b[i] = -T_others_out[i] - strain_terms_out[i];
    }
    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, NRHS, a, LDA, ipiv, b, LDB );

    for(int i = 0; i < 106; i++)
    {
        dxdt[i] = init_condition[i + 106];
        dxdt[i + 106] = b[i];
    }
}

void write_ode(const double& y, const double t)
{
    for(int i = 1; i < 212; i++)
    {
        if(i < 211)
        {
            std::cout << y[i] << ',';
        }
        else
        {
            std::cout << y[i] << std::endl;
        }
    }
}