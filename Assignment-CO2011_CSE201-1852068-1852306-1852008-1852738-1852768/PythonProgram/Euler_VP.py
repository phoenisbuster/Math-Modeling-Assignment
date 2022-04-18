from dataclasses import dataclass
import numpy as np
import csv
import matplotlib.pyplot as plt
import math
from random import seed
from random import random

#? check the variable again
cap_Co2_air = 4.7
cap_Co2_top = 0.4
n_heat_Co2 = 0.057 #! in equation 3,                                        done data in  
U_blow = 0.52#random() #! in equation 3,                                     float type: random from 0.0 -> 1.0
P_blow = 0.0 #0.5 * (10.0 ** 6.0) #! in equation 3,                              done data in
A_flr = 7.8 * (10.0 ** 4.0) #! in equation 3, 4, 5, 14,17                   done data in
U_ext_Co2 = 0.39#random() #! in equation 4,                                  float type: random from 0.0 -> 1.0
phi_ext_Co2 = 4.3 * (10.0 ** 5.0) #! in quation 4,                          done data in
U_pad = 0.6#random() #! in equation 5,                                      float type: random from 0.0 -> 1.0
phi_pad =  0.0 #! in equation 5                                             done data in
Co2_out = 668.0#! in equation 5,9,15                                        done data in (assume value)
Co2_air = 945.005#! in equation 5,9,15                                        data in (float type: random)
U_Th_Scr = 0.4#random() #!in equation 7,13                                  float type: random from 0.0 -> 1.0
K_Th_Scr = 0.25 * (10.0 ** (-3.0)) #! in equation 7                         done data in
T_air = 291.15 #! in equation 7,10,17                                          assume value in ex3 (haven't done yet)
T_top = 295.15 #! in equation 7,10,17                                          assume value in ex3 (haven't done yet)
g = 9.81 #! in equation 7                                                   constan variable
p_air_average_0 = 1.2 #! top cal p_air_averrage in equation 7               done data in
M_air = 28.96#! top cal p_air_averrage in equation 7                        done data in
h_elevation = 1470.0#! top cal p_air_averrage in equation 7                   done data in
R = 8.3145 * (10.0 ** 3)#! top cal p_air_averrage in equation 7 and help equation 18      constan value       
Co2_top = 918.005 #! in equation 6                       data in (float type: random)
dimensionless_ins_scr = 1.0 #! in equation 11                               done data in
c_leak_age = 10.0 ** (-4.0) #! in equation 12                               done data in
v_wind = 2.9 #! in equation 10,12,13,17                                     done data in
C_d = 0.65#!equation to help 10,13                                          done data in
U_side = 2.0#!equation to help 10,13                                        data in (because A_side = 0 so...)
A_side = 0.0#!equation to help 13                                           data in = 0 ??????
C_w = 0.09#!equation to help 10,13,17                                       done data in
U_roof = 0.5#random()#! in equation 10,17                                   float type: random from 0.0 -> 1.0
T_out = 293.15#! in equation 10                                             done data in
T_average_air = 293.15#! in equation 10                                     done data in
A_roof = 14040.0#! in equation 10                                           done data in
h_side_roof = 0.0#! in equation 10                                          done data in (because A_side = 0 so...)
n_side = 0.0 #! in equation 13                                              done data in
n_side_thr = 0.9 #! in equation 13                                          done data in (= n_roof_thr)
U_Vent_forced = 0.6#random() #! in equation 14                              float type: random from 0.0 -> 1.0
the_wind_speed = 0.0 #! in equation 14                                      done data in
h_Vent = 0.97 #! in equation 17                                             done data in
n_roof = 2.0/3.0#!in equation 16                                            done data in
n_roof_thr = 0.9#!in equation 16                                            done data in
M_ch2o = 30.0 * (10.0 ** (-3.0))#! in equation 18                           done data in
h_C_Buf = 1.0 #!To simplify the assignment, hCBuf will always have a value of 1, meaning that CBuf will have no effect on the CO2 fluctuation.  
                           
def MC_blow_air(): #* equation 3
    return (n_heat_Co2 * U_blow * P_blow)/A_flr
#print("MV_blow_air",MC_blow_air())


def MC_ext_air(): #* equation 4
    return (U_ext_Co2 * phi_ext_Co2)/A_flr
#print("MC_ext_air",MC_ext_air())


def MC_pad_air(): #* equation 5
    f_pad = (U_pad * phi_pad)/A_flr
    return f_pad * (Co2_out - Co2_air)
#print("MC_pad_air",MC_pad_air())


#!helper function 6
def p_air():
    temp = (g * M_air * h_elevation)/(293.15 * R)
    return (p_air_average_0 * math.exp(temp))
#print("p_air",p_air())
def p_top():
    temp1 = (g * M_air * (h_elevation + 4.7))/(293.15 * R)
    return (p_air_average_0 * math.exp(temp1))
#print("p_top",p_top())
def p_air_average():
    return (p_top() + p_air())/2

def f_Th_Scr(): #* equation 7
    diff_T = abs(T_air - T_top) ** (2.0/3.0)
    diff_p = abs(p_air() - p_top())
    left_part = (((g * (1.0 - U_Th_Scr)) / (p_air_average() * 2.0)) * diff_p) ** (1.0/2.0)
    return (U_Th_Scr * K_Th_Scr * diff_T) + ((1.0 - U_Th_Scr) * left_part)
#print("f_Th_Scr",f_Th_Scr())
def MC_air_top(): #* equation 6
    return f_Th_Scr() * (Co2_air - Co2_top)
#print("MC_air_top",MC_air_top())

#!helper equation 9
def f_Vent_roof_side(): #* equation 10
    first_part = ((U_roof * U_side * A_side * A_roof) ** 2.0)/(((U_roof**2.0) * (A_roof**2)) + ((U_side**2.0) * (A_side**2)))
    second_part = (2.0 * g * h_side_roof * (T_air-T_out))/T_average_air
    third_part = (((U_roof * A_roof) + (U_side * A_side))/2.0) ** 2.0
    return (C_d * (((first_part * second_part) + (third_part * C_w * (v_wind ** 2.0))) ** (1.0/2.0)))/A_flr
#print("f_Vent_roof_side",f_Vent_roof_side())

def n_ins_scr(): #*equation 11
    return dimensionless_ins_scr * (2.0 - dimensionless_ins_scr)
#print("n_ins_scr",n_ins_scr())

def f_leak_age(): #*equation 12
    return 0.25 * c_leak_age if v_wind < 0.25 else v_wind * c_leak_age
#print("f_leak_age",f_leak_age())

def f_2time_derivative_Vent_side(): #*equation to help 13
    return (C_d * U_side * A_side * v_wind * math.sqrt(C_w))/(2.0 * A_flr)
#print("f_2time_derivative_Vent_side",f_2time_derivative_Vent_side())

def f_Vent_side(): #*equation 13
    if n_side >= n_side_thr :
        return (n_ins_scr() * f_2time_derivative_Vent_side() + 0.5 * f_leak_age())
    else :
        return (n_ins_scr() * (U_Th_Scr * f_2time_derivative_Vent_side() + (1.0 - U_Th_Scr) * f_Vent_roof_side()) + 0.5 * f_leak_age())
#print("f_Vent_side",f_Vent_side())

def f_Vent_forced(): #*equation 14
    return (n_ins_scr() * U_Vent_forced * the_wind_speed)/A_flr
#print("f_Vent_forced",f_Vent_forced())
#!end helper equation 9

def MC_air_out(): #* equation 9
    return ((f_Vent_side() + f_Vent_forced()) * (Co2_air - Co2_out))
#print("MC_air_out",MC_air_out())

#! helper equation 15
def f_2time_derivative_Vent_roof(): #* equation 17
    first = (C_d * U_roof * A_roof) / (2.0 * A_flr)
    second = (g * h_Vent * (T_air - T_out)) / (2.0 * T_average_air)
    return first * ((second + C_w * (v_wind ** 2.0)) ** (1.0/2.0))
#print("f_2time_derivative_Vent_roof",f_2time_derivative_Vent_roof())

def f_Vent_roof(): #* equation 16
    if n_roof >= n_roof_thr:
        return (n_ins_scr() * f_2time_derivative_Vent_roof() + 0.5 * f_leak_age())
    else:
        return (n_ins_scr() * (U_Th_Scr * f_2time_derivative_Vent_roof() + (1.0 - U_Th_Scr) * f_Vent_roof_side() * n_side) + 0.5 * f_leak_age())
#print("f_Vent_roof",f_Vent_roof())
#! end helper equation 15
def MC_top_out(): #*equation 15
    return (f_Vent_roof() * (Co2_top - Co2_out))
#print("MC_top_out",MC_top_out())
#!helper function 18
T_25_k = 298.15
T_can_k = 292.85
H = 22.0 * (10.0 ** 4.0)
E_j = 37.0 * (10.0 ** 3.0)
S = 710.0
L_a_i = 2.0
J_max_25_leaf = 210.0
p_flr = 0.5
K1 = 0.7
K2 = 0.7
p_can = 0.07
t_gh = 0.78
I_Glob = 295.14
n_Glob_PAR = 2.3 
conversion_factor = 0.385
degree_of_curvature = 0.7
n_Co2_air_stom = 0.67
c_Co2_compensation_point = 1.7
def PAR_Gh(): #* equation 9.19
    return t_gh * n_Glob_PAR * I_Glob
#print("PAR_Gh",PAR_Gh())
def PAR_Gh_can(): #* equation 9.18
    return PAR_Gh() * (1.0 - p_can) * (1.0 - math.exp(-K1 * L_a_i))
#print("PAR_Gh_can",PAR_Gh_can())
def PAR_flr_can(): #* equation 9.20
    return p_flr * PAR_Gh() * (1.0 - p_can) * math.exp(-K1 * L_a_i) * (1.0 - math.exp(-K2 * L_a_i))
#print("PAR_flr_can",PAR_flr_can())
def PAR_can(): #*equation 9.17
    return PAR_flr_can() + PAR_Gh_can()
#print("PAR_can",PAR_can())

R1 = 8.3145
def J_Pot(): #* equation 9.15
    first_power =  (E_j * (T_can_k - T_25_k)) / (R1 * T_25_k * T_can_k)
    second_power = (S * T_25_k - H) / (R1 * T_25_k)
    third_power = (S * T_can_k - H) / (R1 * T_can_k) 
    return (L_a_i * J_max_25_leaf * math.exp(first_power) * (1 + math.exp(second_power)))/(1 + math.exp(third_power))
#print("J_Pot",J_Pot())


def J(): #* equation 9.14
    sqrt_side = ((J_Pot() + (conversion_factor * PAR_can())) ** 2.0) - 4.0 * degree_of_curvature * J_Pot() * PAR_can() * conversion_factor
    return (J_Pot() + conversion_factor * PAR_can() - math.sqrt(sqrt_side))/(2.0 * degree_of_curvature)
#print("J",J())


def Co2_stom(): #* equation 9.21
    return n_Co2_air_stom * Co2_air
#print("Co2_stom",Co2_stom())
def Co2_compensation_point(): #* equation 9.22
    return c_Co2_compensation_point * T_can_k
#print("Co2_compensation_point",Co2_compensation_point())
def P():
    return (J() * (Co2_stom() - Co2_compensation_point()))/(4.0 * (Co2_stom() + (2.0 * Co2_compensation_point())))
#print("P()",P())
def R_():
    return (P() * Co2_compensation_point())/Co2_stom()
#print("R_()",R_())
#!end helper function 18
def MC_air_can(): #*equation 18 
    return M_ch2o * h_C_Buf * (P() - R_())
#print("MC_air_can",MC_air_can())





#! ##############################################
h_air = 4.7#! in equation 8.25                                  done data in
h_top = 0.4#! in equation 8.25                                  done data in
s_r_s = -1.0 #!in equation 8.51                                 done data in
R_can = 2.064927554#!in equation 8.51                                  data in
R_can_sp = 5.0#!equation 8.51                                   done data in
c_evap_night_3 = 1.1 * (10.0 ** (-11))#! in equation 8.52       done data in
c_evap_night_4 = 5.2 * (10.0 ** (-6))#! in equation 8.52          done data in
c_evap_1 = 4.3#! in equation 8.50                               done data in
c_evap_2 = 0.54#! in equation 8.50                              done data in
n_mg_ppm = 0.554#! in equation 8.50                             done data in
Vp_can = 2300.0 #! in equation 8.47                                data in
Vp_air = 1700.0 #! in equation 8.47                                data in (unknown)
r_s_min = 82.0#! in equation 8.49                               done data in (maybe a constant)
c_p_air = 10.0 ** 3.0 #! in equation 8.48                       done data in
delta_H = 2.45 * (10.0 ** 6.0)#!in equation 8.48                done data in
gamma = 65.8#! in equation 8.48                                 done data in (a constant)
r_b = 275.0 #! in equation 8.48                                 done data in (maybe a constant)
n_pad = 0.0#!in equation 8.58                                   done data in
x_pad = 5.0#!in equation 8.58                                   data in(no need because n_pad = 0)
x_out = 2.0#!in equation 8.58                                   data in(no need because n_pad = 0)
n_heat_vap = 4.43 * (10.0 ** (-8.0))#!in equation 8.55          done data in
M_water = 18.0#! in equation 8.62                               data in (a constant)
U_frog = 0.6#random() #! in equation 8.64                       float type: random from 0.0 -> 1.0
phi_frog = 2.0#! in equation 8.64                               data in (maybe a constant)
Vp_top = 1650.0#! in equation 8.45                               data in
Vp_out = 1000.0#! in equation 8.45                               data in
c_HEC_in = 1.86#! help equation 8.43                            done data in
T_cov_in = 294.05#! help equation 8.43                             data in
A_cov = 9.0 * (10.0 ** 4.0)#! help equation 8.43                done data in
T_Th_Scr = 294.05#! help equation 8.43                             data in
U_mech_cool = 0.6#random()#!help equation 8.63                  float type
Cop_mech_cool = 2.0#!in equation 8.63                           data in
P_mech_cool = 10000.0#! in equation 8.63                        data in
T_mech_cool = 1.0#!in equation 8.63                             data in
Vp_mach_cool = 1.0#!in equation 8.63                            data in 

def cap_Vp_air(): #* equation 8.25
    return (M_water * h_air)/(R * T_air)
#print("cap_Vp_air",cap_Vp_air())
def cap_Vp_top(): #* equation 8.25
    return (M_water * h_top)/(R * T_top)
#print("cap_Vp_top",cap_Vp_top())
#!helper function 8.47
def S_r_s(): #* equation 8.51 (maybe this equation is redundant)
    return 1.0/(1.0 + math.exp(s_r_s * (R_can - R_can_sp)))
#print("S_r_s",S_r_s())

def c_evap_3(): #* equation 8.52
    return c_evap_night_3 #! explain in detail in report
#print("c_evap_3",c_evap_3())
def c_evap_4(): #* equation 8.52
    return c_evap_night_4 #! read van11 on page 218
#print("c_evap_4",c_evap_4())
def rf_R_can(): #* equation 8.50
    return (R_can + c_evap_1)/(R_can + c_evap_2)
#print("rf_R_can",rf_R_can())
def rf_Co2_air(): #* equation 8.50
    return 1.0 + c_evap_3() * ((n_mg_ppm * Co2_air - 200.0) ** 2.0)
#print("rf_Co2_air",rf_Co2_air())
def rf_Vpcan_Vpair(): #* equation 8.50
    return 1.0 + c_evap_4() * ((Vp_can - Vp_air) ** 2.0)
#print("rf_Vpcan_Vpair",rf_Vpcan_Vpair())

def r_s(): #* equation 8.49
    return r_s_min * rf_R_can() * rf_Co2_air() * rf_Vpcan_Vpair()
#print("r_s",r_s())
def VEC_can_air(): #* equation 8.48
    return (2.0 * p_air() * c_p_air * L_a_i)/(delta_H * gamma * (r_b + r_s()))
#print("VEC_can_air",VEC_can_air())
#!end helper function 8.47
def MV_can_air(): #* equation 8.47
    return (VEC_can_air() * (Vp_can - Vp_air))
#print("MV_can_air",MV_can_air())

def MV_blow_air(): #* equation 8.55
    return (n_heat_vap * U_blow * P_blow)/A_flr
#print("MV_blow_air",MV_blow_air())

def f_pad(): #* to help equation 8.58
    return (U_pad * phi_pad)/A_flr
#print("f_pad",f_pad())
def MV_pad_air(): #* equation 8.58
    return p_air() * f_pad() * (n_pad * (x_pad - x_out) + x_out)
#print("MV_pad_air",MV_pad_air())

def MV_air_out_pad(): #* equation 8.62
    #return (f_pad() * M_water * Vp_air)/(R * (T_air + 273.15))
    return 0
#print("MV_air_out_pad",MV_air_out_pad())

def MV_frog_air(): #* equation 8.64
    #return (U_frog * phi_frog)/A_flr
    return 0
#print("MV_frog_air",MV_frog_air())
#! MV_air_top, MV_air_out and MV_top_out are described analogously to equation 8.45
def MV_air_top(): #* equation 8.45 and according to van11 on page 216, f_air_top = f_Th_Scr 
    return (M_water/R) * f_Th_Scr() * (Vp_air/T_air - Vp_top/T_top)
#print("MV_air_top",MV_air_top())
def MV_air_out(): #* equation 8.45 and according to van11 on page 216, f_air_out = f_vent_side + f_vent_forced
    return (M_water/R) * (f_Vent_side() + f_Vent_forced()) * (Vp_air/T_air - Vp_out/T_out)
#print("MV_air_out",MV_air_out())
def MV_top_out(): #* equation 8.45 and according to van11 on page 216, f_top_out = f_vent_roof
   return (M_water/R) * f_Vent_roof() * (Vp_top/T_top - Vp_out/T_out)
#print("MV_top_out",MV_top_out())

#! MV_air_Th_Scr, MV_top_cov_in and MV_air_mech are described analogously to equation 8.43
def HEC_top_cov_in(): #* to cal MV_top_cov_in
    return (c_HEC_in * A_cov * ((T_top - T_cov_in) ** 0.33))/A_flr #! be careful of negative number power a float number
#print("HEC_top_cov_in",HEC_top_cov_in())

def HEC_air_Th_Scr(): #* to cal MV_air_Th_Scr
    return 1.7 * U_Th_Scr * (abs(T_air - T_Th_Scr) ** 0.33)
#print("HEC_air_Th_Scr",HEC_air_Th_Scr())
 
def HEC_air_mech():#*equation 8.63 to cal MV_air_mech
    numerator = (U_mech_cool * Cop_mech_cool * P_mech_cool)/A_flr
    denominator = T_air - T_mech_cool + 6.4 * (10.0 **(-9.0)) * delta_H * (Vp_air - Vp_mach_cool)
    return numerator/denominator
#print("HEC_air_mech",HEC_air_mech())
Vp_mech = 2400.0
def MV_air_mech(): #* equation 8.43 page 221
    if Vp_air <= Vp_mech :
        return 0
    else :
        return 6.4 * (10.0 ** (-9.0)) * HEC_air_mech() * (Vp_air - Vp_mech)
#print("MV_air_mech",MV_air_mech())
Vp_Th_Scr = 2400.0
def MV_air_Th_Scr(): #* equation 8.43 page 215
    if Vp_air <= Vp_Th_Scr :
        return 0
    else :
        return 6.4 * (10.0 ** (-9.0)) * HEC_air_Th_Scr() * (Vp_air - Vp_Th_Scr) 
Vp_cov_in = 2400.0
#print("MV_air_Th_Scr",MV_air_Th_Scr())
def MV_top_cov_in(): #* equation 8.43 page 215
    if Vp_top <= Vp_cov_in :
        return 0
    else :
        return 6.4 * (10.0 ** (-9.0)) * HEC_top_cov_in() * (Vp_air - Vp_cov_in)
#?print("MV_top_cov_in",MV_top_cov_in())

#! ###############################

def dx_Co2():
    return (MC_blow_air() + MC_ext_air() + MC_pad_air() - MC_air_can() - MC_air_top() - MC_air_out())/cap_Co2_air, (MC_air_top() - MC_top_out())/cap_Co2_top

def dx_Vp():
    return (MV_can_air() + MV_pad_air() + MV_frog_air() + MV_blow_air() - MV_air_Th_Scr() - MV_air_top() - MV_air_out() - MV_air_out_pad() - MV_air_mech())/cap_Vp_air(), (MV_air_top() - MV_top_cov_in() - MV_top_out())/cap_Vp_top()
 
# ! CLASSES
@dataclass
class Euler:
    x_0: float = 0.0
    y_0: float = 0.0
    t: float = 0.0
    h: float = 0.0

    def __init__(self, x0, y0, t0, h):
        self.x_0 = x0
        self.y_0 = y0
        self.t = t0
        self.h = h

    def calculateNewValue(self, x, y):
        x_new = x + self.dx_dt(x, y, self.t) * h
        y_new = y + self.dy_dt(x, y, self.t) * h

        return x_new, y_new

    def calculateNext(self, x, y):
        self.t = self.t + h

        yield self.calculateNewValue(x, y)

    def dx_dt(self, x, y, t):
        rhs = (
            VEC_can_air() * (Vp_can - x) #! MV_can_air
            + MV_pad_air()
            + MV_frog_air()
            + MV_blow_air()
            - 6.4 * (10.0 ** (-9.0)) * HEC_air_Th_Scr() * (x - Vp_Th_Scr) #! MV_air_Th_scr
            - (M_water/R) * f_Th_Scr() * (x/T_air - y/T_top) #! MV_air_top
            - (M_water/R) * (f_Vent_side() + f_Vent_forced()) * (x/T_air - Vp_out/T_out) #! MV_air_out
            - MV_air_out_pad()
            - MV_air_mech()
        )

        return rhs / cap_Vp_air()

    def dy_dt(self, x, y, t):
        rhs = (
            (M_water/R) * f_Th_Scr() * (x/T_air - y/T_top) #! MV_air_top
            - (M_water/R) * f_Vent_roof() * (y/T_top - Vp_out/T_out) #! MV_top_out
            - 6.4 * (10.0 ** (-9.0)) * HEC_top_cov_in() * (x - Vp_cov_in) #! MV_top_cov_in
        )

        return rhs / cap_Vp_top()


print("RUNNING EULER")
x = 1700.0 #! starting point of VPair
y = 1650.0 #! starting point of VPtop
t = 157.0  
h = 1
euler = Euler(x, y, t, h)  # x_0, y_0, t_0, h

list_x = [x]
list_t = [t]
for n in range(5):
    t = t + h
    list_t.append(t)
list_y = [y]
for n in range(5):
    gen_new_val = euler.calculateNext(x, y)
    x, y = next(gen_new_val)

    print(f"step: {n + 1}")
    print(f"new x: {x}\nnew y: {y}\n\n")
    list_x.append(x)
    list_y.append(y)


plt.plot(list_t, list_x, "-b", label = "VP_air")
plt.plot(list_t, list_y, "-r", label = "VP_top")
plt.title("The Vapour pressure below and above the thermal screen from DOY 157-162")
plt.xlabel("time(days)")
plt.ylabel("Vapour pressure")
plt.legend()
plt.show()