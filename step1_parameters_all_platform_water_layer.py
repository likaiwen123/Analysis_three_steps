#!/usr/bin/env python

# import some modules
import math
import os
from string import Template

'''
'boron concentration',
'clad density',
'clad radius',
'coolant density',
'coolant temperature',
'fuel density',
'fuel radius',
'fuel temperature',
'number of the fuel parts',
'oxygen density'
'oxygen radius',
'pin pitch',
'u-235 enrichment',
'water layers',
'calculation parameters',
'energy grid',
'tally'
'''
# global parameters
parameter = ['fuel temperature',
             'coolant temperature',
             'fuel radius',
             'pin pitch',
             'u-235 enrichment',
             'boron concentration',
             'clad radius',
             'oxygen radius',
             'oxygen density',
             'number of the fuel parts',
             'clad density',
             'fuel density',
             'coolant density',
             'water layers',
             'tally',
             'calculation parameters',
             'energy grid',
             ]
# file related
default_file_name_prefix = 'series_'
# for windows
# sign = '\\'
# for linux, also windows
sign = '/'

# geometry related
max_fuel_parts = 10  # maximum supported(because of the format of numbers) 91
pin_length = 1.0  # changed to 1.0cm to be equal to the model of MPACT
# material related
minimum_boron_concentration = 1.0E-6  # ppm
Avogadro_constant = 6.022140857E+23
density_transform = Avogadro_constant / 1.0E+24  # density *density_transform == M_atom_average * number_density
b10_number_percent = 19.9
num_percent_zr = {'Zr-90': 0.5145, 'Zr-91': 0.1122, 'Zr-92': 0.1715,
                  'Zr-94': 0.1738, 'Zr-96': 0.0280}  # total:100.00%
zr_isotope = ['90', '91', '92', '94', '96']
atom_weight = {'H1': 1.007825, 'B10': 10.012937, 'B11': 11.009305,
               'O16': 15.994915, 'Zr90': 89.904704, 'Zr91': 90.905645,
               'Zr92': 91.905040, 'Zr94': 93.906316, 'Zr96': 95.908276,
               'U235': 235.043923, 'U238': 238.050783,}
# temperature related
cross_section_T = ([293.6, '60c'], [300.0, '42c'], [600.0, '81c'],
                   [900.0, '82c'], [3000.1, '65c'])
lwtr_T = ([294.0, '60t'], [300.0, '01t'], [400.0, '61t'],
          [500.0, '03t'], [600.0, '62t'], [800.0, '63t'],
          [1000.0, '64t'])
temperature_l_limit = 50.0
temperature_h_limit = 5000.0
temp_energy_const = 8.617E-11  # Boltzmann Constant (MeV/K)

# water density list
'''
water_density = {278.15: 1.0000, 283.15: 0.9998, 288.15: 0.9992, 293.15: 0.9983, 298.15: 0.9971, 303.15: 0.9957,
                 308.15: 0.9941, 313.15: 0.9923, 318.15: 0.9902, 323.15: 0.9880, 328.15: 0.9860, 333.15: 0.9830,
                 338.15: 0.9800, 343.15: 0.9780, 348.15: 0.9750, 353.15: 0.9720, 358.15: 0.9680, 363.15: 0.9650,
                 368.15: 0.9620, 373.15: 0.9580, 378.15: 0.9540, 383.15: 0.9510, 388.15: 0.9470, 393.15: 0.9430,
                 398.15: 0.9390, 403.15: 0.9350, 408.15: 0.9310, 413.15: 0.9260, 418.15: 0.9220, 423.15: 0.9180,
                 428.15: 0.9120, 433.15: 0.9070, 438.15: 0.9020, 443.15: 0.8970, 448.15: 0.8930, 453.15: 0.8870,
                 458.15: 0.8820, 463.15: 0.8760, 468.15: 0.8700, 473.15: 0.8640, 493.15: 0.8400, 498.15: 0.8340,
                 513.15: 0.8140, 523.15: 0.7990, 533.15: 0.7840, 548.15: 0.7560, 573.15: 0.7140, 598.15: 0.6540,
                 623.15: 0.5750, 633.15: 0.5280
                 }
'''
water_density = {293.6: 1.0, 565.0: 0.743, 600.0: 0.661}


# modified
# all parameters from the input file
class Parameters:
    # just to be more clear
    def __init__(self, fuel_temp=None,
                 coolant_temp=None,
                 fuel_radius=None,
                 pin_pitch=None,
                 u235enrichment=None,
                 boron_concentration=None,
                 calculation_parameters=None,
                 clad_radius=None,
                 oxygen_radius=None,
                 number_of_the_fuel_parts=None,
                 fuel_density=None,
                 clad_density=None,
                 coolant_density=None,
                 oxygen_density=None,
                 energy_grid=None,
                 tally=None,
                 water_layers=[]
                 ):
        self.parameters = {
            "fuel temperature": fuel_temp,
            "coolant temperature": coolant_temp,
            "fuel radius": fuel_radius,
            "pin pitch": pin_pitch,
            "u-235 enrichment": u235enrichment,
            "boron concentration": boron_concentration,
            "clad radius": clad_radius,
            "oxygen radius": oxygen_radius,
            "number of the fuel parts": number_of_the_fuel_parts,
            "calculation parameters": calculation_parameters,
            "fuel density": fuel_density,
            "clad density": clad_density,
            "coolant density": coolant_density,
            "oxygen density": oxygen_density,
            "energy grid": energy_grid,
            "tally": tally,
            "water layers": water_layers
        }
        self.count = [0, 0]
        pass

    # modified
    def create_new_input_files(self, loop):
        self.pre_check()
        self.count[1] = loop
        os.system("mkdir Fold_" + str(loop))
        info = open(('Fold_' + str(loop) + '%sinfo_list.txt') % sign, 'w')
        info.write("The parameters used in each file are listed as follows:\n")
        prompt1 = ' ' * 4 + 'file' + ' ' * 38 + 'geometry' + ' ' * 70 + 'density' + ' ' * 32 + 'temperature(K)' \
                  + ' ' * 12 + 'enrichment\n'
        prompt2 = ' ' * 4 + 'name' + ' ' * 10 + '|' + ' ' * 6 + 'pitch' + ' ' * 7 + 'clad_radius' + ' ' * 5 \
                  + 'O2_radius' + ' ' * 5 + 'fuel_radius' + ' ' * 5 + 'fuel_parts' + ' ' * 3 + '|' + ' ' * 7 \
                  + 'fuel' + ' ' * 10 + 'oxygen' + ' ' * 10 + 'clad' + ' ' * 10 + 'water' + ' ' * 7 + '|' + ' ' * 3 \
                  + 'fuel' + ' ' * 3 + 'coolant' + ' ' * 3 + '|' + ' ' * 3 + 'U-235(%)' + ' ' * 3 + 'boron(ppm)' \
                  + ' ' * 3 + '|\n'
        info.write(prompt1)
        info.write(prompt2)
        info.close()
        for B_C in self.parameters["boron concentration"]:
            for Clad_D in self.parameters["clad density"]:
                for Clad_R in self.parameters["clad radius"]:
                    for Cool_T in self.parameters["coolant temperature"]:
                        for F_D in self.parameters["fuel density"]:
                            for F_R in self.parameters["fuel radius"]:
                                for F_T in self.parameters["fuel temperature"]:
                                    for N_F_P in self.parameters["number of the fuel parts"]:
                                        for O_R in self.parameters["oxygen radius"]:
                                            for O_D in self.parameters["oxygen density"]:
                                                for P_P in (self.parameters["pin pitch"]):
                                                    for U235 in self.parameters["u-235 enrichment"]:
                                                        if check_parameters(B_C, Clad_D, Clad_R,
                                                                            self.parameters["water layers"], Cool_T,
                                                                            F_D, F_R, F_T, N_F_P, O_D, O_R, P_P,
                                                                            U235):
                                                            self.count[0] += 1
                                                            if not self.parameters["coolant density"]:
                                                                cool_D = line_interpolation(Cool_T, water_density)
                                                            else:
                                                                cool_D = self.parameters["coolant density"][
                                                                    self.parameters["coolant temperature"]
                                                                        .index(Cool_T)]
                                                            create_a_file(B_C, Clad_D, Clad_R, cool_D, Cool_T,
                                                                          F_D, F_R, F_T, N_F_P, O_D, O_R, P_P, U235,
                                                                          self.parameters["calculation parameters"],
                                                                          self.parameters["energy grid"],
                                                                          self.parameters["tally"],
                                                                          self.parameters["water layers"],
                                                                          self.count)
        info = open(('Fold_' + str(loop) + '%sinfo_list.txt') % sign, 'a')
        info.write('\n')
        info.writelines(self.parameters['calculation parameters'])
        info.write('\ntally\n')
        info.write('absorption  fission  scatter    n2n      n3n\n')
        info.write('%7s%10s%9s%9s%9s\n' % tuple(self.parameters['tally']))
        info.write('\nenergy grid\n')
        info.writelines(self.parameters['energy grid'])
        info.close()
        print 'Done!'
        pass  # done

    def pre_check(self):
        items = self.parameters.keys()
        for item in items:
            if item == "tally" or "calculation parameters":
                continue
            else:
                # deduplication
                self.parameters[item] = sorted(set(self.parameters[item]), key=self.parameters[item].index)
                # delete negative figures
                index = 0
                while index < len(self.parameters[item]):
                    if self.parameters[item][index] < 0:
                        del self.parameters[item][index]
                        print "[Warning] Values for " + item + "are not all legal, thus some have been deleted"
                    else:
                        index += 1
        if not self.parameters["number of the fuel parts"]:
            print "[Warning] Input problem for the number of fuel parts, the value will be set as 3"
            self.parameters["number of the fuel parts"] = [3]
        if not self.parameters["clad density"]:
            print "[Warning] Input problem for the density of clad(Zr), the value will be set as 6.44389g/cm^-3"
            self.parameters["clad density"] = [6.44389]
        if not self.parameters["fuel density"]:
            print "[Warning] Input problem for the density of fuel(UO2), the value will be set as 10.97g/cm^-3"
            self.parameters["fuel density"] = [10.97]
        if not self.parameters["oxygen density"]:
            print "[Warning] Input problem for the density of the air(oxygen)," \
                  + " the value will be set as 7.13709E-04g/cm^-3"
            self.parameters["oxygen density"] = [7.13709E-04]
        for item in items:
            if (item != "number of the fuel parts") and (not ("density" in item)) and (item != "tally") and (
                        item != "calculation parameters") and (item != "water layers"):
                if not self.parameters[item]:
                    print "[Error] No legal parameters for " + item + " have been got, the process will be ended"
                    exit(1)
        if not self.parameters["coolant density"]:
            print "[Warning] Input problem for the density of coolant(water), line interpolation will be called"
        elif len(self.parameters["coolant density"]) != len(self.parameters["coolant temperature"]):
            print "[Error] The legal input for coolant density and temperature can not be linked, the script " \
                  + "will be ended"
            exit(1)
        else:
            print "[Warning] The links between density and temperature will be done. If the values are not all legal," \
                  + " the links can be worry."
        if len(self.parameters["tally"]) != 5:
            print "[Error] The parameters for tally should be 5 integers, 0 or 1."
            exit(1)
        if self.parameters["water layers"]:
            self.parameters["water layers"].sort()
        pass  # not checked


def main(screen_width=60):
    welcome(screen_width)
    loop = 0
    while True:
        choice = process(loop)
        if choice == 3:
            loop += 1
        elif choice == 0:
            break
    goodbye(screen_width)
    pass


# modifying needed
def welcome(screen_width):
    print
    print '-' * screen_width
    information = "Welcome to the Script for Creating MCNP Input Files"
    blanks = (screen_width - len(information)) / 2
    print ' ' * blanks + information + ' ' * blanks
    print '-' * screen_width


# checked
# main process
def process(loop):
    while True:
        print "\nPlease choose one operation:"
        print "(1) Create a sample parameter file"
        print "(2) Create the tutorial for this script"
        print "(3) Enter an input file and begin creating"
        print "(4) Quit"
        response = raw_input("\nEnter 1~4 here:")
        if response.isdigit():
            choice = int(response)
            if choice in range(1, 5):
                break
            else:
                print "\nPlease enter 1~4"
        else:
            print "\nPlease enter 1~4"
    if choice == 1:
        create_a_sample_file()
        print 'Done!\n'
    elif choice == 2:
        create_tutorial_file()
        print 'Done!\n'
    elif choice == 3:
        para = Parameters()
        loop += 1
        while True:
            in_file_name = raw_input("\nPlease enter the name of the input file:\n")
            in_file_name.strip()
            try:
                temp = open(in_file_name)
                temp.close()
            except IOError:
                print 'File "%s" does not exist!' % in_file_name
            else:
                break
        get_info_from_file(in_file_name, para)
        para.create_new_input_files(loop)
    else:
        if exit_program():
            return 0
    return choice
    pass


# modifying needed
# after the modification, this tutorial should be changed in the future
def create_tutorial_file():
    input_tutorial = open("tutorial.txt", 'w')
    tutorial = ['Tutorial for the script\n',
                '1.Description\n',
                '  A script for creating input files for MCNP.\n',
                '2.Function: Given enough infomation, it can create a series of input files for MCNP\n',
                '3.Case Description\n',
                '  A pin case, see the picture.\n',
                '4.Input file\n',
                '  This script needs an input file to tell it the basic figures of the case.\n',
                '  The input file should follow certain structure, as described below:\n',
                '  **The blank line at the head of the input file will be ignored.\n',
                '  **Using "//" to hide the explanation\n',
                '  **The script provides a sample input file, what you need to do is ',
                'just filling it and give it to the script\n',
                '  **The xsdir file is also needed. Or you can provide the "temp_info" file that contains',
                'the corresponding relationship between temperatures and crosssection signs.',
                '  **The parameter should be listed in the next line to the label\n',
                '    For example:\n',
                '	  Below is a wrong input:\n',
                '	    Fuel Temperature 900K\n',
                '	  You should enter the parameters in the next line:\n',
                '	    Fuel Temperature\n',
                '	    900K 1000K\n',
                '	  Or\n',
                '	    Fuel Temperature\n',
                '		900K\n',
                '		1000K\n',
                '	  You can use a blank or a new line to seperate different parameters.\n',
                '  **A blank line is needed when the parameters for an item end\n',
                '  **The unit of parameters should not enter, and the default units are listed as follows:\n',
                '    temperature       K\n',
                '	length            cm\n',
                '	enrichment        %\n',
                '	concentration     ppm\n',
                '  **The sample input file contains the default calculating parameters, including:\n',
                '      particle number\n',
                '	  inactive generation number\n',
                '	  total generation number\n',
                '	  location of the source\n',
                '	  energy grid\n',
                '	  tally\n',
                '	All input files created by the script contains the same content in this section,\n',
                '	you can modify it if needed.\n',
                '5.Output\n',
                '  The script will create many input files and a file describe the parameters used in each file.\n'
                ]
    input_tutorial.writelines(tutorial)
    input_tutorial.close()
    pass


# checked
def exit_program():
    response = raw_input("\nSure to quit?(y/n)").lower()
    if (response == 'y') or (response == 'yes'):
        return True
    else:
        return False


# checked
def goodbye(screen_width):
    print
    print '-' * screen_width
    info1 = "The program is going to close. Goodbye!"
    info2 = '<Press any key to exit>'
    number1 = (screen_width - len(info1)) / 2
    number2 = (screen_width - len(info2)) / 2
    print ' ' * number1 + info1 + ' ' * number1
    print ' ' * number2 + info2 + ' ' * number2
    print '-' * screen_width
    raw_input()


# modified
# a sample file for parameters inputting
def create_a_sample_file():
    input_sample = open("sample.txt", 'w')
    samples = ['//This is a sample input file, you can fill the parameters into it and submit to the script\n',
               '//The density of water can be calculated if you leave a blank in the domain of coolant density\n',
               'Boron Concentration\n',
               '\n',
               '\n',
               'Clad Radius\n',
               '\n',
               '\n',
               'Coolant Temperature\n',
               '\n',
               '\n',
               'Fuel Radius\n',
               '\n',
               '\n',
               'Fuel Temperature\n',
               '\n',
               '\n',
               'Number of the Fuel Parts\n',
               '3\n',
               '\n',
               'Oxygen Radius\n',
               '\n',
               '\n',
               'Pin Pitch\n',
               '\n',
               '\n',
               'U-235 Enrichment\n',
               '\n',
               '\n',
               'Water layers\n',
               '\n',
               '\n',
               'Calculation parameters\n',
               'mode    n\n',
               'kcode 50000  1.0 100 4000\n',
               'ksrc 0 0 5.00000E-01\n',
               '\n',
               'Tally\n',
               '//absorption fission scatter n2n n3n\n',
               '1 1 1 1 1\n',
               '\n',
               'Clad Density\n',
               '6.44389\n',
               '\n',
               'Fuel Density\n',
               '10.97\n',
               '\n',
               'Coolant Density\n',
               '\n',
               '\n',
               'Oxygen Density\n',
               '7.13709E-04\n',
               '\n',
               'Energy Grid\n',
               'E0    1.0000E-08\n',
               '      3.0000E-08\n',
               '      4.0000E-08\n',
               '      6.0000E-08\n',
               '      8.0000E-08\n',
               '      1.0000E-07\n',
               '      1.5000E-07\n',
               '      2.0000E-07\n',
               '      2.7500E-07\n',
               '      3.5000E-07\n',
               '      5.0000E-07\n',
               '      6.2500E-07\n',
               '      7.5000E-07\n',
               '      9.2500E-07\n',
               '      9.7500E-07\n',
               '      1.0100E-06\n',
               '      1.0800E-06\n',
               '      1.1300E-06\n',
               '      1.1750E-06\n',
               '      1.2500E-06\n',
               '      1.4500E-06\n',
               '      1.8600E-06\n',
               '      2.4700E-06\n',
               '      3.7300E-06\n',
               '      4.7000E-06\n',
               '      5.0000E-06\n',
               '      5.4000E-06\n',
               '      6.2500E-06\n',
               '      7.1500E-06\n',
               '      8.1000E-06\n',
               '      1.1900E-05\n',
               '      1.4400E-05\n',
               '      3.0000E-05\n',
               '      4.8300E-05\n',
               '      7.6000E-05\n',
               '      1.4300E-04\n',
               '      3.0500E-04\n',
               '      9.5000E-04\n',
               '      2.2500E-03\n',
               '      9.5000E-03\n',
               '      2.0000E-02\n',
               '      5.0000E-02\n',
               '      7.3000E-02\n',
               '      2.0000E-01\n',
               '      4.9200E-01\n',
               '      8.2000E-01\n',
               '      1.3560E+00\n',
               '      2.3540E+00\n',
               '      4.3040E+00\n',
               '      6.4340E+00\n',
               '      2.0000E+01\n',
               ]
    input_sample.writelines(samples)
    input_sample.close()


# checked
# get all the parameters from the input file, would be better to be a member function of the class Parameter
def get_info_from_file(filename, para):
    input_file = open(filename)
    # ignore the head
    while True:
        title = input_file.readline()
        test = title[:2]
        if not ((title == "\n") or (test == '//')):
            break
    item = title
    while item:
        item = item.strip('\n').lower()
        if item in parameter:
            # each parameter gets its list of values from this line
            para.parameters[item] = get_parameters(input_file, item)
        item = input_file.readline()
        while (item is "\n") or (item[0:2] is "//"):
            item = input_file.readline()
    input_file.close()
    pass


# checked
# each parameter gets its list of values from this function
def get_parameters(input_file, item):
    number = parameter.index(item)
    parameters = []
    while True:
        # get a line that is not beginning with '//'
        content = input_file.readline().strip('\n')
        if content[:2] == '//':
            continue
        elif content:
            # not calculation parameters or energy grid or tally
            if number < (len(parameter) - 2):
                remain = content
                while remain:
                    position = remain.find(' ')
                    if position >= 0:
                        para = remain[:position]
                        remain = remain[position + 1:]
                    else:
                        para = remain
                        remain = []
                    if para is '':
                        continue
                    parameters.append(float(para))
            # calculation parameters and energy grid, just copy part of the content of the input file
            else:
                parameters.append(content + '\n')
        else:
            break
    # choice -- tally, absorption, nufission, scatter, n2n, n3n
    # 0 means no tally for the item, others means that tally for this item is needed
    if number == len(parameter) - 3:
        for ch in range(5):
            if int(parameters[ch]):
                parameters[ch] = True
            else:
                parameters[ch] = False
    return parameters  # a list of values
    pass


# check the rationality of the values for one MCNP input file
def check_parameters(b_c, clad_d, clad_r, w_layers, cool_t, f_d, f_r,
                     f_t, n_f_p, o_d, o_r, p_p, u235):
    # single value check
    if (b_c < 0) or (clad_d <= 0) or (cool_t < 0) or (f_d <= 0) or (f_t < 0) or (n_f_p < 1) or (
                p_p <= 0) or (u235 < 0) or (o_d <= 0):
        return False
    # relationship check
    if (f_r <= 0) or (f_r >= o_r) or (o_r >= clad_r) or ((2 * clad_r) >= p_p):
        return False
    # special check
    # n_f_p max: 10 change can be done at the head of the file
    if n_f_p > max_fuel_parts:
        print "The maximum of fuel parts is now " + str(max_fuel_parts) + ", you need to change it for this case."
        return False
    if w_layers and ((w_layers[0] <= clad_r) or (2 * w_layers[len(w_layers) - 1] >= p_p)):
        return False
    return True
    pass


# not used in this mode
# get database suffix from the dictionary defined at the head of this script
def database_choose(fuel_t, database):
    if fuel_t < database[0][0]:
        return database[0][1]
    else:
        for index in range(1, len(database)):
            if fuel_t < ((database[index - 1][0] + database[index][0]) / 2):
                return database[index - 1][1]
        return database[index][1]
    pass


# lib_type: c d y t p u e m g
# continuous-energy neutron tables : 'c',
# discrete-reaction neutron tables : 'd',
# dosimetry tables : 'y',
# S(alpha,beta) thermal : 't',
# continuous-energy photoatomic : 'p',
# continuous-energy photonuclear tables : 'u',
# continuous-energy electron tables : 'e',
# multigroup neutron tables : 'm',
# multigroup photon tables : 'g',
def library_choose_from_xsdir(z_aid, temp, lib_type='c', xsdir_name='xsdir'):
    libs = [''] * len(temp)
    index_list = range(len(temp))
    temp_info = open("temp_info", 'a')
    temp_info.close()
    temp_info = open("temp_info")
    for index_no in range(len(temp)):
        while True:
            piece = temp_info.readline().strip('\n')
            if not piece:
                break
            if ('%7s%7.1f' % (z_aid, temp[index_no])) in piece:
                libs[index_no] = piece[15:]
                index_list.remove(index_no)
                break
    temp_info.close()
    if not index_list:
        return libs
    try:
        xsdir_file = open(xsdir_name)
    except IOError:
        print '[Error] XS directory file "' + xsdir_name + '" does not exist!'
        exit(1)
    temp_dictionary = [[0.0], ['']]
    temp_dictionary[0] = []
    temp_dictionary[1] = []
    line_text = xsdir_file.readline().strip(' \n')
    while not line_text:
        line_text = xsdir_file.readline().strip(" \n")
    while line_text:
        if line_text == 'directory':
            while True:
                item = xsdir_file.readline().strip(' \n')
                if item:
                    if item[:len(z_aid)] == z_aid:
                        dot = item.find('.')
                        blank = item.find(' ')
                        suffix = item[dot + 1:blank]
                        # if the format of suffix is not '.xxn' one day, this should be changed
                        if not (suffix[2] == lib_type):
                            continue
                        item = item[blank:]
                        temp_mev = ''
                        while True:
                            new_blank = item.find(' ')
                            if new_blank < 0:
                                sub_str = item
                                item = ''
                            else:
                                sub_str = item[:new_blank]
                                item = item[(new_blank + 1):]
                            if sub_str and (not (sub_str == 'ptable') and (not (sub_str == "+"))):
                                temp_mev = sub_str
                            elif sub_str == "+":
                                item = xsdir_file.readline().strip(" \n")
                            if not item:
                                break
                        temp_energy = float(temp_mev)
                        # temperature = temp_energy / temp_energy_const
                        temperature = int(temp_energy * 10.0 / temp_energy_const) / 10.0
                        if (temperature > temperature_l_limit) and (temperature < temperature_h_limit):
                            if temperature in temp_dictionary[0]:
                                continue
                            else:
                                temp_dictionary[0].append(temperature)
                                temp_dictionary[1].append(suffix)
                else:
                    line_text = ''
                    break
        else:
            line_text = xsdir_file.readline().strip(' \n')
    temp_info = open("temp_info", 'a')
    for number in index_list:
        error = []
        for index in range(len(temp_dictionary[0])):
            error.append(math.fabs(temp_dictionary[0][index] - temp[number]))
        temp1 = temperature_h_limit
        min_index = -1
        for index in range(len(error)):
            if temp1 > error[index]:
                temp1 = error[index]
                min_index = index
        suffix = temp_dictionary[1][min_index]
        libs[number] = suffix
        info_piece = '%7s%7.1f %s\n' % (z_aid, temp[number], suffix)
        temp_info.write(info_piece)
    xsdir_file.close()
    temp_info.close()
    return libs


# checked
# create one MCNP input file
def create_a_file(b_c, clad_d, clad_r, cool_d, cool_t, f_d, f_r,
                  f_t, n_f_p, o_d, o_r, p_p, u235, cal, e_g, tally, w_layers, count):
    file_name = default_file_name_prefix + str(count[0]) + '.inp'
    info = open(('Fold_' + str(count[1]) + '%sinfo_list.txt') % sign, 'a')
    n_fuel_parts = int(n_f_p)
    # information written into the info file
    # format:
    # information = '%s18' + '|' + '%15.6E' * 4 + '%8i' + ' ' * 8 + '|' + '%15.6E' * 4 + ' ' * 3 + '|' + '%7i' + '8i' \
    #              + ' ' * 5 + '|' + '%9.2f11.2f' + ' ' * 7 + '|\n'
    information = '%-18s|%15.6E%15.6E%15.6E%15.6E%9i' + ' ' * 8 + \
                  '|%15.6E%15.6E%15.6E%15.6E   |%7.1f%8.1f     |%9.2f%11.2f' + ' ' * 7 + '|\n'
    number_info = (file_name, p_p, clad_r, o_r, f_r, n_fuel_parts, f_d, o_d, clad_d, cool_d, f_t, cool_t, u235, b_c)
    info.write(information % number_info)
    info.close()
    f_t_mev = f_t * temp_energy_const
    cool_t_mev = cool_t * temp_energy_const
    output = open('Fold_' + str(count[1]) + sign + file_name, 'w')
    output.write('c ' + file_name + '\n')
    output.write('c\n')
    # transform density to number density
    uo2_average = (atom_weight['U235'] * u235 + atom_weight['U238'] * (100 - u235)) / 100 + atom_weight['O16'] * 2
    uo2_number_d = f_d * density_transform / uo2_average  # number density of U (UO2)
    o_number_density = o_d * density_transform / atom_weight['O16']  # number density of O (H2O)
    clad_average = 0.0
    for atom in zr_isotope:
        clad_average += atom_weight['Zr' + atom] * num_percent_zr['Zr-' + atom]
    clad_number_density = clad_d * density_transform / clad_average
    coolant_average = atom_weight['H1'] * 2 + atom_weight['O16']
    coolant_o_num_density = cool_d * density_transform / coolant_average

    inner0 = '%-2i' + '%4i' + '%14.6E' + ' ' * 4 + '%4i' * 3 + ' ' * 17 + 'imp:n=%1i  tmp=%8.2E $%4.1fK fuel\n'
    layer_num = len(w_layers)
    # part_number,material_number,number_density,geometry_define,importance,temperature_MeV,temperature_K

    # begin the MCNP input file

    # cell card
    # fuel
    # first cell
    output.write(inner0 % (1, 1, 3 * uo2_number_d, -1, -max_fuel_parts - 3 - layer_num, max_fuel_parts + 4 + layer_num,
                           1, f_t_mev, f_t))
    # other cells containing fuel
    for parts in range(1, n_fuel_parts):
        inner_fuel = '%-2i' + '%4i' + '%14.6E' + '%4i' * 4 + ' ' * 17 + 'imp:n=%1i  tmp=%8.2E $%4.1fK fuel\n'
        paras = (parts + 1, 1, uo2_number_d * 3, parts, -parts - 1,
                 -max_fuel_parts - 3 - layer_num, max_fuel_parts + 4 + layer_num, 1, f_t_mev, f_t)
        output.write(inner_fuel % paras)
    # oxygen
    oxygen = '%-2i' + '%4i' + '%14.6E' + '%4i' * 4 + ' ' * 17 + 'imp:n=%1i  tmp=%8.2E $%4.1fK oxygen\n'
    para_oxygen = (max_fuel_parts + 1, 2, o_number_density, n_fuel_parts, -max_fuel_parts - 1,
                   -max_fuel_parts - 3 - layer_num, max_fuel_parts + 4 + layer_num, 1, cool_t_mev, cool_t)
    output.write(oxygen % para_oxygen)
    # clad
    clad = '%-2i' + '%4i' + '%14.6E' + '%4i' * 4 + ' ' * 17 + 'imp:n=%1i  tmp=%8.2E $%4.1fK clad\n'
    para_clad = (max_fuel_parts + 2, 3, clad_number_density, max_fuel_parts + 1, -max_fuel_parts - 2,
                 -max_fuel_parts - 3 - layer_num, max_fuel_parts + 4 + layer_num, 1, cool_t_mev, cool_t)
    output.write(clad % para_clad)
    # coolant
    coolant = "%-2i" + "%4i" + "%14.6E" + "%4i" * 4 + " " * 17 + 'imp:n=%1i  tmp=%8.2E $%4.1fK coolant water\n'
    for index in range(layer_num):
        para_coolant = (
            max_fuel_parts + 3 + index, 4, coolant_o_num_density * 3, max_fuel_parts + 2 + index,
            -max_fuel_parts - 3 - index, -max_fuel_parts - 3 - layer_num, max_fuel_parts + 4 + layer_num,
            1, cool_t_mev, cool_t
        )
        output.write(coolant % para_coolant)
    coolant = '%-2i' + '%4i' + '%14.6E' + '%4i' * 7 + ' ' * 5 + 'imp:n=%1i  tmp=%8.2E $%4.1fK coolant water\n'
    para_coolant = (
        max_fuel_parts + 3 + layer_num, 4, coolant_o_num_density * 3, max_fuel_parts + 2 + layer_num,
        -max_fuel_parts - 3 - layer_num, max_fuel_parts + 4 + layer_num, -max_fuel_parts - 5 - layer_num,
        max_fuel_parts + 6 + layer_num, -max_fuel_parts - 7 - layer_num, max_fuel_parts + 8 + layer_num,
        1, cool_t_mev, cool_t)
    output.write(coolant % para_coolant)
    # outer space vacuum
    outer = '%-2i' + '%4i' + ' ' * 19 + '%3i:' * 5 + '%3i' + ' ' * 5 + 'imp:n=%1i  tmp=%8.2E $%4.1fK outside\n'
    para_outer = (max_fuel_parts + 4 + layer_num, 0, max_fuel_parts + 3 + layer_num, -max_fuel_parts - 4 - layer_num,
                  max_fuel_parts + 5 + layer_num, -max_fuel_parts - 6 - layer_num, max_fuel_parts + 7 + layer_num,
                  -max_fuel_parts - 8 - layer_num, 0, cool_t_mev, cool_t)
    output.write(outer % para_outer)
    output.write('\n')
    # cell card end

    # face card
    face_cylinder = '%-2i    cz%14.6E\n'
    # cylindrical surfaces
    for index_fuel in range(1, n_fuel_parts + 1):
        face_radius = f_r / math.sqrt(float(n_fuel_parts) / float(index_fuel))
        para_face_fuel = (index_fuel, face_radius)
        output.write(face_cylinder % para_face_fuel)
    other_radius = [o_r, clad_r] + w_layers
    for index_others in range(2 + layer_num):
        para_others = (max_fuel_parts + index_others + 1, other_radius[index_others])
        output.write(face_cylinder % para_others)
    # planes
    box_axis = ['z', 'z', 'x', 'x', 'y', 'y']
    half_pin_pitch = p_p / 2.0
    box_para = [pin_length, 0.0, half_pin_pitch, -half_pin_pitch, half_pin_pitch, -half_pin_pitch]
    for index_box in range(6):
        face_box = '*%-2i   p' + box_axis[index_box] + '%14.6E\n'
        para_face_box = (max_fuel_parts + 3 + layer_num + index_box, box_para[index_box])
        output.write(face_box % para_face_box)
    output.write('\n')

    # material
    # suffix_fuel = database_choose(f_t, cross_section_T)
    # suffix_coolant = database_choose(cool_t, cross_section_T)
    # suffix_lwtr = database_choose(cool_t, lwtr_T)
    h1_suffix = library_choose_from_xsdir('1001', [cool_t])  # H1
    h1_lwtr_suffix = library_choose_from_xsdir('lwtr', [cool_t], 't')  # H1   S(alpha, beta)
    b_suffix = [[], []]
    b_suffix[0] = library_choose_from_xsdir('5010', [cool_t])  # B10
    b_suffix[1] = library_choose_from_xsdir('5011', [cool_t])  # B11
    o16_suffix = library_choose_from_xsdir('8016', [f_t, cool_t])  # O16
    zr_suffix = {}
    for weight in zr_isotope:
        zr_suffix[weight] = library_choose_from_xsdir('400' + weight, [cool_t])
    u235_suffix = library_choose_from_xsdir('92235', [f_t])
    u238_suffix = library_choose_from_xsdir('92238', [f_t])
    # the format of the material card, each line
    material_format = Template('${prefix}%5s.${suffix}%14.6E\n')
    # fuel
    u235_num = uo2_number_d * u235 / 100
    u238_num = uo2_number_d * (100 - u235) / 100
    fuel_o_num = uo2_number_d * 2
    # output.write('m1    92235.' + u235_suffix + '%14.6E\n' % u235_num)
    output.write(material_format.substitute(prefix='m1    ', suffix=u235_suffix[0]) % ('92235', u235_num))
    # output.write(' ' * 6 + '92238.' + suffix_fuel + '%14.6E\n' % u238_num)
    output.write(material_format.substitute(prefix=' ' * 6, suffix=u238_suffix[0]) % ('92238', u238_num))
    # output.write(' ' * 7 + '8016.' + suffix_fuel + '%14.6E\n' % fuel_o_num)
    output.write(material_format.substitute(prefix=' ' * 6, suffix=o16_suffix[0]) % ('8016', fuel_o_num))

    # oxygen
    output.write(material_format.substitute(prefix='m2    ', suffix=o16_suffix[1]) % ('8016', o_number_density))

    # clad
    # number percentage
    zr_num = []
    for index_num in range(5):
        zr_num.append(clad_number_density * num_percent_zr['Zr-' + zr_isotope[index_num]])
    # output.write('m3    40090.' + suffix_coolant + '%14.6E\n' % zr_num[0])
    output.write(material_format.substitute(prefix='m3    ', suffix=zr_suffix['90'][0]) % ('40090', zr_num[0]))
    for index_num in range(1, 5):
        # output.write(' ' * 6 + '400' + zr_isotope[index_num] + '.' + suffix_coolant + '%14.6E\n' % zr_num[index_num])
        output.write(material_format.substitute(prefix=' ' * 6, suffix=zr_suffix[zr_isotope[index_num]][0])
                     % ('400' + zr_isotope[index_num], zr_num[index_num]))

    # water
    # output.write('m4     8016.' + suffix_coolant + '%14.6E\n' % (coolant_o_num_density * 2))
    output.write(material_format.substitute(prefix='m4    ', suffix=o16_suffix[1]) %
                 ('8016', coolant_o_num_density))
    # output.write(' ' * 7 + '1001.' + suffix_coolant + '%14.6E\n' % coolant_o_num_density)
    output.write(material_format.substitute(prefix=' ' * 6, suffix=h1_suffix[0]) % ('1001', coolant_o_num_density * 2))
    # boron
    # 10B content may be as low as 19.1% and as high as 20.3% in natural samples.
    if math.fabs(b_c) > minimum_boron_concentration:
        total_boron_num = cool_d * density_transform * b_c / (
            atom_weight['B10'] * b10_number_percent + atom_weight['B11'] * (100 - b10_number_percent)) / (10 ** 4)
        b_num = [total_boron_num * b10_number_percent / 100, total_boron_num * (100 - b10_number_percent) / 100]
        for counts in range(2):
            # output.write(' ' * 7 + '501' + str(counts) + '.' + suffix_coolant + '%14.6E\n' % b_num[counts])
            output.write(material_format.substitute(prefix=' ' * 6, suffix=(b_suffix[counts])[0])
                         % ('501' + str(counts), b_num[counts]))
    output.write('mt4    lwtr.' + h1_lwtr_suffix[0] + '\n')
    # material finished

    # materials spectially for tally
    tally_m = "%-4s%7s.%3s  1.000000E+00\n"
    isotope_tally = ['92235', '92238', '8016', '8016', '40090', '40091', '40092', '40094', '40096', '8016', '1001']
    temp_tally_isotope = [f_t] * 3 + [cool_t] * 8
    for i in range(11):
        output.write(tally_m % ('m' + str(i + 5), isotope_tally[i],
                                library_choose_from_xsdir(isotope_tally[i], [temp_tally_isotope[i]])[0]))
    output.write('mt15   lwtr.%s\n' % library_choose_from_xsdir('lwtr', [cool_t], 't')[0])
    if b_c > minimum_boron_concentration:
        output.write(tally_m % ('m16', '5010', library_choose_from_xsdir('5010', [cool_t])[0]))
        output.write(tally_m % ('m17', '5011', library_choose_from_xsdir('5011', [cool_t])[0]))
    # material card end

    # calculation parameters and energy grid
    output.writelines(cal)
    # the volume that MCNP can not calculate
    volume_un_calculated = (p_p * p_p - math.pi * clad_r * clad_r) * pin_length
    output.write('vol %iJ%13.6E\n' % (n_fuel_parts + 2, volume_un_calculated))
    output.writelines(e_g)
    # tally
    absorption = '(-2:-6)'
    fission = '(-6 -7)'
    scatter = '2:51:52:53:54:55:56:57:58:59:60:61:62:63:64:65:66:67:68:69:70:71:72:73:74:' + \
              '75:76:77:78:79:80:81:82:83:84:85:86:87:88:89:90:91'
    n2n = '16'
    n3n = '17'
    tally_bin = [absorption, fission, scatter, n2n, n3n]
    bin_template = Template("(${coefficient} %s ${reaction_type}) ")
    bins = ''
    count_bins = 0
    for i in range(5):
        if tally[i]:
            # coefficient is the first valuse in the brackets
            bins += bin_template.substitute(coefficient='1.0', reaction_type=tally_bin[i])
            count_bins += 1
    bins = bins.strip(' ')
    bins += '\n'
    output.write("f14:n")
    for num in range(n_fuel_parts):
        output.writelines(' %i' % (num + 1))
    for num in range(3 + layer_num):
        output.writelines(' %i' % (max_fuel_parts + num + 1))
    output.write('\n')

    output.write("f24:n")
    for num in range(n_fuel_parts):
        output.writelines(' %i' % (num + 1))
    output.write('\n')
    head = ['fm24:n   '] + [' ' * 9] * 3
    for mat_num in range(5, 8):
        fm24_bins = bins % ((mat_num,) * count_bins)
        fm24_bins = head[mat_num - 5] + fm24_bins
        sign_break = ':'
        while True:
            position = fm24_bins.find(sign_break, 65, len(fm24_bins) - 2) + 1
            if position > 0:
                output.write(fm24_bins[:position])
                output.write(' &\n  ')
                fm24_bins = fm24_bins[position:]
                sign_break = ':'
            elif (len(fm24_bins) > 70) and (sign_break == ':'):
                sign_break = ')'
            else:
                break
        output.write(fm24_bins)

    output.write('f34:n %i\n' % (max_fuel_parts + 1))
    head = ['fm34:n   '] + [' ' * 9] * 5
    for mat_num in range(8, 9):
        fm34_bins = bins % ((mat_num,) * count_bins)
        fm34_bins = head[mat_num - 8] + fm34_bins
        sign_break = ':'
        while True:
            position = fm34_bins.find(sign_break, 65, len(fm34_bins) - 2) + 1
            if position > 0:
                output.write(fm34_bins[:position])
                output.write(' &\n  ')
                fm34_bins = fm34_bins[position:]
                sign_break = ':'
            elif (len(fm34_bins) > 70) and (sign_break == ':'):
                sign_break = ')'
            else:
                break
        output.write(fm34_bins)

    output.write('f44:n %i\n' % (max_fuel_parts + 2))
    head = ['fm44:n   '] + [' ' * 9] * 5
    for mat_num in range(9, 14):
        fm44_bins = bins % ((mat_num,) * count_bins)
        fm44_bins = head[mat_num - 9] + fm44_bins
        sign_break = ':'
        while True:
            position = fm44_bins.find(sign_break, 65, len(fm44_bins) - 2) + 1
            if position > 0:
                output.write(fm44_bins[:position])
                output.write(' &\n  ')
                fm44_bins = fm44_bins[position:]
                sign_break = ':'
            elif (len(fm44_bins) > 70) and (sign_break == ':'):
                sign_break = ')'
            else:
                break
        output.write(fm44_bins)

    f54_template = "f54:n" + " %i" * (1 + layer_num) + '\n'
    cell_nums = tuple(range(max_fuel_parts + 3, max_fuel_parts + 4 + layer_num))
    output.write(f54_template % cell_nums)
    addition = 0
    if b_c > minimum_boron_concentration:
        addition = 2
    head = ['fm54:n   '] + [' ' * 9] * (2 + addition)
    for mat_num in range(14, 16 + addition):
        fm54_bins = bins % ((mat_num,) * count_bins)
        fm54_bins = head[mat_num - 14] + fm54_bins
        sign_break = ':'
        while True:
            position = fm54_bins.find(sign_break, 65, len(fm54_bins) - 2) + 1
            if position > 0:
                output.write(fm54_bins[:position])
                output.write(' &\n  ')
                fm54_bins = fm54_bins[position:]
                sign_break = ':'
            elif (len(fm54_bins) > 70) and (sign_break == ':'):
                sign_break = ')'
            else:
                break
        output.write(fm54_bins)
    # file end
    output.close()
    pass


# line interpolation
def line_interpolation(key, dic):
    keys = dic.keys()
    keys.sort()
    number = len(keys)
    result = 0.0
    if key < keys[0]:
        print "Warning: the key is too small"
        result = dic[keys[0]] + (dic[keys[1]] - dic[keys[0]]) * (key - keys[0]) / (keys[1] - keys[0])
    elif key > keys[number - 1]:
        print "Warning: the key is too big"
        result = dic[keys[number - 1]] + (dic[keys[number - 1]] - dic[keys[number - 2]]) * (key - keys[number - 1]) / (
            keys[number - 1] - keys[number - 2])
    else:
        for index in range(1, number):
            if key <= keys[index]:
                result = dic[keys[index]] - (dic[keys[index]] - dic[keys[index - 1]]) * (keys[index] - key) / (
                    keys[index] - keys[index - 1])
                break
    return result


main()
