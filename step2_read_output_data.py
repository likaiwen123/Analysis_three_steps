#!usr/bin/env python

# import modules
import glob
import numpy as np
import h5py
import sys


# the class Cell, Face, Material, Tally are all members of the class Inputdata, mainly serve as containers of the data
class Cell:
    def __init__(self, content):
        self.cell_id = 0
        self.material_id = 0
        self.number_density = 0.0E+00
        self.isotope_num = 0
        self.raw_geometry = ''
        self.importance = 0
        self.temperature_mev = 0.0E+00
        self.volume = 0.0E+00
        self.other_info = ""
        self.total_flux = []
        self.tally_data = []
        self.get_info_from_line(content)

    def get_info_from_line(self, content):
        count = 0
        coefficient = []
        while count < 3:
            position = content.find(' ')
            if position >= 0:
                para = content[:position]
                content = content[position + 1:]
            else:
                para = content
                content = ""
            if para == '':
                continue
            if count < 2:
                value = int(para)
                coefficient.append(value)
                if value == 0:
                    break
            else:
                value = float(para)
                coefficient.append(value)
            count += 1
        self.cell_id = coefficient[0]
        self.material_id = coefficient[1]
        if not (value == 0):
            self.number_density = coefficient[2]
        content = content.strip(" ")
        length = len(content)
        for index in range(length):
            if content[index:index + 1].isalpha():
                break
        self.raw_geometry = content[:index].strip(" ")
        self.other_info = content[index:]
        content = self.other_info
        imp_pos = content.find("imp")
        imp_pos = content.find("=", imp_pos)
        imp_content = content[imp_pos + 1:]
        while imp_content:
            position = imp_content.find(' ')
            if position >= 0:
                para = imp_content[:position]
                imp_content = imp_content[position + 1:]
            else:
                para = imp_content
                imp_content = []
            if para == '':
                continue
            else:
                self.importance = int(para)
                break
        tmp_pos = content.find("tmp")
        tmp_pos = content.find("=", tmp_pos)
        tmp_content = content[tmp_pos + 1:]
        while tmp_content:
            position = tmp_content.find(' ')
            if position >= 0:
                para = tmp_content[:position]
                tmp_content = tmp_content[position + 1:]
            else:
                para = tmp_content
                tmp_content = []
            if para == '':
                continue
            else:
                self.temperature_mev = float(para)
                break
        pass


class Face:
    def __init__(self, content):
        self.face_id = 0
        self.face_type = ""
        self.parameters = []
        self.reflection = False
        self.get_info_from_line(content)

    def get_info_from_line(self, content):
        parameters = []
        while content:
            position = content.find(' ')
            if position >= 0:
                para = content[:position]
                content = content[position + 1:]
            else:
                para = content
                content = []
            if para is '':
                continue
            parameters.append(para)
        if "*" in parameters[0]:
            self.reflection = True
            parameters[0] = parameters[0].strip("*")
        self.face_id = int(parameters[0])
        self.face_type = parameters[1]
        for value in parameters[2:]:
            self.parameters.append(float(value))
        pass


class Material:
    def __init__(self, content):
        self.mat_id = 0
        self.isotopes = []
        self.sab = []  # maybe no use for this problem
        self.get_info_from_line(content)

    def get_info_from_line(self, content):
        if "mt" in content:
            pos = content.find("mt")
            sab = content[pos:]
            content = content[:pos].strip(" ")
            while sab:
                position = sab.find(' ')
                if position >= 0:
                    para = sab[:position]
                    sab = sab[position + 1:]
                else:
                    para = sab
                    sab = []
                if para is '':
                    continue
                self.sab.append(para)
        parameters = []
        while content:
            position = content.find(' ')
            if position >= 0:
                para = content[:position]
                content = content[position + 1:]
            else:
                para = content
                content = []
            if para is '':
                continue
            parameters.append(para)
        self.mat_id = int(parameters[0][1:])
        number = len(parameters)
        for index in range(1, number, 2):
            isotope = [0, 0.0, ""]
            position = parameters[index].find(".")
            isotope[0] = int(parameters[index][:position])
            isotope[1] = float(parameters[index + 1])
            isotope[2] = parameters[index][position + 1:]
            self.isotopes.append(isotope)
        pass


class Tally:
    def __init__(self, content):
        self.tally_id = 0
        self.type = []
        self.cells = []
        self.multiplier = []
        self.database = []
        self.total_number = 0
        self.get_info_from_line(content)

    def get_info_from_line(self, content):
        if "fm" in content:
            pos = content.find("fm")
            multiplier = content[pos:]
            content = content[:pos].strip(" ")
            start_point = multiplier.find("(")
            multiplier = multiplier[start_point:]
            temporary = list(multiplier)
            count = 0
            bracket_list = [0]
            number = len(temporary)
            for index in range(number):
                item = temporary[index]
                if item == "(":
                    count += 1
                elif item == ")":
                    count -= 1
                    if count == 0:
                        bracket_list.append(index)
            number = len(bracket_list)
            temp = multiplier
            multiplier = []
            for index in range(1, number):
                multiplier.append(temp[bracket_list[index - 1] + 1: bracket_list[index]].strip(" ("))
            for item in multiplier:
                times = 0
                parameters = []
                while times < 2:
                    position = item.find(' ')
                    if position >= 0:
                        para = item[:position]
                        item = item[position + 1:].strip(" ")
                    else:
                        para = item
                        item = []
                    if para is '':
                        continue
                    parameters.append(para)
                    times += 1
                if item[:1] == "(":
                    item = item.strip("() ")
                parameters[0] = float(parameters[0])
                parameters[1] = int(parameters[1])
                parameters.append(item)
                self.multiplier.append(parameters)
        parameters = []
        while content:
            position = content.find(' ')
            if position >= 0:
                para = content[:position]
                content = content[position + 1:]
            else:
                para = content
                content = []
            if para is '':
                continue
            parameters.append(para)
        type_info = parameters[0][1:]
        position = type_info.find(":")
        self.tally_id = int(type_info[:position])
        self.type.append(int(type_info[position - 1:position]))
        self.type.append(type_info[position + 1:])
        number = len(parameters)
        for index in range(1, number):
            self.cells.append(int(parameters[index]))
        self.total_number = len(self.cells) * max(len(self.multiplier), 1)
        pass


class Inputdata:
    def __init__(self, content):
        self.cell_list = []
        self.cell_index = {}
        self.void_cells_list = []
        self.material_list = []
        self.mat_index = {}
        self.isotopes_max_number = 0
        self.face_list = []
        self.ksrc = []
        self.energy_grid = []
        self.tally = []
        self.get_data_from_output(content)
        pass

    def get_data_from_output(self, content):
        cards = self.split_cards(content)
        self.get_cell_data(cards[0])
        self.get_face_data(cards[1])
        self.get_material_data(cards[2])
        self.get_ksrc_data(cards[3])
        self.get_energy_grid_data(cards[4])
        self.get_tally_data(cards[5])
        self.modification()
        pass

    # split the input file part of MCNP output file into several parts(cards), then get data from the corresponding part
    def split_cards(self, content):
        number = len(content)
        cards = [[]] * 6
        board1 = content.index("")
        cards[0] = content[:board1]
        board2 = content.index("", board1 + 1)
        cards[1] = content[board1 + 1:board2]
        for index in range(board2 + 1, number):
            test = list(content[index][:3])
            if test[0].isdigit() or (test[0] == "m" and (test[1].isdigit() or (test[1] == "t"))):
                continue
            else:
                board3 = index
                break
        cards[2] = content[board2 + 1:board3]
        for index in range(board3 + 1, number):
            if content[index][0] == 'e':
                board4 = index
                break
        cards[3] = content[board3:board4]
        for index in range(board4 + 1, number):
            if content[index][0] == 'f':
                board5 = index
                break
        cards[4] = content[board4:board5]
        cards[5] = content[board5:]
        return cards

    # add data to the member object of Class Cell in Class Inputdata
    def get_cell_data(self, content):
        # count = 0
        for index in content:
            self.cell_list.append(Cell(index))
            # self.cell_index[self.cell_list[count].cell_id] = count
            # count += 1
        pass

    # add face data, those data may be used in the future
    def get_face_data(self, content):
        for index in content:
            self.face_list.append(Face(index))
        pass

    # add material data
    def get_material_data(self, content):
        material_index = []
        number = len(content)
        for index in range(number):
            if content[index][:1] == "m":
                material_index.append(index)
        material = []
        for count in range(len(material_index) - 1):
            material.append(" ".join(content[material_index[count]:material_index[count + 1]]))
        material.append(" ".join(content[material_index[len(material_index) - 1]:]))
        index = 0
        while index < len(material):
            if material[index][:2] == "mt":
                material[index - 1] += (" " + material[index])
                del material[index]
            else:
                index += 1
        material_count = 0
        for index in material:
            self.material_list.append(Material(index))
            isotope_number = len(self.material_list[material_count].isotopes)
            if isotope_number > self.isotopes_max_number:
                self.isotopes_max_number = isotope_number
            self.mat_index[self.material_list[material_count].mat_id] = material_count
            material_count += 1
        pass

    # add ksrc data, this data is just a copy of the input file part, we can get kcode and ksrc data from this
    def get_ksrc_data(self, content):
        self.ksrc = content
        pass

    # add tally lists, that is, the tally part of the MCNP input file
    def get_tally_data(self, content):
        tally_index = []
        number = len(content)
        for index in range(number):
            if content[index][:1] == "f":
                tally_index.append(index)
        tally = []
        for count in range(len(tally_index) - 1):
            tally.append(" ".join(content[tally_index[count]:tally_index[count + 1]]))
        tally.append(" ".join(content[tally_index[len(tally_index) - 1]:]))
        index = 0
        while index < len(tally):
            if tally[index][:2] == "fm":
                tally[index - 1] += (" " + tally[index])
                del tally[index]
            else:
                index += 1
        for index in tally:
            self.tally.append(Tally(index))
        pass

    # add energy grid data
    def get_energy_grid_data(self, content):
        content[0] = content[0][content[0].find("e0") + 2:].strip(" ")
        data_sheet = " ".join(content)
        while data_sheet:
            position = data_sheet.find(' ')
            if position >= 0:
                para = data_sheet[:position]
                data_sheet = data_sheet[position + 1:]
            else:
                para = data_sheet
                data_sheet = []
            if para is '':
                continue
            self.energy_grid.append(float(para))
        pass

    # organize the tally data and save the data in the corresponding cell of the cell list
    def organize_tally(self):
        # total flux, tally[0]
        for index in range(len(self.tally[0].cells)):
            inside_index = self.cell_index[self.tally[0].cells[index]]
            self.cell_list[inside_index].total_flux = self.tally[0].database[index]
        # tally data,tally[1],[2],[3]...
        for tally_index in range(1, len(self.tally)):
            group_count = 0
            for index in range(len(self.tally[tally_index].cells)):
                multiplier_num = len(self.tally[tally_index].multiplier)
                inside_index = self.cell_index[self.tally[tally_index].cells[index]]
                self.cell_list[inside_index].tally_data = self.tally[tally_index].database[
                                                          group_count:group_count + multiplier_num]
                group_count += multiplier_num
        # for index in range(len(self.tally)):
        #    self.tally[index].database = []
        pass

    def modification(self):
        # number of the isotopes added
        for index in self.cell_list:
            if index.material_id == 0:
                index.isotope_num = 0
            else:
                mat_id = self.mat_index[index.material_id]
                index.isotope_num = len(self.material_list[mat_id].isotopes)
        # the link between cells index in the cell list and the number of the cells
        index = 0
        while index < len(self.cell_list):
            if self.cell_list[index].importance == 0:
                self.void_cells_list.append(self.cell_list[index])
                del self.cell_list[index]
            else:
                index += 1
        self.cell_list.reverse()
        # cell index
        for count in range(len(self.cell_list)):
            self.cell_index[self.cell_list[count].cell_id] = count


# main function
def get_data_from_all_output_files():
    if len(sys.argv) == 1:
        file_list = glob.glob("./*.out")
    else:
        file_list = sys.argv[1:]
    total_number = len(file_list)
    # hdf5 database
    data_base = [[]] * total_number
    # loop begin
    for index in range(total_number):
        data_base[index] = [None, None]
        get_data_from_a_file(file_list[index], data_base[index])
        create_an_hdf5_file(file_list[index], data_base[index])
    # database written into the hdf5 file
    pass


# modify the input file part of the MCNP output file to be easy to get data from this part
# delete the warning and the content after "$" or "c" and the first line
# link the lines broken by "&"
def input_initial_process(content):
    flag = 0
    former = ""
    count = 0
    del content[0]
    for index in content:
        if index[:8] == "warning.":
            content.remove(index)
    while count < len(content):
        content[count] = content[count][content[count].find("-") + 1:].strip(" ")
        if content[count][:1] == "c":
            content.remove(content[count])
            continue
        if "$" in content[count]:
            content[count] = content[count][:content[count].find("$")].strip(" ")
        if flag == 1:
            content[count] = former + content[count]
            flag = 0
        if content[count][len(content[count]) - 1:] == "&":
            former = content[count].strip(" &") + " "
            content.remove(content[count])
            count -= 1
            flag = 1
        count += 1
    pass


# split the file into three parts (in fact 5 parts)
# all are transformed to lower case
def split_file(content):
    number = len(content)
    for index in range(number):
        content[index] = content[index].strip(" \n").lower()
    board1 = content.index("")
    head_content = content[:board1]
    board2 = content.index("", board1 + 1)
    name_content = content[board1 + 1:board2]
    board3 = content.index("", board2 + 1)
    input_content = content[board2 + 1:board3]
    input_initial_process(input_content)
    for index in range(board3 + 1, number):
        if "1tally" in content[index]:
            break
    keff_content = content[board3 + 1:index - 1]
    tally_content = content[index:]
    result = [input_content, keff_content, tally_content]
    return result


# get keff value and also the error of keff
def get_keff_from_output(content):
    keff = []
    for index in content:
        if index[:11] == "| the final":
            break
    while index:
        position = index.find(' ')
        if position >= 0:
            para = index[:position]
            index = index[position + 1:]
        else:
            para = index
            index = []
        if para == '':
            continue
        try:
            value = float(para)
        except ValueError:
            continue
        else:
            keff.append(value)
    return keff


# get dat from one MCNP output file
def get_data_from_a_file(file_name, data_base):
    input_file = open(file_name)
    input_content = input_file.readlines()
    input_file.close()
    # slice_input = [[], [], []]
    slice_input = split_file(input_content)
    input_data = Inputdata(slice_input[0])
    keff = get_keff_from_output(slice_input[1])
    total_length = len(slice_input[2])
    # cross section part, volumes should be added in this part
    count = 0
    cell_volumes = 0
    tally_count = 0
    end_flag = 0
    while tally_count < len(input_data.tally):
        # find "1tally"
        while True:
            lines = slice_input[2][count]
            count += 1
            if "1tally" in lines:
                break
            elif count >= total_length:
                end_flag = 1
                break
        if end_flag == 1:
            break
        # fill the volumes of cells
        while True:
            lines = slice_input[2][count]
            count += 1
            if "volumes" in lines:
                break
        if cell_volumes < len(input_data.cell_list):
            items = ''
            count -= 1
            while True:
                count += 2
                lines = slice_input[2][count]
                if "cell" in lines:
                    break
                items += (lines + " ")
            volumes = []
            while items:
                position = items.find(' ')
                if position >= 0:
                    para = items[:position]
                    items = items[position + 1:]
                else:
                    para = items
                    items = []
                if para is '':
                    continue
                volumes.append(float(para))
            for index in range(len(volumes)):
                # better refer to the id of the cells
                cell_index = input_data.cell_index[input_data.tally[tally_count].cells[index]]
                input_data.cell_list[cell_index].volume = volumes[index]
                cell_volumes += 1

        # extract the data of reaction rates, the error values are also extracted to the object of Class Inputdata
        while True:
            while True:
                count += 1
                lines = slice_input[2][count]
                if "energy" in lines:
                    break
            count += 1
            lines = slice_input[2][count]
            count += 1
            cross_sections_data = ""
            while lines:
                cross_sections_data += (" " + lines)
                lines = slice_input[2][count]
                count += 1
            cell_data_groups = []
            cell_data_group = []
            member_count = 0
            cross_sections_data.strip(" ")
            while cross_sections_data:
                position = cross_sections_data.find(' ')
                if position >= 0:
                    para = cross_sections_data[:position]
                    cross_sections_data = cross_sections_data[position + 1:]
                else:
                    para = cross_sections_data
                    cross_sections_data = []
                if para is '':
                    continue
                try:
                    data = float(para)
                except ValueError:
                    if para == "total":
                        data = 0.0
                    else:
                        print "Error"
                        exit()
                cell_data_group.append(data)
                member_count += 1
                if member_count == 3:
                    cell_data_groups.append(cell_data_group)
                    member_count = 0
                    cell_data_group = []
            input_data.tally[tally_count].database.append(cell_data_groups)
            if not slice_input[2][count]:
                break
        tally_count += 1
    input_data.organize_tally()
    data_base[0] = input_data
    data_base[1] = keff  # data_base: [0]input_data [1]keff


# create one hdf5 file from the data of the object of Class Inputdata
def create_an_hdf5_file(file_name, data_base):
    # head editing
    # driver="core" will fill the data when the file is closed, which may be more effective
    hdf5_file = h5py.File(file_name + ".h5", "w", driver="core")
    state1 = hdf5_file.create_group("STATE_0001")
    reaction_rate = state1.create_group("reaction_rate_details")
    energy_groups_number = len(data_base[0].energy_grid)
    energy_grid = np.array(list(reversed(data_base[0].energy_grid))) * 1.0E+06
    cells_number = len(data_base[0].cell_list)
    isotopes_number = data_base[0].isotopes_max_number
    # containers
    rr_isoID = np.zeros((1, 1, 1, 1, cells_number, isotopes_number), dtype="i")
    rr_mID = np.zeros((1, 1, 1, 1, cells_number), dtype="i")
    rr_numDen = np.zeros((1, 1, 1, 1, cells_number, isotopes_number))
    tally_shape = (1, 1, 1, 1, cells_number, energy_groups_number, isotopes_number)
    rr_abs = np.zeros(tally_shape)
    rr_scatter = np.zeros(tally_shape)
    rr_nufiss = np.zeros(tally_shape)
    rr_flux = np.zeros((1, 1, 1, 1, cells_number, energy_groups_number))
    rr_volume = np.zeros((1, 1, 1, 1, cells_number))
    rr_n2n = np.zeros(tally_shape)
    rr_n3n = np.zeros(tally_shape)
    keff = np.array(data_base[1])
    # fill the containers part1
    for cel in range(cells_number):
        rr_volume[0][0][0][0][cel] = data_base[0].cell_list[cel].volume
        rr_mID[0][0][0][0][cel] = data_base[0].cell_list[cel].material_id
    for cel in range(cells_number):
        isotope_num = data_base[0].cell_list[cel].isotope_num
        for iso in range(isotope_num):
            inside_mat_id = data_base[0].mat_index[rr_mID[0][0][0][0][cel]]
            temp_mat = data_base[0].material_list[inside_mat_id].isotopes
            rr_isoID[0][0][0][0][cel][iso] = temp_mat[iso][0] * 10
            rr_numDen[0][0][0][0][cel][iso] = temp_mat[iso][1]
    # fill the containers part2
    for cel in range(cells_number):
        for e_range in range(energy_groups_number):
            reverse_e_range = energy_groups_number - e_range - 1
            rr_flux[0][0][0][0][cel][e_range] = data_base[0].cell_list[cel].total_flux[reverse_e_range][1]
            isotope_num = data_base[0].cell_list[cel].isotope_num
            isotope_num = min(isotope_num, len(data_base[0].cell_list[cel].tally_data) / 5)
            for iso in range(isotope_num):
                rr_abs[0][0][0][0][cel][e_range][iso] = data_base[0].cell_list[cel].tally_data[iso * 5 + 0][
                                                            reverse_e_range][1] * rr_numDen[0, 0, 0, 0, cel, iso]
                rr_nufiss[0][0][0][0][cel][e_range][iso] = data_base[0].cell_list[cel].tally_data[iso * 5 + 1][
                                                               reverse_e_range][1] * rr_numDen[0, 0, 0, 0, cel, iso]
                rr_scatter[0][0][0][0][cel][e_range][iso] = data_base[0].cell_list[cel].tally_data[iso * 5 + 2][
                                                                reverse_e_range][1] * rr_numDen[0, 0, 0, 0, cel, iso]
                rr_n2n[0][0][0][0][cel][e_range][iso] = data_base[0].cell_list[cel].tally_data[iso * 5 + 3][
                                                            reverse_e_range][1] * rr_numDen[0, 0, 0, 0, cel, iso]
                rr_n3n[0][0][0][0][cel][e_range][iso] = data_base[0].cell_list[cel].tally_data[iso * 5 + 4][
                                                            reverse_e_range][1] * rr_numDen[0, 0, 0, 0, cel, iso]
    # modification because of (n, 2n) and (n, 3n)
    rr_abs = rr_abs - rr_n2n - 2 * rr_n3n
    rr_scatter = rr_scatter + 2 * rr_n2n + 3 * rr_n3n

    # create datasets in the hdf5 file
    reaction_rate.create_dataset("rr_vol", data=rr_volume)
    reaction_rate.create_dataset("rr_mID", data=rr_mID)
    reaction_rate.create_dataset("rr_isoID", data=rr_isoID)
    reaction_rate.create_dataset("rr_numDen", data=rr_numDen)
    reaction_rate.create_dataset("rr_flux", data=rr_flux)
    reaction_rate.create_dataset("rr_abs", data=rr_abs)
    reaction_rate.create_dataset("rr_scatter", data=rr_scatter)
    reaction_rate.create_dataset("rr_nufiss", data=rr_nufiss)
    reaction_rate.create_dataset("rr_eMesh", data=np.array(energy_grid))
    state1.create_dataset("keff", data=keff[0:1])
    hdf5_file.close()
    pass


# main function
# if there is no parameter other than the name of this script when calling this script, it will detect all *.out files
#  and extract their data into corresponding hdf5 file
# if there is any parameter after the name of the script, those parameters will be regarded as the input files
get_data_from_all_output_files()
