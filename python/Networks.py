
def main():
    # read_biogrid("C:\\Users\\annsb\\OneDrive\\Documents\\PPI-Network-Alignment\\ensembl\\worm_protein_ids106.txt")
    dict = mp_dict("C:\\Users\\annsb\\OneDrive\\Documents\\PPI-Network-Alignment\\test_obo.obo")
    print(dict)



def get_id(line):
    """
    takes in a line of text and reads the id

    :param line: line of text as str
    :return: id as str if contains, else return None
    """

    # checks if the line starts with 'id:'
    if line[:4] == "id: ":
        return line[7:-1]

    else:
        return '!'

def find_inter_id(line):
    """
    function reads a block of text and extracts an interaction id using 'is_a'

    :param line: line of text from file as str
    :return: MPid as str or none
    """
    # strip the line of /n char
    line.strip()

    # check for 'is_a'
    if line[:4] == 'is_a':
        # if true, append the id to the list
        return line[9:13]
    # else, return nonetype
    else:
        return None

def if_term(line):
    """
    function checks if the line contains '[Term]'

    :param line: line of text from file as str
    :return: boolean true if contains or false if not
    """

    # strip the line of any /n char
    line.strip()

    # check if the line is equal to '[Term]'
    if line[:6] == '[Term]':
        return True

    return False

def mp_dict(text_file):
    """
    function iterates through the obo file to make a dictionary of mp ids with interaction ids

    :param text_file: text file read as str
    :return: dictionary of ids with interactions
    """

    # read the text file
    with open(text_file, encoding='utf-8') as f:
        contents = f.readlines()

    # dictionary to be returned
    id_dict = {}

    # temporary list of multiple interaction ids
    temp_lst = []

    mp_id = ''

    # loop through the text file
    for line in contents:
        # assign the id
        if get_id(line) != '!':
            mp_id = get_id(line)

        # assign an interaction id
        temp = find_inter_id(line)

        # check if temp is none
        if temp is not None:
            # append temp to the temp list
            temp_lst.append(temp)

        # if mp_id is not None and temp_lst length > 0 and [Term]
        if len(temp_lst) > 0 and if_term(line):
            for parent in temp_lst:
                if parent not in id_dict:
                    id_dict[parent] = [mp_id]
                else:
                    id_dict[parent].append(mp_id)

            id_dict[mp_id] = temp_lst

            # reset temp_lst to be empty
            temp_lst = []

        elif if_term(line) and len(temp_lst) == 0:
            id_dict[mp_id] = []

            temp_lst = []

    # return the dictionary of ids
    return id_dict
def is_exp(dict, mp_id):
    """
    takes an MP id and checks if it has an experimental interaction

    :param dict: dictionary of mp ids with their connections as the definition
    :param mp_id: the mp id that will be checked
    :return: boolean True if original mp id is contained in the list, false if not
    """

    # base case 1 - if the mp id is '0407' return true
    if mp_id == '0045':
        return True

    # base case 2 - if the input is empty return false
    if mp_id == []:
        return False

    # otherwise, continue recursion
    for code in dict[mp_id]:
        if is_exp(dict, code):
            return True

def is_phys(dict, mp_id):
    """
    function takes in the dictionary and an MP id and recursively checks if it has a physical interaction

    :param dict: dictionary of mp ids with their connections as the definition
    :param mp_id: the mp id that will be checked
    :return: boolean True if original mp id is contained in the list, false if not
    """

    # base case 1 - if the mp id is '0407' return true
    if mp_id == '0407':
        return True

    # base case 2 - if the input is empty return false
    if mp_id == []:
        return False

    # otherwise, continue recursion
    for code in dict[mp_id]:
        if is_phys(dict, code):
            return True

    return False

# def read_biogrid(filename):
#     gene_dict = {}
#
#     with open(filename, 'r') as f1:
#         data = f1.readlines()
#
#         for line in data:
#             line = line.strip()
#             eles = line.split("\t")
#
#             print(eles)
#
#             # create a dictionary with entrezgene as the key
#             entrez = eles.pop(1)
#
#             gene_dict[entrez] = eles


main()