"""ortho_processing.py: takes data from inParanoid and processes it for use on the FE"""
__author__ = "Anna Sheaffer"
__email__ = "asheaffe@iwu.edu"
__date__ = "May 22, 2023"

from classes import Orthogroup
import python.utility as utility

def add_orthogroup(num:str, protein:utility, prev:Orthogroup):
    """Function for checking if the previous orthogroup has the same number as the current,
    and adds an orthogroup to the list of orthogroups if they do, otherwise return the new orthogroup object
    
    :param num: entry number from current line in inParanoid file to be checked
    :param current: current protein object to be checked
    :param prev: previous protein object to
    :return: Orthogroup object which will either be the updated version of the previous or a new Orthogroup object that is created"""
    
    # if the previous Orthogroup object has the same group number as the current
    if prev.is_same_group(num):
        # Add the current ortholog to the previous
        prev.add_ortholog(protein)
        val = prev
    else:
        # Otherwise, create a new Orthogroupcurrent
        val = Orthogroup.Orthogroup(num, protein)
    return val

def main():
    s1_input = 'yeast'
    s2_input = 'worm'

    yeast_ensembl = utility.get_ensembl_data(s1_input)
    worm_ensembl = utility.get_ensembl_data(s2_input)

    # #print(worm_ensembl)

    # used to find associated Protein object for each entry
    swissprot_lookup_s1 = {y: yeast_ensembl[x] for x in yeast_ensembl for y in yeast_ensembl[x].swissprots}
    trembl_lookup_s1 = {y: yeast_ensembl[x] for x in yeast_ensembl for y in yeast_ensembl[x].trembls}

    swissprot_lookup_s2 = {y: worm_ensembl[x] for x in worm_ensembl for y in worm_ensembl[x].swissprots}
    trembl_lookup_s2 = {y: worm_ensembl[x] for x in worm_ensembl for y in worm_ensembl[x].trembls}

    #print(swissprot_lookup_s1['Q3E7C1'])
    # print(trembl_lookup['Q7YWY8'])
    # print(trembl_lookup['Q7YWY8'].gene_id)

    # map each uniprot id to ensembl id using the above dictionaries
    # keep track of any unmappable ids and include in output

    filepath = "../../inParanoid-files/c_elegans-s_cerevisiae_inparanoidb-9.fa"

    orthogroups = {}
    species1 = ''
    species2 = ''
    unmappable = 0  # for counting number of ids that are unmappable to a Protein object
    unmappable_ids = []
    with open(filepath, 'r') as f1:
        content = f1.readlines()

        line1 = content[0].split('\t')

        # set species1 based on content
        species1 = line1[2].split('.')[0]
        orthogroups[species1] = ''

        prev_orthogroup = None

        with open('test_orthology_data.txt', 'w') as f2:

            f2.write("! Orthology data for " + s1_input + " and " + s2_input + "\n")

            for line in content:
                line = line.strip()
                line = line.split('\t')
                
                # defining the current species within the loop for readability
                current_species = line[2].split('.')[0]
                current_protein = line[4]
                entry_num = line[0]

                if species2 == '' and current_species != species1:
                    species2 = current_species
                    orthogroups[species1] = {species2: []}
                
                protein_obj = None
                current_orthogroup = None
                # check for species 1 or species 2
                if current_species == species1:
                    
                    # map uniprot to ensembl id
                    if current_protein in swissprot_lookup_s1:
                        protein_obj = swissprot_lookup_s1[current_protein]

                        if current_orthogroup is None:
                            current_orthogroup = Orthogroup.Orthogroup(entry_num, protein_obj)

                        if prev_orthogroup is not None:
                            current_orthogroup = add_orthogroup(entry_num, protein_obj, prev_orthogroup)

                    elif current_protein in trembl_lookup_s1:
                        protein_obj = trembl_lookup_s1[current_protein]

                        if current_orthogroup is None:
                            current_orthogroup = Orthogroup.Orthogroup(entry_num, protein_obj)

                        if prev_orthogroup is not None:
                            current_orthogroup = add_orthogroup(entry_num, protein_obj, prev_orthogroup)
                    else:
                        unmappable += 1
                        unmappable_ids.append(current_protein)
                        f2.write("! " + current_protein + "\tspecies:" + species1 + "\n")    # 559292 is species 1 (yeast)
                
                elif current_species == species2:
                    
                    # map uniprot to ensembl id
                    if current_protein in swissprot_lookup_s2:
                        protein_obj = swissprot_lookup_s2[current_protein]

                        if current_orthogroup is None:
                            current_orthogroup = Orthogroup.Orthogroup(entry_num, protein_obj)

                        if prev_orthogroup is not None:
                            current_orthogroup = add_orthogroup(entry_num, protein_obj, prev_orthogroup)
                    elif current_protein in trembl_lookup_s2:
                        protein_obj = trembl_lookup_s2[current_protein]

                        if current_orthogroup is None:
                            current_orthogroup = Orthogroup.Orthogroup(entry_num, protein_obj)

                        if prev_orthogroup is not None:
                            current_orthogroup = add_orthogroup(entry_num, protein_obj, prev_orthogroup)
                    else:
                        unmappable += 1
                        unmappable_ids.append(current_protein)
                        f2.write("! " + current_protein + "\tspecies:" + species2 + "\n")    # 6239 is species 2 (worm)

                #print(prev_orthogroup.is_same_group(current_orthogroup))
                # add orthogroup to dictionary if the prev has a different entry num compared to the current
                if prev_orthogroup is not None and prev_orthogroup.is_same_group(current_orthogroup) is False:
                    orthogroups[species1][species2].append(prev_orthogroup.group)

                prev_orthogroup = current_orthogroup

            f2.write("! " + str(unmappable) + " unmappable IDs, " + str(len(content)) + " total entries\n")
            f2.write("! " + str(unmappable/len(content)*100) + "% of all entries in file are unmappable\n")
            #print(orthogroups)

            for group in orthogroups[species1][species2]:
                for item in group:
                    if item.gene_id == 'WBGene00004187':
                        print(item)
                    f2.write(item.gene_id + "\t")

                    # if i < len(group) - 1:
                    #     f2.write("\t")

                f2.write("\n")

            print("completed")
main()