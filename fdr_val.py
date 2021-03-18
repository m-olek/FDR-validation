#!/usr/bin/env python
import decimal
import argparse
import sys
import os
import collections
import numpy as np
import gemmi
import mrcfile
#import pickle
#import pandas as pd
import math
from gemmi import cif

class getMap():
    def __init__(self, confmap):

        with mrcfile.open(confmap, 'r') as mrc:
            self.map = mrc.data
            self.apix = mrc.voxel_size.item()
            self.orix = mrc.header.origin.x
            self.oriy = mrc.header.origin.y
            self.oriz = mrc.header.origin.z
            self.map_dim_x = mrc.header.nx
            self.map_dim_y = mrc.header.ny
            self.map_dim_z = mrc.header.nz

def structure_in(path):
    if(
            path[-4:] == ".cif"
            or path[-6:] == ".mmcif"
    ):
        doc = gemmi.cif.read(path)[0]
        structure = gemmi.make_structure_from_block(doc)
    else:
        structure = gemmi.read_structure(path)

    structure.remove_empty_chains()
    structure.remove_hydrogens()
    return structure

def structure_out(path, structure):
    if(
            path[-4:] == ".cif"
    ):


        structure.make_mmcif_document().write_file(path[:-4]+'validated_NEW.cif')
    elif(
            path[-6:] == ".mmcif"
    ):
        structure.make_mmcif_document().write_file(path[:-6]+'validated_NEW.mmcif')
    else:
        structure.write_pdb(str(path[:-4]+'_validated_NEW.pdb'))
        #structure = gemmi.read_structure(path)

def validate_CAs(path_pdb, path_conf_map):
    print('Validating based on CA positions...')
    with mrcfile.open(path_conf_map, 'r') as mrc:
        map = mrc.data
        apix = mrc.voxel_size.item()
        orix = mrc.header.origin.x
        oriy = mrc.header.origin.y
        oriz = mrc.header.origin.z
        map_dim_x = mrc.header.nx
        map_dim_y = mrc.header.ny
        map_dim_z = mrc.header.nz
        temp_map = np.zeros(map_dim_x*map_dim_y*map_dim_z, dtype=np.float16).reshape(map_dim_x, map_dim_y, map_dim_z)

    chains2remove2 = []
    residues2remove2 = []
    resnames2rem = []
    lig_list = []
    water_list = []
    structure_in = gemmi.read_structure(path_pdb)
    temp_valCA=0
    temp_valxx = 0
    c = 0

    for model in structure_in:

        for chain in model:
            lignd = chain.get_ligands()
            #print(lignd)
            for lig in lignd:
                lig_list.append(lig.seqid.num)
            #print(lig_list)
            waters = chain.get_waters()
            for water in waters:
                water_list.append(water.seqid.num)

            for residue in chain:
                if residue.name in ['A','C','T','G','U']:
                    for atom in residue:
                        if (atom.name == 'C1\''):
# for nucleic acids: P, OP1, OP2, O5',C5', C4', C3',C2', C1', o4', o3',o2' -> average
                            x_pos = ((atom.pos.x-orix)/apix[0])
                            y_pos = ((atom.pos.y-oriy)/apix[1])
                            z_pos = ((atom.pos.z-oriz)/apix[2])
                            temp_valCA = map[int(round(z_pos)),int(round(y_pos)),int(round(x_pos))]
                            #print('ID: ' + str(residue.seqid)+' x: '+str(int(x_pos))+' y:'+str(y_pos)+' z:'+str(z_pos))

                elif residue.seqid.num in lig_list:
                    for atom in residue:
                        x_pos = ((atom.pos.x-orix)/apix[0])
                        y_pos = ((atom.pos.y-oriy)/apix[1])
                        z_pos = ((atom.pos.z-oriz)/apix[2])

                        c += 1
                        temp_valxx += map[int(round(z_pos)),int(round(y_pos)),int(round(x_pos))]
                        temp_valCA = temp_valxx/c

                    c= 0
                    temp_valxx = 0

                elif residue.seqid.num in water_list:
                    for atom in residue:
                        x_pos = ((atom.pos.x-orix)/apix[0])
                        y_pos = ((atom.pos.y-oriy)/apix[1])
                        z_pos = ((atom.pos.z-oriz)/apix[2])

                        c += 1
                        temp_valxx += map[int(round(z_pos)),int(round(y_pos)),int(round(x_pos))]
                        temp_valCA = temp_valxx/c

                    c= 0
                    temp_valxx = 0
                else:
                    for atom in residue:
                        if (atom.name == 'CA'):
                            #atom in C, CA, N, O -> average value for a backbone, add as an option
#PDB file CONFMAP -backbone_val -prune

                            x_pos = ((atom.pos.x-orix)/apix[0]) #atom.pos.x * scale
                            y_pos = ((atom.pos.y-oriy)/apix[1]) #atom.pos.y * scale
                            z_pos = ((atom.pos.z-oriz)/apix[2]) #atom.pos.z * scale

                            #xpos_rem = ((atom.pos.x-orix)%apix[0])
                            #if (xpos_rem < 0.2):
                            #    x_pos = x_pos + 1
                            #elif (xpos_rem > 0.8):

                            temp_valCA = map[int(round(z_pos)),int(round(y_pos)),int(round(x_pos))]

                for atom in residue:
                    atom.b_iso = temp_valCA

                temp_valCA = 0

    split_name = os.path.splitext(os.path.basename(path_pdb));

    #map_basename = os.path.splitext(os.path.basename(path_pdb))[0]
    #self.confidence_map = os.path.join(splitFilename[0] + '_confidenceMap.mrc')

    structure_in.write_pdb(split_name[0]+'_CA_validated.pdb')

#    structure_in.write_pdb(str(path_pdb[:-4]+'_rounded.pdb'))
#    with mrcfile.new(path_pdb[:-4]+'_map.mrc',overwrite=True) as new_map:
#        new_map.set_data(nowa2)
    print('Validation based on CA done.')

def pruneCAs(path_pdb,path_conf_map):
    print('Pruning...')
    x123 = getMap(path_conf_map)

    chains2remove2 = []
    residues2remove2 = []
    resnames2rem = []

    structure_in = gemmi.read_structure(path_pdb)

    for model in structure_in:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if (atom.name == "CA"): 
                        x_pos = ((atom.pos.x-x123.orix)/x123.apix[0]) #atom.pos.x * scale
                        y_pos = ((atom.pos.y-x123.oriy)/x123.apix[1]) #atom.pos.y * scale
                        z_pos = ((atom.pos.z-x123.oriz)/x123.apix[2]) #atom.pos.z * scale
                        if (x123.map[int(round(z_pos)),int(round(y_pos)),int(round(x_pos))] != 1.0):

                            chains2remove2.append(str(chain.name))
                            residues2remove2.append(str(int(str(residue.seqid.num))-1))
                            resnames2rem.append(str(residue.name))

                            chains2remove2.append(str(chain.name))
                            residues2remove2.append(str(residue.seqid.num))
                            resnames2rem.append(str(residue.name))

                            chains2remove2.append(str(chain.name))
                            residues2remove2.append(str(int(str(residue.seqid.num))+1))
                            resnames2rem.append(str(residue.name))

    unique_res = collections.OrderedDict.fromkeys(zip(chains2remove2,residues2remove2,resnames2rem))

    split_name = os.path.splitext(os.path.basename(path_pdb));
    path_out = split_name[0]+'_BACKBONE_average.pdb'
    structout = gemmi.read_structure(path_out)
    removed_residues = []
    res_iter =0
    for model in structout:
        for chain in model:
            for residue in chain:
                for k in unique_res:
                    if((str(chain.name) == str(k[0])) and (str(residue.seqid.num) == str(k[1]) and (str(residue.seqid.icode)==' '))):
                        removed_residues.append(str(k[0])+'\t\t'+str(k[1])+'\t\t'+str(k[2])+'\n')
                        res_iter = (len(residue))
                        for i in range(res_iter):
                            del residue[0]
                        res_iter = 0
                        break
    structout.write_pdb(split_name[0]+'_PRUNED.pdb')

    out_path = split_name[0] + '_FDR_pruned_removed_residues.txt'
    with open(out_path, 'w') as fp:
        fp.write('=============================== \
        \nResidues removed from the model \
        \n=============================== \
        \nChain ID - Residue ID - Residue name\n')
        fp.write(''.join('%s' % x for x in removed_residues))

    print('Removed ' + str(len(removed_residues))+ ' residues')
    print('Pruning done.')

def validate_backbone(path_pdb, path_conf_map):
    print('Backbone validation started...')
    x123 = getMap(path_conf_map)

    chains2remove2 = []
    residues2remove2 = []
    resnames2rem = []
    lig_list = []
    water_list = []
    structure_inn = structure_in(path_pdb)
    structure_inn2 = structure_in(path_pdb)
    temp_valCA=0
    temp_valxx = 0
    c = 0
    for_csv = []
    for model in structure_inn:
        for chain in model:
            lignd = chain.get_ligands()
            for lig in lignd:
                lig_list.append(lig.seqid.num)
            waters = chain.get_waters()
            for water in waters:
                water_list.append(water.seqid.num)
            for residue in chain:
                if residue.name in ['A','C','T','G','U']:
                    for atom in residue:
                        if (atom.name in ['C1\'','C2\'','C3\'','C4\'', 'C5\'', 'P', 'O3\'', 'O4\'', 'O5\'']):
                            x_pos = ((atom.pos.x-x123.orix)/x123.apix[0])
                            y_pos = ((atom.pos.y-x123.oriy)/x123.apix[1])
                            z_pos = ((atom.pos.z-x123.oriz)/x123.apix[2])
                            c += 1
                            temp_valxx += x123.map[int(round(z_pos)),int(round(y_pos)),int(round(x_pos))]
                            temp_valCA = temp_valxx/c
                    c=0
                    temp_valxx = 0
                elif residue.seqid.num in lig_list:
                    for atom in residue:
                        x_pos = ((atom.pos.x-x123.orix)/x123.apix[0])
                        y_pos = ((atom.pos.y-x123.oriy)/x123.apix[1])
                        z_pos = ((atom.pos.z-x123.oriz)/x123.apix[2])
                        c += 1
                        temp_valxx += x123.map[int(round(z_pos)),int(round(y_pos)),int(round(x_pos))]
                        temp_valCA = temp_valxx/c
                    c= 0
                    temp_valxx = 0
                elif residue.seqid.num in water_list:
                    for atom in residue:
                        x_pos = ((atom.pos.x-x123.orix)/x123.apix[0])
                        y_pos = ((atom.pos.y-x123.oriy)/x123.apix[1])
                        z_pos = ((atom.pos.z-x123.oriz)/x123.apix[2])
                        c += 1
                        temp_valxx += x123.map[int(round(z_pos)),int(round(y_pos)),int(round(x_pos))]
                        temp_valCA = temp_valxx/c
                    c= 0
                    temp_valxx = 0
                else:
                    for atom in residue:
                        if (atom.name in ['CA', 'N', 'C']):
                            x_pos = ((atom.pos.x-x123.orix)/x123.apix[0]) #atom.pos.x * scale
                            y_pos = ((atom.pos.y-x123.oriy)/x123.apix[1]) #atom.pos.y * scale
                            z_pos = ((atom.pos.z-x123.oriz)/x123.apix[2]) #atom.pos.z * scale
                            temp_valxx +=  x123.map[int(round(z_pos)),int(round(y_pos)),int(round(x_pos))]
                            c += 1
                            temp_valCA = temp_valxx/c
                    c=0
                    temp_valxx = 0
                for atom in residue:
                    atom.b_iso = temp_valCA
                for_csv.append(str(chain.name)+','+str(residue.seqid)+','+str(residue.name)+','+str(atom.b_iso)+'\n') 
                c=0
                temp_valCA = 0
    split_name = os.path.splitext(os.path.basename(path_pdb));
    out_path_csv = split_name[0] + '_FDR_ranked_residues.csv'
    with open(out_path_csv, 'w') as fc:
        fc.write('Chain_name'+','+'residue_id'+','+'residue_name'+','+'conf_score'+'\n')
        fc.write(''.join('%s' % x for x in for_csv))

    structure_inn.write_pdb(split_name[0]+'_BACKBONE_average.pdb')
    path_out = split_name[0]+'_BACKBONE_average.pdb'
    print('Backbone validation done.')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input_pdb",help="Input pdb file",type = str)
    parser.add_argument("-input_map",help="FDR map to validate the model",\
            type = str)
    parser.add_argument("-prune", "--prune", help="pruning mode", action = "store_true")
    parser.add_argument("-valCA", "--valCA", help="validate by CA position", action = "store_true")
    args = parser.parse_args()
    path_pdb = args.input_pdb
    path_conf_map = args.input_map
    if (args.prune and args.valCA):
        validate_backbone(path_pdb,path_conf_map)
        pruneCAs(path_pdb,path_conf_map)
        validate_CAs(path_pdb,path_conf_map)
    elif args.prune:
        validate_backbone(path_pdb,path_conf_map)
        pruneCAs(path_pdb,path_conf_map)
    elif args.valCA:
        validate_backbone(path_pdb,path_conf_map)
        validate_CAs(path_pdb,path_conf_map)
    else:
        validate_backbone(path_pdb,path_conf_map) 

if __name__ == '__main__':
    sys.exit(main())
