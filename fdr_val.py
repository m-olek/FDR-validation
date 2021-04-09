#!/usr/bin/env python
import decimal
import argparse
import sys
import os
import collections
import numpy as np
import gemmi
import mrcfile
import pickle
import math
from gemmi import cif
import time

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

def structure_in(path: str) -> gemmi.Structure:
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
    structure.remove_alternative_conformations()

    return structure

def structure_out(path: str, structure: gemmi.Structure) -> None:
    if(
            path[-4:] == ".cif"
    ):
        structure.make_mmcif_document().write_file(path[:-4]+'validated.cif')
    
    elif(
            path[-6:] == ".mmcif"
    ):
        structure.make_mmcif_document().write_file(path[:-6]+'validated.mmcif')
    
    else:

        structure.write_pdb(str(path[:-4]+'_validated.pdb'))

def remap_FDR(x,y,z,x123):
    x_pos = ((x-x123.orix)/round(x123.apix[0],3)) 
    y_pos = ((y-x123.oriy)/round(x123.apix[1],3))
    z_pos = ((z-x123.oriz)/round(x123.apix[2],3))
    z_pos_min = int(round(z_pos-(1.0/round(x123.apix[2],3))))
    z_pos_min = min(z_pos_min,int(z_pos))
    z_pos_max = int(round(z_pos+(1.0/round(x123.apix[2],3))))
    z_pos_max = max(z_pos_max,int(z_pos))
    list_zpos = list(range(int(z_pos),int(round(z_pos))+1))
    y_pos_min = int(round(y_pos-(1.0/round(x123.apix[2],3))))
    y_pos_min = min(y_pos_min,int(y_pos))
    y_pos_max = int(round(y_pos+(1.0/round(x123.apix[2],3))))
    y_pos_max = max(y_pos_max,int(y_pos))
    list_ypos = list(range(int(y_pos),int(round(y_pos))+1))
    x_pos_min = int(round(x_pos-(1.0/round(x123.apix[2],3))))
    x_pos_min = min(x_pos_min,int(x_pos))
    x_pos_max = int(round(x_pos+(1.0/round(x123.apix[2],3))))
    x_pos_max = max(x_pos_max,int(x_pos))
    list_xpos = list(range(int(x_pos),int(round(x_pos))+1))

    list_scores = []

    for z in list_zpos:
        for y in list_ypos:
            for x in list_xpos:
                list_scores.append(x123.map[z,y,x])

    temp_valxx =  sum(list_scores)/len(list_scores)
    return temp_valxx

def validate_CAs(path_pdb, path_conf_map):
    temp_map = getMap(path_conf_map)   
    lig_list = []
    water_list = []
    structure_in = gemmi.read_structure(path_pdb)
    temp_valCA=0 
    residue_score = 0
    c = 0
    for_csv = []
    for_attribute=[]
    
    for model in structure_in:
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
                        if (atom.name == 'C1\''):
                            atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map)
                            residue_score += atom_score
                            c += 1
                            temp_valCA = residue_score/c
                    c=0
                    residue_score = 0

                elif residue.seqid.num in water_list:
                    for atom in residue:
                        atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map)
                        residue_score += atom_score
                        c += 1
                        temp_valCA = residue_score/c
                    c=0
                    residue_score = 0
                
                elif residue.seqid.num in lig_list:
                    for atom in residue:
                        atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map)
                        residue_score += atom_score
                        c += 1
                        temp_valCA = residue_score/c
                    c= 0
                    residue_score = 0

                else:
                    for atom in residue:
                        if (atom.name == 'CA'):
                            atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map)
                            residue_score += atom_score
                            c += 1
                            temp_valCA = residue_score/c
                    c=0
                    residue_score = 0

                for_csv.append(str(chain.name)+','+str(residue.seqid)+','+str(residue.name)+','+str(temp_valCA)+'\n')
                for_attribute.append('\t:'+str(residue.seqid)+'.'+str(chain.name)+'\t'+str(temp_valCA)+'\n')
                c=0
                temp_valCA = 0

    out_path2 = path_pdb[:-4] + '_FDR_ranked_residues_CA_only.csv'
    with open(out_path2, 'w') as fp:
        fp.write('Chain_name'+','+'residue_id'+','+'residue_name'+','+'conf_score'+'\n')
        fp.write(''.join('%s' % x for x in for_csv))

    out_path3 = path_pdb[:-4] + 'CA_only_attribute.txt'
    with open(out_path3, 'w') as fp:
        fp.write('attribute: FDRscore\n')
        fp.write('match mode: 1-to-1\n')
        fp.write('recipient: residues\n')
        fp.write(''.join('%s' % x for x in for_attribute))


    #structure_in.write_pdb(str(path_pdb[:-4]+'_CA_validated.pdb'))

def prune_res(path_pdb,path_conf_map):
    temp_map = getMap(path_conf_map)

    chains2remove2 = []
    residues2remove2 = []
    resnames2rem = []

    structure_in = gemmi.read_structure(path_pdb)

    for model in structure_in:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if (atom.name in ['CA', 'N', 'C']):#atom.name == "CA"):
                        
                        atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map)

                        if (atom_score <0.9):

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

    structout = gemmi.read_structure(path_pdb)
    removed_residues = []
    for model in structout:
        for chain in model:
            for residue in chain:
                for k in unique_res:
                    if((str(chain.name) == str(k[0])) and (str(residue.seqid.num) == str(k[1]) and (str(residue.seqid.icode)==' '))):
                        removed_residues.append(str(k[0])+' '+str(k[1])+' '+str(k[2])+'\n')
                        del residue[:]
                        break
    
    print('Removed ' + str(len(removed_residues))+ ' residues')
    out_path = path_pdb[:-4] + '_FDR_pruned_removed_residues.txt'
    with open(out_path, 'w') as fp:
        fp.write(''.join('%s' % x for x in removed_residues))
    
    structout.write_pdb(str(path_pdb[:-4]+'_FDR_pruned.pdb'))

def validate_backbone(path_pdb, path_conf_map):

    temp_map = getMap(path_conf_map)
    chains2remove2 = []
    residues2remove2 = []
    resnames2rem = []
    lig_list = []
    water_list = []
    structure_inn = structure_in(path_pdb)
    temp_valCA = 0 
    residue_score = 0
    c = 0
    for_csv = []
    for_attribute=[]
    list_ypos = []
    list_zpos = []
    list_xpos = []

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
                            atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map)
                            residue_score += atom_score
                            c += 1
                            temp_valCA = residue_score/c
                    c=0
                    residue_score = 0
                
                elif residue.seqid.num in lig_list:
                    for atom in residue:
                        atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map)
                        residue_score += atom_score
                        c += 1
                        temp_valCA = residue_score/c
                    c= 0
                    residue_score = 0
                
                elif residue.seqid.num in water_list:
                    for atom in residue:
                        atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map)
                        residue_score += atom_score
                        c += 1
                        temp_valCA = residue_score/c
                    c= 0
                    residue_score = 0
                
                else:
                    for atom in residue:
                        if (atom.name in ['CA', 'N', 'C']):
                            atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map)
                            residue_score += atom_score
                            c += 1
                            temp_valCA = residue_score/c
                    c=0
                    residue_score = 0
                for atom in residue:
                    atom.b_iso=temp_valCA

                for_csv.append(str(chain.name)+','+str(residue.seqid)+','+str(residue.name)+','+str(temp_valCA)+'\n')
                for_attribute.append('\t:'+str(residue.seqid)+'.'+str(chain.name)+'\t'+str(temp_valCA)+'\n')
                c=0
                temp_valCA = 0

    #structure_out(path_pdb, structure_inn)
            
            
    #print('TO CSV')
    out_path2 = path_pdb[:-4] + '_FDR_ranked_residues.csv'
    with open(out_path2, 'w') as fp:
        fp.write('Chain_name'+','+'residue_id'+','+'residue_name'+','+'conf_score'+'\n')
        fp.write(''.join('%s' % x for x in for_csv))

    out_path3 = path_pdb[:-4] + '_attribute.txt'
    with open(out_path3, 'w') as fp:
        fp.write('attribute: FDRscore\n')
        fp.write('match mode: 1-to-1\n')
        fp.write('recipient: residues\n')
        fp.write(''.join('%s' % x for x in for_attribute))

#    structure_inn.write_pdb(str(path_pdb[:-4]+'_BACKBONE_average.pdb'))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_pdb",help="Input pdb file",type = str)
    parser.add_argument("input_map",help="FDR map to validate the model",\
            type = str)
    parser.add_argument("-prune", "--prune", help="pruning mode", action = "store_true")
    parser.add_argument("-valCA", "--valCA", help="validate by CA position", action = "store_true")
    args = parser.parse_args()
    path_pdb = args.input_pdb
    path_conf_map = args.input_map

    if args.prune:
        prune_res(path_pdb,path_conf_map)
    elif args.valCA:
        validate_CAs(path_pdb,path_conf_map)
    else:
        validate_backbone(path_pdb,path_conf_map)

if __name__ == '__main__':
    sys.exit(main())

