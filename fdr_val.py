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
    '''
    Read in input model with as a gemmi structure
    '''

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
    '''
    Defines output file format to match the input(pdb/cif/mmcif).
    '''
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

def remap_FDR(x,y,z,map_ref,element,mode,resolution):
    '''
        Calculate average confidence score at the atom position

        Args:
            x,y,z (float): atom coordinates
            map_ref (): FDR map used as a reference to get the confidence values at the atomic position
            element (float): Van der Walls radius for atom (from gemmi)
            mode: specifies the radius of the region of interest for FDR score calculation
            resolution (float): map resolution used for confidence score calculation in Sigma Threshold(ST) mode
        Returns:
            temp_val: Confidence score for an atom
    '''
    if mode.lower() == '1a':
        val_radius = 1.0
    elif mode.lower() == 'vdw':
        val_radius = element
    elif mode.lower() == 'st':
        val_radius = 2.5*0.225*resolution
    x_pos = ((x-map_ref.orix)/(map_ref.apix[0]))
    y_pos = ((y-map_ref.oriy)/(map_ref.apix[1]))
    z_pos = ((z-map_ref.oriz)/(map_ref.apix[2]))

    z_pos_min = int(round(z_pos-(val_radius/(map_ref.apix[2]))))
    z_pos_min = min(z_pos_min,int(z_pos))
    z_pos_max = int(round(z_pos+(val_radius/(map_ref.apix[2]))))
    z_pos_max = max(z_pos_max,int(z_pos))
    list_zpos = list(range(int(z_pos_min),int(round(z_pos_max))+1))

    y_pos_min = int(round(y_pos-(val_radius/(map_ref.apix[2]))))
    y_pos_min = min(y_pos_min,int(y_pos))
    y_pos_max = int(round(y_pos+(val_radius/(map_ref.apix[2]))))
    y_pos_max = max(y_pos_max,int(y_pos))
    list_ypos = list(range(int(y_pos_min),int(round(y_pos_max))+1))

    x_pos_min = int(round(x_pos-(val_radius/(map_ref.apix[2]))))
    x_pos_min = min(x_pos_min,int(x_pos))
    x_pos_max = int(round(x_pos+(val_radius/(map_ref.apix[2]))))
    x_pos_max = max(x_pos_max,int(x_pos))
    list_xpos = list(range(int(x_pos_min),int(round(x_pos_max))+1))

    dst = int(round(val_radius/round(map_ref.apix[2],3)))
    list_scores = []
    for z in list_zpos:
        for y in list_ypos:
            for x in list_xpos:
                if np.sqrt((x-x_pos)**2+(y-y_pos)**2+(z-z_pos)**2) <= dst:
                    list_scores.append(map_ref.map[z,y,x])
    temp_val =  sum(list_scores)/len(list_scores)
    return temp_val

def prune_res(path_pdb,path_conf_map,val_coord,map_res):
    '''
    Remove residues if any of the CA, C and N atoms score lower than 0.9
    Write out TXT file with IDs and names of removed residues
    and PDB file with pruned atomic model

    Args:
         path_pdb (str): path to the atomic model to validate
         path_conf_map (str): path to the confidence map used as a reference for validation
         val_coord (str): specify mode of FDR score calculation
         map_res (float): map resolution for Sigma Threshold(ST) mode

    '''
    temp_map = getMap(path_conf_map)

    chains2remove2 = []
    residues2remove2 = []
    resnames2rem = []

    structure_in = gemmi.read_structure(path_pdb)

    for model in structure_in:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if (atom.name in ['CA', 'N', 'C']):

                        atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map,atom.element.vdw_r,val_coord, map_res)

                        if (atom_score <0.9):

                            if residue.seqid.num > 1:
                                chains2remove2.append(str(chain.name))
                                residues2remove2.append(str(int(str(residue.seqid.num))-1))
                                resnames2rem.append(str(chain[residue.seqid.num-2].name))

                            chains2remove2.append(str(chain.name))
                            residues2remove2.append(str(residue.seqid.num))
                            resnames2rem.append(str(chain[residue.seqid.num-1].name))

                            if residue.seqid.num < len(chain):
                                chains2remove2.append(str(chain.name))
                                residues2remove2.append(str(int(str(residue.seqid.num))+1))
                                resnames2rem.append(str(chain[residue.seqid.num].name))

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

def validate_backbone(path_pdb, path_conf_map, val_coord,map_res,val_CA,val_AA):
    '''
    Calculate the FDR scores for residues in the model
    Write out the CSV file with scores and TXT attribute file.

    Args:
        path_pdb (str): path to the atomic model to validate
        path_conf_map (str): path to the confidence map used as a reference for validation
        val_coord (str): specify mode of FDR score calculation
        map_res (float): map resolution for Sigma Threshold(ST) mode
        val_CA (bool): True for additional C-alpha validation, False to skip(default)
        val_AA (bool): True for additional all atoms validation, False to skip(default)

    '''

    temp_map = getMap(path_conf_map)
    chains2remove2 = []
    residues2remove2 = []
    resnames2rem = []
    lig_list = []
    water_list = []
    structure_inn = structure_in(path_pdb)
    temp_valCA = 0
    residue_score_bb = 0
    c = 0
    for_csv = []
    for_attribute=[]
    for_csvCA=[]
    for_attributeCA=[]
    for_csvAA=[]
    for_attributeAA=[]
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
                            atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map,atom.element.vdw_r,val_coord,map_res)
                            residue_score += atom_score
                            c += 1
                            temp_valCA = residue_score/c
                    c=0
                    residue_score = 0

                elif residue.seqid.num in lig_list:
                    for atom in residue:
                        atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map,atom.element.vdw_r,val_coord,map_res)
                        residue_score += atom_score
                        c += 1
                        temp_valCA = residue_score/c
                    c= 0
                    residue_score = 0

                elif residue.seqid.num in water_list:
                    for atom in residue:
                        atom_score = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map,atom.element.vdw_r,val_coord, map_res)
                        residue_score += atom_score
                        c += 1
                        temp_valCA = residue_score/c
                    c= 0
                    residue_score = 0

                else:
                    for atom in residue:
                        if (atom.name in ['CA', 'N', 'C']):
                            atom_score_bb = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map,atom.element.vdw_r,val_coord,map_res)
                            residue_score_bb += atom_score_bb
                            c += 1
                            temp_valCA = residue_score_bb/c
                        if (atom.name in ['CA']) and val_CA == True:
                            atom_score_CA = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map,atom.    element.vdw_r,val_coord,map_res)

                            for_csvCA.append(str(chain.name)+','+str(residue.seqid)+','+str(residue.name)+','+str(atom_score_CA)+'\n')
                            for_attributeCA.append('\t:'+str(residue.seqid)+'.'+str(chain.name)+'\t'+str(atom_score_CA)+'\n')

                        if val_AA == True:
                            atom_score_AA = remap_FDR(atom.pos.x,atom.pos.y,atom.pos.z,temp_map,atom.element.vdw_r,val_coord,map_res)
                            for_csvAA.append(str(chain.name)+','+str(residue.seqid)+','+           str(residue.name)+','+str(temp_valCA)+','+str(atom.name)+','+str(atom_score_AA)+'\n')
                            for_attributeAA.append('\t:'+str(residue.seqid)+'.'+str(chain.name)+'\t'+str(temp_valCA)+'\n')

                    c=0
                    residue_score_bb = 0
                for atom in residue:
                    atom.b_iso=temp_valCA

                for_csv.append(str(chain.name)+','+str(residue.seqid)+','+str(residue.name)+','+str(temp_valCA)+'\n')
                for_attribute.append('\t:'+str(residue.seqid)+'.'+str(chain.name)+'\t'+str(temp_valCA)+'\n')
                c=0
                temp_valCA = 0


    if val_CA == True:
        out_pathCA = path_pdb[:-4] + '_FDR_ranked_CA_'+val_coord+'.csv'
        with open(out_pathCA, 'w') as fp:
            fp.write('Chain_name'+','+'residue_id'+','+'residue_name'+','+'conf_score'+'\n')
            fp.write(''.join('%s' % x for x in for_csvCA))

    if val_AA == True:
        out_pathCA = path_pdb[:-4] + '_FDR_ranked_AA_'+val_coord+'.csv'
        with open(out_pathCA, 'w') as fp:
            fp.write('Chain_name'+','+'residue_id'+','+'residue_name'+','+'conf_score'+'\n')
            fp.write(''.join('%s' % x for x in for_csvAA))


    out_pathBB = path_pdb[:-4] + '_FDR_ranked_residues_'+val_coord+'.csv'
    with open(out_pathBB, 'w') as fp:
        fp.write('Chain_name'+','+'residue_id'+','+'residue_name'+','+'conf_score'+'\n')
        fp.write(''.join('%s' % x for x in for_csv))

    out_pathBB_attr = path_pdb[:-4] + '_attribute_'+val_coord+'.txt'
    with open(out_pathBB_attr, 'w') as fp:
        fp.write('attribute: FDRscore\n')
        fp.write('match mode: 1-to-1\n')
        fp.write('recipient: residues\n')
        fp.write(''.join('%s' % x for x in for_attribute))


def check_res(map_res, val_coord):
    '''
    Terminate with error message if the resolution is not provided
    for Sigma Threshold(ST) mode.
    '''
    if val_coord.lower() == 'st' and map_res is None:
        print('ERROR: Please provide the map resolution for ST mode!')
        sys.exit()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_pdb",help="Input pdb file",type = str)
    parser.add_argument("input_map",help="FDR map to validate the model",\
            type = str)
    parser.add_argument("-coord_mode", "--coord_mode", help="Validation mode: 1A radius(default), VdW: Van der Waals radius, ST: sigma threshold", type = str,default='1A')
    parser.add_argument("-prune", "--prune", help="pruning mode", action = "store_true")
    parser.add_argument("-valCA", "--valCA", help="validate CA positions",default=False ,action = "store_true")
    parser.add_argument("-resolution","--resolution",help = 'Map resolution for ST mode',type=float)
    parser.add_argument("-valAA","--valAA", help = 'validate all atoms',default=False, action="store_true")
    args = parser.parse_args()
    path_pdb = args.input_pdb
    path_conf_map = args.input_map
    val_coord = args.coord_mode
    map_res = args.resolution
    val_CA = args.valCA
    val_AA = args.valAA
    check_res(map_res,val_coord)

    if args.prune:
        validate_backbone(path_pdb,path_conf_map,val_coord,map_res,val_CA,val_AA)
        prune_res(path_pdb,path_conf_map,val_coord,map_res)

    else:
        validate_backbone(path_pdb,path_conf_map,val_coord,map_res,val_CA,val_AA)

if __name__ == '__main__':
    sys.exit(main())

