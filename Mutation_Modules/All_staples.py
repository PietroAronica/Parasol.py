# All Stapler

import Frcmod_creator
import PDBHandler
import Leapy
from parmed.tools.actions import *
from parmed.amber.readparm import *

def parmed_command_AMX(vx1):
	bc = {}
        with open('Param_files/AminoAcid/ASN.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/AMX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vx1), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vx1), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vx1), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vx1), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vx1), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vx1), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vx1), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vx1), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@OD1'.format(vx1), bc['OD1']+((fc['OD1']-bc['OD1'])/10)*i).execute()
                change(parm, 'charge', ':{}@ND2'.format(vx1), bc['ND2']+((fc['ND2']-bc['ND2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vx1), bc['HD21']+((fc['HD2']-bc['HD21'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD22'.format(vx1), bc['HD22']-(bc['HD22']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vx1), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vx1), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def parmed_command_ABX(vx1):
	bc = {}
        with open('Param_files/AminoAcid/ABU.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/ABX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vx1), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vx1), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vx1), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vx1), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vx1), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vx1), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vx1), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vx1), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG1'.format(vx1), bc['HG1']-(bc['HG1']/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vx1), bc['HG2']+((fc['HG2']-bc['HG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vx1), bc['HG3']+((fc['HG3']-bc['HG3'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vx1), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vx1), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def parmed_command_CYX(vx1):
	bc = {}
        with open('Param_files/AminoAcid/CYS.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/CYX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vx1), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vx1), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vx1), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vx1), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vx1), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vx1), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vx1), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@SG'.format(vx1), bc['SG']+((fc['SG']-bc['SG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG'.format(vx1), bc['HG']-(bc['HG']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vx1), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vx1), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def parmed_command_MCX(vx1):
	bc = {}
        with open('Param_files/AminoAcid/MCY.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/MCX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vx1), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vx1), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vx1), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vx1), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vx1), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vx1), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vx1), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@SG'.format(vx1), bc['SG']+((fc['SG']-bc['SG'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vx1), bc['CD']+((fc['CD']-bc['CD'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD1'.format(vx1), bc['HD1']-(bc['HD1']/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vx1), bc['HD2']+((fc['HD2']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vx1), bc['HD3']+((fc['HD3']-bc['HD3'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vx1), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vx1), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def parmed_command_ADX(vx1):
	bc = {}
        with open('Param_files/AminoAcid/ABU.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/ADX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vx1), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vx1), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vx1), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vx1), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vx1), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vx1), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vx1), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vx1), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG1'.format(vx1), bc['HG1']-(bc['HG1']/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vx1), bc['HG2']-(bc['HG2']/10)*i).execute()
                change(parm, 'charge', ':{}@HG'.format(vx1), bc['HG3']+((fc['HG']-bc['HG3'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vx1), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vx1), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def parmed_command_ALX(vx1):
	bc = {}
        with open('Param_files/AminoAcid/ALA.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/ALX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vx1), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vx1), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vx1), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vx1), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vx1), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB1'.format(vx1), bc['HB1']-(bc['HB1']/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vx1), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vx1), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vx1), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vx1), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def parmed_command_I4X(vx1):
	bc = {}
        with open('Param_files/AminoAcid/NLE.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/NL4.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HA :{}@HB21 0 0'.format(vx1, vx1)).execute()
                change(parm, 'charge', ':{}@N'.format(vx1), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vx1), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vx1), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vx1), bc['HA']-(bc['HA']/10)*i).execute()
                change(parm, 'charge', ':{}@CB2'.format(vx1), (fc['CB2']/10)*i).execute()
                change(parm, 'charge', ':{}@HB21'.format(vx1), (fc['HB21']/10)*i).execute()
                change(parm, 'charge', ':{}@HB22'.format(vx1), (fc['HB22']/10)*i).execute()
                change(parm, 'charge', ':{}@HB23'.format(vx1), (fc['HB23']/10)*i).execute()
                change(parm, 'charge', ':{}@CB1'.format(vx1), bc['CB']+((fc['CB1']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB12'.format(vx1), bc['HB2']+((fc['HB12']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB13'.format(vx1), bc['HB3']+((fc['HB13']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vx1), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vx1), bc['HG2']+((fc['HG2']-bc['HG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vx1), bc['HG3']+((fc['HG3']-bc['HG3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vx1), bc['CD']+((fc['CD']-bc['CD'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vx1), bc['HD2']+((fc['HD2']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vx1), bc['HD3']+((fc['HD3']-bc['HD3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CE'.format(vx1), bc['CE']+((fc['CE']-bc['CE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE'.format(vx1), bc['HE1']+((fc['HE']-bc['HE1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE2'.format(vx1), bc['HE2']-(bc['HE2']/10)*i).execute()
                change(parm, 'charge', ':{}@HE3'.format(vx1), bc['HE3']-(bc['HE3']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vx1), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vx1), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def parmed_command_I7X(vx1):
	bc = {}
        with open('Param_files/AminoAcid/NKI.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/I7X.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
		changeLJPair(parm, ':{}@HA :{}@HB21 0 0'.format(vx1, vx1)).execute()
                change(parm, 'charge', ':{}@N'.format(vx1), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vx1), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vx1), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vx1), bc['HA']-(bc['HA']/10)*i).execute()
                change(parm, 'charge', ':{}@CB2'.format(vx1), (fc['CB2']/10)*i).execute()
                change(parm, 'charge', ':{}@HB21'.format(vx1), (fc['HB21']/10)*i).execute()
                change(parm, 'charge', ':{}@HB22'.format(vx1), (fc['HB22']/10)*i).execute()
                change(parm, 'charge', ':{}@HB23'.format(vx1), (fc['HB23']/10)*i).execute()
                change(parm, 'charge', ':{}@CB1'.format(vx1), bc['CB']+((fc['CB1']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB12'.format(vx1), bc['HB2']+((fc['HB12']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB13'.format(vx1), bc['HB3']+((fc['HB13']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vx1), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vx1), bc['HG2']+((fc['HG2']-bc['HG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vx1), bc['HG3']+((fc['HG3']-bc['HG3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vx1), bc['CD']+((fc['CD']-bc['CD'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vx1), bc['HD2']+((fc['HD2']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vx1), bc['HD3']+((fc['HD3']-bc['HD3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CE'.format(vx1), bc['CE']+((fc['CE']-bc['CE'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE2'.format(vx1), bc['HE2']+((fc['HE2']-bc['HE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE3'.format(vx1), bc['HE3']+((fc['HE3']-bc['HE3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CZ'.format(vx1), bc['CZ']+((fc['CZ']-bc['CZ'])/10)*i).execute()
                change(parm, 'charge', ':{}@HZ2'.format(vx1), bc['HZ2']+((fc['HZ2']-bc['HZ2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HZ3'.format(vx1), bc['HZ3']+((fc['HZ3']-bc['HZ3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CH'.format(vx1), bc['CH']+((fc['CH']-bc['CH'])/10)*i).execute()
                change(parm, 'charge', ':{}@HH2'.format(vx1), bc['HH2']+((fc['HH2']-bc['HH2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HH3'.format(vx1), bc['HH3']+((fc['HH3']-bc['HH3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CT'.format(vx1), bc['CT']+((fc['CT']-bc['CT'])/10)*i).execute()
                change(parm, 'charge', ':{}@HT'.format(vx1), bc['HT1']+((fc['HT']-bc['HT1'])/10)*i).execute()
                change(parm, 'charge', ':{}@HT2'.format(vx1), bc['HT2']-(bc['HT2']/10)*i).execute()
                change(parm, 'charge', ':{}@HT3'.format(vx1), bc['HT3']-(bc['HT3']/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vx1), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vx1), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def parmed_command_HSX(vx1):
	bc = {}
        with open('Param_files/AminoAcid/TZ2.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/HSX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vx1), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vx1), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vx1), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vx1), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vx1), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vx1), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vx1), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vx1), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vx1), bc['HG2']+((fc['HG2']-bc['HG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vx1), bc['HG3']+((fc['HG3']-bc['HG3'])/10)*i).execute()
                change(parm, 'charge', ':{}@ND'.format(vx1), bc['ND']+((fc['ND']-bc['ND'])/10)*i).execute()
                change(parm, 'charge', ':{}@NE1'.format(vx1), bc['NE1']+((fc['NE1']-bc['NE1'])/10)*i).execute()
                change(parm, 'charge', ':{}@NZ1'.format(vx1), bc['NZ1']+((fc['NZ1']-bc['NZ1'])/10)*i).execute()
                change(parm, 'charge', ':{}@CZ2'.format(vx1), bc['CZ2']+((fc['CZ2']-bc['CZ2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HZ2'.format(vx1), bc['HZ2']-(bc['HZ2']/10)*i).execute()
                change(parm, 'charge', ':{}@CE2'.format(vx1), bc['CE2']+((fc['CE2']-bc['CE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HE2'.format(vx1), bc['HE2']+((fc['HE2']-bc['HE2'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vx1), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vx1), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def parmed_command_VLX(vx1):
	bc = {}
        with open('Param_files/AminoAcid/NVA.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		bc[key] = float(value)
        b.close()
	fc = {}
        with open('Param_files/AminoAcid/VLX.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
		key, value = line.split()
		fc[key] = float(value)
        b.close()
	for i in range(11):
		a = i*10
		parm = AmberParm('Solv_{}_{}.prmtop'.format(a, 100-a))
                change(parm, 'charge', ':{}@N'.format(vx1), bc['N']+((fc['N']-bc['N'])/10)*i).execute()
                change(parm, 'charge', ':{}@H'.format(vx1), bc['H']+((fc['H']-bc['H'])/10)*i).execute()
                change(parm, 'charge', ':{}@CA'.format(vx1), bc['CA']+((fc['CA']-bc['CA'])/10)*i).execute()
                change(parm, 'charge', ':{}@HA'.format(vx1), bc['HA']+((fc['HA']-bc['HA'])/10)*i).execute()
                change(parm, 'charge', ':{}@CB'.format(vx1), bc['CB']+((fc['CB']-bc['CB'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB2'.format(vx1), bc['HB2']+((fc['HB2']-bc['HB2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HB3'.format(vx1), bc['HB3']+((fc['HB3']-bc['HB3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CG'.format(vx1), bc['CG']+((fc['CG']-bc['CG'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG2'.format(vx1), bc['HG2']+((fc['HG2']-bc['HG2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HG3'.format(vx1), bc['HG3']+((fc['HG3']-bc['HG3'])/10)*i).execute()
                change(parm, 'charge', ':{}@CD'.format(vx1), bc['CD']+((fc['CD']-bc['CD'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD1'.format(vx1), bc['HD1']-(bc['HD1']/10)*i).execute()
                change(parm, 'charge', ':{}@HD2'.format(vx1), bc['HD2']+((fc['HD2']-bc['HD2'])/10)*i).execute()
                change(parm, 'charge', ':{}@HD3'.format(vx1), bc['HD3']+((fc['HD3']-bc['HD3'])/10)*i).execute()
                change(parm, 'charge', ':{}@C'.format(vx1), bc['C']+((fc['C']-bc['C'])/10)*i).execute()
                change(parm, 'charge', ':{}@O'.format(vx1), bc['O']+((fc['O']-bc['O'])/10)*i).execute()
		setOverwrite(parm).execute()
		parmout(parm, 'Solv_{}_{}.prmtop'.format(a, 100-a)).execute()

def makevxi_HSX(struct, out, resid, vx1):
        struct.residue_dict[resid].set_resname(vx1)
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
                       	pdb.write(atom.formatted())
	                try:
        	                pdb.write(struct.other_dict[atom.get_number()].ter())
                	except:
                        	pass
        for oth in struct.other_dict:
                try:
                        if oth.startswith('Conect'):
                                pdb.write(struct.other_dict[oth].formatted())
                except:
                        pass
        pdb.write('END\n')

def makevxi_VLX(struct, out, resid, vx1):
        struct.residue_dict[resid].set_resname(vx1)
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
                       	pdb.write(atom.formatted())
	                try:
        	                pdb.write(struct.other_dict[atom.get_number()].ter())
                	except:
                        	pass
        for oth in struct.other_dict:
                try:
                        if oth.startswith('Conect'):
                                pdb.write(struct.other_dict[oth].formatted())
                except:
                        pass
        pdb.write('END\n')

def makevxi_ADX(struct, out, resid, vx1):
        struct.residue_dict[resid].set_resname(vx1)
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HG3' and res.get_resnumber() == resid:
                                pdb.write(atom.change_name('HG'))
			else:
	                       	pdb.write(atom.formatted())
	                try:
        	                pdb.write(struct.other_dict[atom.get_number()].ter())
                	except:
                        	pass
        for oth in struct.other_dict:
                try:
                        if oth.startswith('Conect'):
                                pdb.write(struct.other_dict[oth].formatted())
                except:
                        pass
        pdb.write('END\n')

def makevxi_MCX(struct, out, resid, vx1):
        struct.residue_dict[resid].set_resname(vx1)
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
                       	pdb.write(atom.formatted())
	                try:
        	                pdb.write(struct.other_dict[atom.get_number()].ter())
                	except:
                        	pass
        for oth in struct.other_dict:
                try:
                        if oth.startswith('Conect'):
                                pdb.write(struct.other_dict[oth].formatted())
                except:
                        pass
        pdb.write('END\n')

def makevxi_CYX(struct, out, resid, vx1):
        struct.residue_dict[resid].set_resname(vx1)
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
                       	pdb.write(atom.formatted())
	                try:
        	                pdb.write(struct.other_dict[atom.get_number()].ter())
                	except:
                        	pass
        for oth in struct.other_dict:
                try:
                        if oth.startswith('Conect'):
                                pdb.write(struct.other_dict[oth].formatted())
                except:
                        pass
        pdb.write('END\n')

def makevxi_ABX(struct, out, resid, vx1):
        struct.residue_dict[resid].set_resname(vx1)
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
                       	pdb.write(atom.formatted())
	                try:
        	                pdb.write(struct.other_dict[atom.get_number()].ter())
                	except:
                        	pass
        for oth in struct.other_dict:
                try:
                        if oth.startswith('Conect'):
                                pdb.write(struct.other_dict[oth].formatted())
                except:
                        pass
        pdb.write('END\n')

def makevxi_ALX(struct, out, resid, vx1):
        struct.residue_dict[resid].set_resname(vx1)
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
                       	pdb.write(atom.formatted())
	                try:
        	                pdb.write(struct.other_dict[atom.get_number()].ter())
                	except:
                        	pass
        for oth in struct.other_dict:
                try:
                        if oth.startswith('Conect'):
                                pdb.write(struct.other_dict[oth].formatted())
                except:
                        pass
        pdb.write('END\n')

def makevxi_AMX(struct, out, resid, vx1):
        struct.residue_dict[resid].set_resname(vx1)
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HD21' and res.get_resnumber() == resid:
                                pdb.write(atom.change_name('HD2'))
			else:
                        	pdb.write(atom.formatted())
	                try:
        	                pdb.write(struct.other_dict[atom.get_number()].ter())
                	except:
                        	pass
        for oth in struct.other_dict:
                try:
                        if oth.startswith('Conect'):
                                pdb.write(struct.other_dict[oth].formatted())
                except:
                        pass
        pdb.write('END\n')

def makevxi_I4X(struct, out, resid, vx1):
        struct.residue_dict[resid].set_resname(vx1)
        CA = struct.residue_dict[resid].atom_dict['CA']
        HA = struct.residue_dict[resid].atom_dict['HA']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HA' and res.get_resnumber() == resid:
                        	pdb.write(atom.formatted())
				pdb.write(atom.halfway_between('CB2', CA, HA))
                                pdb.write(atom.superimposed1('HB21', HA))
                                pdb.write(atom.superimposed2('HB22', HA))
                                pdb.write(atom.superimposed3('HB23', HA))
			elif atom.get_name() == 'CB' and res.get_resname() == vx1:
                                pdb.write(atom.change_name('CB1'))
			elif atom.get_name() == 'HB2' and res.get_resname() == vx1:
                                pdb.write(atom.change_name('HB12'))
			elif atom.get_name() == 'HB3' and res.get_resname() == vx1:
                                pdb.write(atom.change_name('HB13'))
			elif atom.get_name() == 'HE1' and res.get_resname() == vx1:
                                pdb.write(atom.change_name('HE'))
			else:
                        	pdb.write(atom.formatted())
	                try:
        	                pdb.write(struct.other_dict[atom.get_number()].ter())
                	except:
                        	pass
        for oth in struct.other_dict:
                try:
                        if oth.startswith('Conect'):
                                pdb.write(struct.other_dict[oth].formatted())
                except:
                        pass
        pdb.write('END\n')

def makevxi_I7X(struct, out, resid, vx1):
        struct.residue_dict[resid].set_resname(vx1)
        CA = struct.residue_dict[resid].atom_dict['CA']
        HA = struct.residue_dict[resid].atom_dict['HA']
	pdb = open(out, 'w')
        try:
                pdb.write(struct.other_dict['Cryst1'].formatted())
        except KeyError:
                pass
        for res in struct.residue_list:
                for atom in res.atom_list:
			if atom.get_name() == 'HA' and res.get_resnumber() == resid:
                        	pdb.write(atom.formatted())
				pdb.write(atom.halfway_between('CB2', CA, HA))
                                pdb.write(atom.superimposed1('HB21', HA))
                                pdb.write(atom.superimposed2('HB22', HA))
                                pdb.write(atom.superimposed3('HB23', HA))
			elif atom.get_name() == 'CB' and res.get_resname() in (vx1):
                                pdb.write(atom.change_name('CB1'))
			elif atom.get_name() == 'HB2' and res.get_resname() in (vx1):
                                pdb.write(atom.change_name('HB12'))
			elif atom.get_name() == 'HB3' and res.get_resname() in (vx1):
                                pdb.write(atom.change_name('HB13'))
			elif atom.get_name() == 'HT1' and res.get_resname() == vx1:
                                pdb.write(atom.change_name('HT'))
			else:
                        	pdb.write(atom.formatted())
	                try:
        	                pdb.write(struct.other_dict[atom.get_number()].ter())
                	except:
                        	pass
        for oth in struct.other_dict:
                try:
                        if oth.startswith('Conect'):
                                pdb.write(struct.other_dict[oth].formatted())
                except:
                        pass
        pdb.write('END\n')

def lib_make_VLX(ff, outputfile, vxi, newcar='dc', gonhyd='dh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/VLX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "C"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "C"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "C"\n'%vxi)
	ctrl.write('set %s.1.16 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG3"\n'%vxi)
	ctrl.write('set %s.1.11 name "CD"\n'%vxi)
	ctrl.write('set %s.1.12 name "HD1"\n'%vxi)
	ctrl.write('set %s.1.13 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.14 name "HD3"\n'%vxi)
	ctrl.write('set %s.1.15 name "C"\n'%vxi)
	ctrl.write('set %s.1.16 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "CT"\n'%vxi)
	ctrl.write('set %s.1.9 type "HC"\n'%vxi)
	ctrl.write('set %s.1.10 type "HC"\n'%vxi)
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, newcar))
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.13 type "HC"\n'%vxi)
	ctrl.write('set %s.1.14 type "HC"\n'%vxi)
	ctrl.write('set %s.1.15 type "C"\n'%vxi)
	ctrl.write('set %s.1.16 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.16\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.11\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def lib_make_HSX(ff, outputfile, vxi, gonhyd='gh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/HSX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "C"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "N"\n'%vxi)
	ctrl.write('set %s.1.12 element "N"\n'%vxi)
	ctrl.write('set %s.1.13 element "N"\n'%vxi)
	ctrl.write('set %s.1.14 element "C"\n'%vxi)
	ctrl.write('set %s.1.15 element "H"\n'%vxi)
	ctrl.write('set %s.1.16 element "C"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "C"\n'%vxi)
	ctrl.write('set %s.1.19 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG3"\n'%vxi)
	ctrl.write('set %s.1.11 name "ND"\n'%vxi)
	ctrl.write('set %s.1.12 name "NE1"\n'%vxi)
	ctrl.write('set %s.1.13 name "NZ1"\n'%vxi)
	ctrl.write('set %s.1.14 name "CZ2"\n'%vxi)
	ctrl.write('set %s.1.15 name "HZ2"\n'%vxi)
	ctrl.write('set %s.1.16 name "CE2"\n'%vxi)
	ctrl.write('set %s.1.17 name "HE2"\n'%vxi)
	ctrl.write('set %s.1.18 name "C"\n'%vxi)
	ctrl.write('set %s.1.19 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "CT"\n'%vxi)
	ctrl.write('set %s.1.9 type "H1"\n'%vxi)
	ctrl.write('set %s.1.10 type "H1"\n'%vxi)
	ctrl.write('set %s.1.11 type "NX"\n'%vxi)
	ctrl.write('set %s.1.12 type "NW"\n'%vxi)
	ctrl.write('set %s.1.13 type "NZ"\n'%vxi)
	ctrl.write('set %s.1.14 type "CB"\n'%vxi)
	ctrl.write('set %s.1.15 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.16 type "CC"\n'%vxi)
	ctrl.write('set %s.1.17 type "H4"\n'%vxi)
	ctrl.write('set %s.1.18 type "C"\n'%vxi)
	ctrl.write('set %s.1.19 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.11 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.14 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.16 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.19\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.11\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def lib_make_ADX(ff, outputfile, vxi, newcar='dc', newhyd='dh', gonhyd='gh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ADX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "C"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "C"\n'%vxi)
	ctrl.write('set %s.1.13 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG1"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.11 name "HG"\n'%vxi)
	ctrl.write('set %s.1.12 name "C"\n'%vxi)
	ctrl.write('set %s.1.13 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, newcar))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.11 type "%s"\n'%(vxi, newhyd))
	ctrl.write('set %s.1.12 type "C"\n'%vxi)
	ctrl.write('set %s.1.13 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.8\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def lib_make_MCX(ff, outputfile, vxi, newcar='dc', gonhyd='gh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/MCX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "S"\n'%vxi)
	ctrl.write('set %s.1.9 element "C"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "C"\n'%vxi)
	ctrl.write('set %s.1.14 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "SG"\n'%vxi)
	ctrl.write('set %s.1.9 name "CD"\n'%vxi)
	ctrl.write('set %s.1.10 name "HD1"\n'%vxi)
	ctrl.write('set %s.1.11 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.12 name "HD3"\n'%vxi)
	ctrl.write('set %s.1.13 name "C"\n'%vxi)
	ctrl.write('set %s.1.14 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "H1"\n'%vxi)
	ctrl.write('set %s.1.7 type "H1"\n'%vxi)
	ctrl.write('set %s.1.8 type "S"\n'%vxi)
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, newcar))
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.11 type "H1"\n'%vxi)
	ctrl.write('set %s.1.12 type "H1"\n'%vxi)
	ctrl.write('set %s.1.13 type "C"\n'%vxi)
	ctrl.write('set %s.1.14 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.9\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def lib_make_CYX(ff, outputfile, vxi, newsul='ds', gonhyd='dh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/CYX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "S"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "C"\n'%vxi)
	ctrl.write('set %s.1.11 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "SG"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG"\n'%vxi)
	ctrl.write('set %s.1.10 name "C"\n'%vxi)
	ctrl.write('set %s.1.11 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "H1"\n'%vxi)
	ctrl.write('set %s.1.7 type "H1"\n'%vxi)
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, newsul))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.10 type "C"\n'%vxi)
	ctrl.write('set %s.1.11 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.11\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.8\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def lib_make_ABX(ff, outputfile, vxi, newcar='dc', gonhyd='dh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ABX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "C"\n'%vxi)
	ctrl.write('set %s.1.9 element "H"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "C"\n'%vxi)
	ctrl.write('set %s.1.13 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG"\n'%vxi)
	ctrl.write('set %s.1.9 name "HG1"\n'%vxi)
	ctrl.write('set %s.1.10 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.11 name "HG3"\n'%vxi)
	ctrl.write('set %s.1.12 name "C"\n'%vxi)
	ctrl.write('set %s.1.13 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, newcar))
	ctrl.write('set %s.1.9 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.10 type "HC"\n'%vxi)
	ctrl.write('set %s.1.11 type "HC"\n'%vxi)
	ctrl.write('set %s.1.12 type "C"\n'%vxi)
	ctrl.write('set %s.1.13 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.8\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def lib_make_AMX(ff, outputfile, vxi, newnit='gn', gonhyd='gh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/AMX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "C"\n'%vxi)
	ctrl.write('set %s.1.9 element "O"\n'%vxi)
	ctrl.write('set %s.1.10 element "N"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "H"\n'%vxi)
	ctrl.write('set %s.1.13 element "C"\n'%vxi)
	ctrl.write('set %s.1.14 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.8 name "CG"\n'%vxi)
	ctrl.write('set %s.1.9 name "OD1"\n'%vxi)
	ctrl.write('set %s.1.10 name "ND2"\n'%vxi)
	ctrl.write('set %s.1.11 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.12 name "HD22"\n'%vxi)
	ctrl.write('set %s.1.13 name "C"\n'%vxi)
	ctrl.write('set %s.1.14 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "CT"\n'%vxi)
	ctrl.write('set %s.1.6 type "HC"\n'%vxi)
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "C"\n'%vxi)
	ctrl.write('set %s.1.9 type "O"\n'%vxi)
	ctrl.write('set %s.1.10 type "%s"\n'%(vxi, newnit))
	ctrl.write('set %s.1.11 type "H"\n'%vxi)
	ctrl.write('set %s.1.12 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.13 type "C"\n'%vxi)
	ctrl.write('set %s.1.14 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.8 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.10 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.13 %s.1.14\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect2 %s.1.10\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def lib_make_ALX(ff, outputfile, vxi, newcar='dc', gonhyd='dh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/ALX.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "H"\n'%vxi)
	ctrl.write('set %s.1.9 element "C"\n'%vxi)
	ctrl.write('set %s.1.10 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB1"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB2"\n'%vxi)
	ctrl.write('set %s.1.8 name "HB3"\n'%vxi)
	ctrl.write('set %s.1.9 name "C"\n'%vxi)
	ctrl.write('set %s.1.10 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "H1"\n'%vxi)
	ctrl.write('set %s.1.5 type "%s"\n'%(vxi, newcar))
	ctrl.write('set %s.1.6 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.7 type "HC"\n'%vxi)
	ctrl.write('set %s.1.8 type "HC"\n'%vxi)
	ctrl.write('set %s.1.9 type "C"\n'%vxi)
	ctrl.write('set %s.1.10 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.5\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def lib_make_I4X(ff, outputfile, vxi, newcar='dc', newhyd='dh', gonhyd='gh', metcar='mc', methyd='mh', hydhyd='hh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/I4X.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "H"\n'%vxi)
	ctrl.write('set %s.1.9 element "C"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "C"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "C"\n'%vxi)
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "C"\n'%vxi)
	ctrl.write('set %s.1.19 element "H"\n'%vxi)
	ctrl.write('set %s.1.20 element "H"\n'%vxi)
	ctrl.write('set %s.1.21 element "H"\n'%vxi)
	ctrl.write('set %s.1.22 element "C"\n'%vxi)
	ctrl.write('set %s.1.23 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB2"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB21"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB22"\n'%vxi)
	ctrl.write('set %s.1.8 name "HB23"\n'%vxi)
	ctrl.write('set %s.1.9 name "CB1"\n'%vxi)
	ctrl.write('set %s.1.10 name "HB12"\n'%vxi)
	ctrl.write('set %s.1.11 name "HB13"\n'%vxi)
	ctrl.write('set %s.1.12 name "CG"\n'%vxi)
	ctrl.write('set %s.1.13 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.14 name "HG3"\n'%vxi)
	ctrl.write('set %s.1.15 name "CD"\n'%vxi)
	ctrl.write('set %s.1.16 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.17 name "HD3"\n'%vxi)
	ctrl.write('set %s.1.18 name "CE"\n'%vxi)
	ctrl.write('set %s.1.19 name "HE"\n'%vxi)
	ctrl.write('set %s.1.20 name "HE2"\n'%vxi)
	ctrl.write('set %s.1.21 name "HE3"\n'%vxi)
	ctrl.write('set %s.1.22 name "C"\n'%vxi)
	ctrl.write('set %s.1.23 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "%s"\n'%(vxi, hydhyd))
	ctrl.write('set %s.1.5 type "%s"\n'%(vxi, metcar))
	ctrl.write('set %s.1.6 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.7 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.9 type "CT"\n'%vxi)
	ctrl.write('set %s.1.10 type "HC"\n'%vxi)
	ctrl.write('set %s.1.11 type "HC"\n'%vxi)
	ctrl.write('set %s.1.12 type "CT"\n'%vxi)
	ctrl.write('set %s.1.13 type "HC"\n'%vxi)
	ctrl.write('set %s.1.14 type "HC"\n'%vxi)
	ctrl.write('set %s.1.15 type "CT"\n'%vxi)
	ctrl.write('set %s.1.16 type "HC"\n'%vxi)
	ctrl.write('set %s.1.17 type "HC"\n'%vxi)
	ctrl.write('set %s.1.18 type "%s"\n'%(vxi, newcar))
	ctrl.write('set %s.1.19 type "%s"\n'%(vxi, newhyd))
	ctrl.write('set %s.1.20 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.21 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.22 type "C"\n'%vxi)
	ctrl.write('set %s.1.23 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.22 %s.1.23\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.18\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def lib_make_I7X(ff, outputfile, vxi, newcar='dc', newhyd='dh', gonhyd='gh', metcar='mc', methyd='mh', hydhyd='hh'):
        ctrl = open('lyp.in', 'w')
        ctrl.write("source %s\n"%ff)
	ctrl.write("%s=loadpdb Param_files/LibPDB/I7X.pdb\n"%vxi)
	ctrl.write('set %s.1.1 element "N"\n'%vxi)
	ctrl.write('set %s.1.2 element "H"\n'%vxi)
	ctrl.write('set %s.1.3 element "C"\n'%vxi)
	ctrl.write('set %s.1.4 element "H"\n'%vxi)
	ctrl.write('set %s.1.5 element "C"\n'%vxi)
	ctrl.write('set %s.1.6 element "H"\n'%vxi)
	ctrl.write('set %s.1.7 element "H"\n'%vxi)
	ctrl.write('set %s.1.8 element "H"\n'%vxi)
	ctrl.write('set %s.1.9 element "C"\n'%vxi)
	ctrl.write('set %s.1.10 element "H"\n'%vxi)
	ctrl.write('set %s.1.11 element "H"\n'%vxi)
	ctrl.write('set %s.1.12 element "C"\n'%vxi)
	ctrl.write('set %s.1.13 element "H"\n'%vxi)
	ctrl.write('set %s.1.14 element "H"\n'%vxi)
	ctrl.write('set %s.1.15 element "C"\n'%vxi)
	ctrl.write('set %s.1.16 element "H"\n'%vxi)
	ctrl.write('set %s.1.17 element "H"\n'%vxi)
	ctrl.write('set %s.1.18 element "C"\n'%vxi)
	ctrl.write('set %s.1.19 element "H"\n'%vxi)
	ctrl.write('set %s.1.20 element "H"\n'%vxi)
	ctrl.write('set %s.1.21 element "C"\n'%vxi)
	ctrl.write('set %s.1.22 element "H"\n'%vxi)
	ctrl.write('set %s.1.23 element "H"\n'%vxi)
	ctrl.write('set %s.1.24 element "C"\n'%vxi)
	ctrl.write('set %s.1.25 element "H"\n'%vxi)
	ctrl.write('set %s.1.26 element "H"\n'%vxi)
	ctrl.write('set %s.1.27 element "C"\n'%vxi)
	ctrl.write('set %s.1.28 element "H"\n'%vxi)
	ctrl.write('set %s.1.29 element "H"\n'%vxi)
	ctrl.write('set %s.1.30 element "H"\n'%vxi)
	ctrl.write('set %s.1.31 element "C"\n'%vxi)
	ctrl.write('set %s.1.32 element "O"\n'%vxi)
	ctrl.write('set %s.1.1 name "N"\n'%vxi)
	ctrl.write('set %s.1.2 name "H"\n'%vxi)
	ctrl.write('set %s.1.3 name "CA"\n'%vxi)
	ctrl.write('set %s.1.4 name "HA"\n'%vxi)
	ctrl.write('set %s.1.5 name "CB2"\n'%vxi)
	ctrl.write('set %s.1.6 name "HB21"\n'%vxi)
	ctrl.write('set %s.1.7 name "HB22"\n'%vxi)
	ctrl.write('set %s.1.8 name "HB23"\n'%vxi)
	ctrl.write('set %s.1.9 name "CB1"\n'%vxi)
	ctrl.write('set %s.1.10 name "HB12"\n'%vxi)
	ctrl.write('set %s.1.11 name "HB13"\n'%vxi)
	ctrl.write('set %s.1.12 name "CG"\n'%vxi)
	ctrl.write('set %s.1.13 name "HG2"\n'%vxi)
	ctrl.write('set %s.1.14 name "HG3"\n'%vxi)
	ctrl.write('set %s.1.15 name "CD"\n'%vxi)
	ctrl.write('set %s.1.16 name "HD2"\n'%vxi)
	ctrl.write('set %s.1.17 name "HD3"\n'%vxi)
	ctrl.write('set %s.1.18 name "CE"\n'%vxi)
	ctrl.write('set %s.1.19 name "HE2"\n'%vxi)
	ctrl.write('set %s.1.20 name "HE3"\n'%vxi)
	ctrl.write('set %s.1.21 name "CZ"\n'%vxi)
	ctrl.write('set %s.1.22 name "HZ2"\n'%vxi)
	ctrl.write('set %s.1.23 name "HZ3"\n'%vxi)
	ctrl.write('set %s.1.24 name "CH"\n'%vxi)
	ctrl.write('set %s.1.25 name "HH2"\n'%vxi)
	ctrl.write('set %s.1.26 name "HH3"\n'%vxi)
	ctrl.write('set %s.1.27 name "CT"\n'%vxi)
	ctrl.write('set %s.1.28 name "HT"\n'%vxi)
	ctrl.write('set %s.1.29 name "HT2"\n'%vxi)
	ctrl.write('set %s.1.30 name "HT3"\n'%vxi)
	ctrl.write('set %s.1.31 name "C"\n'%vxi)
	ctrl.write('set %s.1.32 name "O"\n'%vxi)
	ctrl.write('set %s.1.1 type "N"\n'%vxi)
	ctrl.write('set %s.1.2 type "H"\n'%vxi)
	ctrl.write('set %s.1.3 type "CT"\n'%vxi)
	ctrl.write('set %s.1.4 type "%s"\n'%(vxi, hydhyd))
	ctrl.write('set %s.1.5 type "%s"\n'%(vxi, metcar))
	ctrl.write('set %s.1.6 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.7 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.8 type "%s"\n'%(vxi, methyd))
	ctrl.write('set %s.1.9 type "CT"\n'%vxi)
	ctrl.write('set %s.1.10 type "HC"\n'%vxi)
	ctrl.write('set %s.1.11 type "HC"\n'%vxi)
	ctrl.write('set %s.1.12 type "CT"\n'%vxi)
	ctrl.write('set %s.1.13 type "HC"\n'%vxi)
	ctrl.write('set %s.1.14 type "HC"\n'%vxi)
	ctrl.write('set %s.1.15 type "CT"\n'%vxi)
	ctrl.write('set %s.1.16 type "HC"\n'%vxi)
	ctrl.write('set %s.1.17 type "HC"\n'%vxi)
	ctrl.write('set %s.1.18 type "CT"\n'%vxi)
	ctrl.write('set %s.1.19 type "HC"\n'%vxi)
	ctrl.write('set %s.1.20 type "HC"\n'%vxi)
	ctrl.write('set %s.1.21 type "CT"\n'%vxi)
	ctrl.write('set %s.1.22 type "HC"\n'%vxi)
	ctrl.write('set %s.1.23 type "HC"\n'%vxi)
	ctrl.write('set %s.1.24 type "CT"\n'%vxi)
	ctrl.write('set %s.1.25 type "HC"\n'%vxi)
	ctrl.write('set %s.1.26 type "HC"\n'%vxi)
	ctrl.write('set %s.1.27 type "%s"\n'%(vxi, newcar))
	ctrl.write('set %s.1.28 type "%s"\n'%(vxi, newhyd))
	ctrl.write('set %s.1.29 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.30 type "%s"\n'%(vxi, gonhyd))
	ctrl.write('set %s.1.31 type "C"\n'%vxi)
	ctrl.write('set %s.1.32 type "O"\n'%vxi)
	ctrl.write('bond %s.1.1 %s.1.2\n'%(vxi, vxi))
	ctrl.write('bond %s.1.1 %s.1.3\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.4\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.5\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.9\n'%(vxi, vxi))
	ctrl.write('bond %s.1.3 %s.1.31\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.6\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.7\n'%(vxi, vxi))
	ctrl.write('bond %s.1.5 %s.1.8\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.10\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.11\n'%(vxi, vxi))
	ctrl.write('bond %s.1.9 %s.1.12\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.13\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.14\n'%(vxi, vxi))
	ctrl.write('bond %s.1.12 %s.1.15\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.16\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.17\n'%(vxi, vxi))
	ctrl.write('bond %s.1.15 %s.1.18\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.19\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.20\n'%(vxi, vxi))
	ctrl.write('bond %s.1.18 %s.1.21\n'%(vxi, vxi))
	ctrl.write('bond %s.1.21 %s.1.22\n'%(vxi, vxi))
	ctrl.write('bond %s.1.21 %s.1.23\n'%(vxi, vxi))
	ctrl.write('bond %s.1.21 %s.1.24\n'%(vxi, vxi))
	ctrl.write('bond %s.1.24 %s.1.25\n'%(vxi, vxi))
	ctrl.write('bond %s.1.24 %s.1.26\n'%(vxi, vxi))
	ctrl.write('bond %s.1.24 %s.1.27\n'%(vxi, vxi))
	ctrl.write('bond %s.1.27 %s.1.28\n'%(vxi, vxi))
	ctrl.write('bond %s.1.27 %s.1.29\n'%(vxi, vxi))
	ctrl.write('bond %s.1.27 %s.1.30\n'%(vxi, vxi))
	ctrl.write('bond %s.1.31 %s.1.32\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect0 %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.C\n'%(vxi, vxi))
	ctrl.write('set %s.1 connect1 %s.1.27\n'%(vxi, vxi))
	ctrl.write('set %s name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s.1 name "%s"\n'%(vxi, vxi))
	ctrl.write('set %s head %s.1.N\n'%(vxi, vxi))
	ctrl.write('set %s tail %s.1.C\n'%(vxi, vxi))
	ctrl.write('saveoff %s %s.lib\n'%(vxi, vxi))
	ctrl.write("quit\n") 
        ctrl.close()
	Leapy.run('lyp.in', outputfile)

def all_make():
	for i in range(0,110,10):
		Frcmod_creator.make ('{}_{}.frcmod'.format(i, 100-i))

def cal(x, y, i):
        num = x+((y-x)/10)*i
        return num

def alc(x, y, i):
        num = x+((y-x)/10)*i*i/10
        return num

def lac(y, x, i):
        num = x+((y-x)/10)*i
        return num

def stock_add_to_all_AMX(dist, newnit='gn', gonhyd='gh'):
        p = {}
        with open('Param_files/Stock/Stock.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
                p[line.split()[0]] = []
                for point in line.split()[1:]:
                        p[line.split()[0]].append(float(point))
        b.close()
        for i in range(11):
                a = i*10
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newnit, cal(p['NA'][0], p['NA'][0], i), cal(p['NA'][1], p['NA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['H'][0], p['0_H'][0], i), cal(p['H'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newnit, gonhyd), cal(p['NA_H'][0], p['H_sCT'][0], i), cal(p['NA_H'][1], p['H_sCT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newnit, 'H'), cal(p['NA_H'][0], p['NA_H'][0], i), cal(p['NA_H'][1], p['NA_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newnit, 'dc'), cal(p['0_0'][0], p['CT_N'][0], i), cal(dist, p['CT_N'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newnit, 'C '), cal(p['C_N'][0], p['C_N'][0], i), cal(p['C_N'][1], p['C_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(gonhyd, newnit, 'dc'), cal(p['X_CM_CM_0'][0], p['Close'][0], i), cal(p['X_CM_CM_0'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H', newnit, 'dc'), cal(p['X_CM_CM_0'][0], p['C_N_H'][0], i), cal(p['X_CM_CM_0'][1], p['C_N_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C ', newnit, 'dc'), cal(p['X_CM_CM_0'][0], p['C_N_CT'][0], i), cal(p['X_CM_CM_0'][1], p['C_N_CT'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newnit, 'dc', 'dh'), cal(p['X_CM_CM_0'][0], p['Close'][0], i), cal(p['X_CM_CM_0'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newnit, 'dc', 'HC'), cal(p['X_CM_CM_0'][0], p['C_C_H'][0], i), cal(p['X_CM_CM_0'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newnit, 'dc', 'CT'), cal(p['X_CM_CM_0'][0], p['CT_CT_N'][0], i), cal(p['X_CM_CM_0'][1], p['CT_CT_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H ', newnit, gonhyd), cal(p['H_N_H'][0], p['X_CM_CM_0'][0], i), cal(p['H_N_H'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C ', newnit, gonhyd), cal(p['C_N_H'][0], p['X_CM_CM_0'][0], i), cal(p['C_N_H'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C ', newnit, 'H '), cal(p['C_N_H'][0], p['C_N_H'][0], i), cal(p['C_N_H'][1], p['C_N_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'C ', newnit), cal(p['C_C_N'][0], p['C_C_N'][0], i), cal(p['C_C_N'][1], p['C_C_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('O ', 'C ', newnit), cal(p['O_C_N'][0], p['O_C_N'][0], i), cal(p['O_C_N'][1], p['O_C_N'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newnit, 'dc', 'dh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newnit, 'dc', 'CT'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newnit, 'dc', 'HC'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H ', newnit, 'dc', 'dh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', newnit, 'dc', 'dh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', newnit, 'dc', 'HC'), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', newnit, 'dc', 'CT'), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H ', newnit, 'dc', 'HC'), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H ', newnit, 'dc', 'CT'), cal(p['Ring_Dihe_2'][0], p['Ring_Dihe_2'][0], i), cal(p['Ring_Dihe_2'][1], p['Ring_Dihe_2'][1], i), cal(p['Ring_Dihe_2'][2], p['Ring_Dihe_2'][2], i), cal(p['Ring_Dihe_2'][3], p['Ring_Dihe_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H ', newnit, 'C ', 'O '), cal(p['O_C_N_H_1'][0], p['O_C_N_H_1'][0], i), cal(p['O_C_N_H_1'][1], p['O_C_N_H_1'][1], i), cal(p['O_C_N_H_1'][2], p['O_C_N_H_1'][2], i), cal(p['O_C_N_H_1'][3], p['O_C_N_H_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H ', newnit, 'C ', 'O '), cal(p['O_C_N_H_2'][0], p['O_C_N_H_2'][0], i), cal(p['O_C_N_H_2'][1], p['O_C_N_H_2'][1], i), cal(p['O_C_N_H_2'][2], p['O_C_N_H_2'][2], i), cal(p['O_C_N_H_2'][3], p['O_C_N_H_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newnit, 'C ', 'O '), cal(p['O_C_N_H_1'][0], p['0_3'][0], i), cal(p['O_C_N_H_1'][1], p['0_3'][1], i), cal(p['O_C_N_H_1'][2], p['0_3'][2], i), cal(p['O_C_N_H_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newnit, 'C ', 'O '), cal(p['O_C_N_H_2'][0], p['0_11'][0], i), cal(p['O_C_N_H_2'][1], p['0_11'][1], i), cal(p['O_C_N_H_2'][2], p['0_11'][2], i), cal(p['O_C_N_H_2'][3], p['0_11'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('dc', newnit, 'C ', 'O '), cal(p['Ring_0'][0], p['X_C_N_X'][0], i), cal(p['Ring_0'][1], p['X_C_N_X'][1], i), cal(p['Ring_0'][2], p['X_C_N_X'][2], i), cal(p['Ring_0'][3], p['X_C_N_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('dc', newnit, 'C ', 'CT'), cal(p['Ring_0'][0], p['X_C_N_X'][0], i), cal(p['Ring_0'][1], p['X_C_N_X'][1], i), cal(p['Ring_0'][2], p['X_C_N_X'][2], i), cal(p['Ring_0'][3], p['X_C_N_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newnit, 'C ', 'CT'), cal(p['X_C_N_X'][0], p['Ring_0'][0], i), cal(p['X_C_N_X'][1], p['Ring_0'][1], i), cal(p['X_C_N_X'][2], p['Ring_0'][2], i), cal(p['X_C_N_X'][3], p['Ring_0'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H ', newnit, 'C ', 'CT'), cal(p['X_C_N_X'][0], p['X_C_N_X'][0], i), cal(p['X_C_N_X'][1], p['X_C_N_X'][1], i), cal(p['X_C_N_X'][2], p['X_C_N_X'][2], i), cal(p['X_C_N_X'][3], p['X_C_N_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newnit, 'dc', 'CT', 'CT'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newnit, 'dc', 'CT', 'HC'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newnit, 'dc', 'CT', 'H1'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newnit, 'dc', 'CT', 'C '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newnit, 'dc', 'CT', 'N'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'H ', newnit, gonhyd), cal(p['Ami_imp'][0], p['Imp_0'][0], i), cal(p['Ami_imp'][1], p['Imp_0'][1], i), cal(p['Ami_imp'][2], p['Imp_0'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'H ', newnit, 'dc'), cal(p['Imp_0'][0], p['Ami_imp'][0], i), cal(p['Imp_0'][1], p['Ami_imp'][1], i), cal(p['Imp_0'][2], p['Ami_imp'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', gonhyd, newnit, 'dc'), cal(p['Imp_0'][0], p['Imp_0'][0], i), cal(p['Imp_0'][1], p['Imp_0'][1], i), cal(p['Imp_0'][2], p['Imp_0'][2], i))
		Frcmod_creator.IMPROPER_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H ', gonhyd, newnit, 'dc'), cal(p['Imp_0'][0], p['Imp_0'][0], i), cal(p['Imp_0'][1], p['Imp_0'][1], i), cal(p['Imp_0'][2], p['Imp_0'][2], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newnit, cal(p['NA'][2], p['NA'][2], i), cal(p['NA'][3], p['NA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['H'][2], p['0_H'][2], i), cal(p['H'][3], p['0_H'][3], i))

def stock_add_to_all_HSX(dist, gonhyd='gh'):
        p = {}
        with open('Param_files/Stock/Stock.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
                p[line.split()[0]] = []
                for point in line.split()[1:]:
                        p[line.split()[0]].append(float(point))
        b.close()
        for i in range(11):
                a = i*10
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['H4'][0], p['0_H'][0], i), cal(p['H4'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CB', gonhyd), cal(p['CB_H4'][0], p['H_sCC'][0], i), cal(p['CB_H4'][1], p['H_sCC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CB', 'dc'), cal(p['0_0'][0], p['CB_CT'][0], i), cal(p['0_0'][1], p['CB_CT'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CC', 'CB', gonhyd), cal(p['CC_CB_H4'][0], p['X_CM_CM_0'][0], i), cal(p['CC_CB_H4'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('NZ', 'CB', gonhyd), cal(p['NZ_CB_H4'][0], p['X_CM_CM_0'][0], i), cal(p['NZ_CB_H4'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CB', 'dc', 'CT'), cal(p['X_CM_CM_0'][0], p['CB_CT_CT'][0], i), cal(p['X_CM_CM_0'][1], p['CB_CT_CT'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CB', 'dc', 'HC'), cal(p['X_CM_CM_0'][0], p['CB_CT_HC'][0], i), cal(p['X_CM_CM_0'][1], p['CB_CT_HC'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CB', 'dc', 'dh'), cal(p['X_CM_CM_0'][0], p['Close'][0], i), cal(p['X_CM_CM_0'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(gonhyd, 'CB', 'dc'), cal(p['X_CM_CM_0'][0], p['Close'][0], i), cal(p['X_CM_CM_0'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('NZ', 'CB', 'dc'), cal(p['X_CM_CM_0'][0], p['NZ_CB_CT'][0], i), cal(p['X_CM_CM_0'][1], p['NZ_CB_CT'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CC', 'CB', 'dc'), cal(p['X_CM_CM_0'][0], p['CC_CB_CT'][0], i), cal(p['X_CM_CM_0'][1], p['CC_CB_CT'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, 'CB', 'CC', 'NX'), cal(p['X_CB_CC_X'][0], p['0_1'][0], i), cal(p['X_CB_CC_X'][1], p['0_1'][1], i), cal(p['X_CB_CC_X'][2], p['0_1'][2], i), cal(p['X_CB_CC_X'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CC', 'CB', 'dc', 'HC'), cal(p['X_CT_N_X'][0], p['X_CT_N_X'][0], i), cal(p['X_CT_N_X'][1], p['X_CT_N_X'][1], i), cal(p['X_CT_N_X'][2], p['X_CT_N_X'][2], i), cal(p['X_CT_N_X'][3], p['X_CT_N_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CC', 'CB', 'dc', 'CT'), cal(p['X_CT_N_X'][0], p['X_CT_N_X'][0], i), cal(p['X_CT_N_X'][1], p['X_CT_N_X'][1], i), cal(p['X_CT_N_X'][2], p['X_CT_N_X'][2], i), cal(p['X_CT_N_X'][3], p['X_CT_N_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CC', 'CB', 'dc', 'dh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('NZ', 'CB', 'dc', 'HC'), cal(p['X_CT_N_X'][0], p['X_CT_N_X'][0], i), cal(p['X_CT_N_X'][1], p['X_CT_N_X'][1], i), cal(p['X_CT_N_X'][2], p['X_CT_N_X'][2], i), cal(p['X_CT_N_X'][3], p['X_CT_N_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('NZ', 'CB', 'dc', 'CT'), cal(p['X_CT_N_X'][0], p['X_CT_N_X'][0], i), cal(p['X_CT_N_X'][1], p['X_CT_N_X'][1], i), cal(p['X_CT_N_X'][2], p['X_CT_N_X'][2], i), cal(p['X_CT_N_X'][3], p['X_CT_N_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('NZ', 'CB', 'dc', 'dh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, 'CB', 'dc', 'dh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, 'CB', 'dc', 'CT'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, 'CB', 'dc', 'HC'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CB', 'dc', 'CT', 'CT'), cal(p['0_12'][0], p['C_C_C_C_1'][0], i), cal(p['0_12'][1], p['C_C_C_C_1'][1], i), cal(p['0_12'][2], p['C_C_C_C_1'][2], i), cal(p['0_12'][3], p['C_C_C_C_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CB', 'dc', 'CT', 'CT'), cal(p['0_11'][0], p['C_C_C_C_2'][0], i), cal(p['0_11'][1], p['C_C_C_C_2'][1], i), cal(p['0_11'][2], p['C_C_C_C_2'][2], i), cal(p['0_11'][3], p['C_C_C_C_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CB', 'dc', 'CT', 'CT'), cal(p['0_2'][0], p['C_C_C_C_3'][0], i), cal(p['0_2'][1], p['C_C_C_C_3'][1], i), cal(p['0_2'][2], p['C_C_C_C_3'][2], i), cal(p['0_2'][3], p['C_C_C_C_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CB', 'dc', 'CT', 'HC'), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['H4'][2], p['0_H'][2], i), cal(p['H4'][3], p['0_H'][3], i))

def stock_add_to_all_VLX(dist, newcar='dc', gonhyd='dh'):
        p = {}
        with open('Param_files/Stock/Stock.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
                p[line.split()[0]] = []
                for point in line.split()[1:]:
                        p[line.split()[0]].append(float(point))
        b.close()
        for i in range(11):
                a = i*10
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['CT'][0], p['CT'][0], i), cal(p['CT'][1], p['CT'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, gonhyd), cal(p['CT_HC'][0], p['H_sCT'][0], i), cal(p['CT_HC'][1], p['H_sCT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'CT'), cal(p['CT_CT'][0], p['CT_CT'][0], i), cal(p['CT_CT'][1], p['CT_CT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'HC'), cal(p['CT_HC'][0], p['CT_HC'][0], i), cal(p['CT_HC'][1], p['CT_HC'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', newcar), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', newcar, 'HC'), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', newcar, 'HC'), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, 'CT', 'HC'), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', newcar, gonhyd), cal(p['H_C_H'][0], p['X_CM_CM_0'][0], i), cal(p['H_C_H'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', newcar, gonhyd), cal(p['C_C_H'][0], p['X_CM_CM_0'][0], i), cal(p['C_C_H'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'CT', 'HC'), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'CT', 'CT'), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', newcar, 'CT', 'HC'), cal(p['H_C_C_H'][0], p['H_C_C_H'][0], i), cal(p['H_C_C_H'][1], p['H_C_C_H'][1], i), cal(p['H_C_C_H'][2], p['H_C_C_H'][2], i), cal(p['H_C_C_H'][3], p['H_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', newcar, 'CT', 'CT'), cal(p['C_C_C_H'][0], p['C_C_C_H'][0], i), cal(p['C_C_C_H'][1], p['C_C_C_H'][1], i), cal(p['C_C_C_H'][2], p['C_C_C_H'][2], i), cal(p['C_C_C_H'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['CT'][2], p['CT'][2], i), cal(p['CT'][3], p['CT'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))

def stock_add_to_all_CYX(dist, newsul='ds', gonhyd='dh'):
        p = {}
        with open('Param_files/Stock/Stock.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
                p[line.split()[0]] = []
                for point in line.split()[1:]:
                        p[line.split()[0]].append(float(point))
        b.close()
        for i in range(11):
                a = i*10
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newsul, cal(p['SH'][0], p['SH'][0], i), cal(p['SH'][1], p['SH'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HS'][0], p['0_H'][0], i), cal(p['HS'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newsul, 'ds'), cal(p['0_0'][0], p['S_S'][0], i), cal(p['0_0'][1], p['S_S'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newsul, gonhyd), cal(p['SH_HS'][0], p['H_sSS'][0], i), cal(p['SH_HS'][1], p['H_sSS'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newsul, 'CT'), cal(p['CT_SH'][0], p['CT_S'][0], i), cal(p['CT_SH'][1], p['CT_S'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', newsul), cal(p['C_C_SH'][0], p['C_C_S'][0], i), cal(p['C_C_SH'][1], p['C_C_S'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newsul, 'CT', 'H1'), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', newsul, gonhyd), cal(p['C_SH_H'][0], p['X_CM_CM_0'][0], i), cal(p['C_SH_H'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newsul, 'ds', 'CT'), cal(p['X_CM_CM_0'][0], p['C_S_S'][0], i), cal(p['X_CM_CM_0'][1], p['C_S_S'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(gonhyd, newsul, 'ds'), cal(p['X_CM_CM_0'][0], p['X_CM_CM_0'][0], i), cal(p['X_CM_CM_0'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newsul, 'CT', 'H1'), cal(p['X_C_SH_X'][0], p['0_5'][0], i), cal(p['X_C_SH_X'][1], p['0_5'][1], i), cal(p['X_C_SH_X'][2], p['0_5'][2], i), cal(p['X_C_SH_X'][3], p['0_5'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newsul, 'CT', 'CT'), cal(p['X_C_SH_X'][0], p['0_5'][0], i), cal(p['X_C_SH_X'][1], p['0_5'][1], i), cal(p['X_C_SH_X'][2], p['0_5'][2], i), cal(p['X_C_SH_X'][3], p['0_5'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', newsul, 'ds', 'CT'), cal(p['0_1'][0], p['C_S_S_C_2'][0], i), cal(p['0_1'][1], p['C_S_S_C_2'][1], i), cal(p['0_1'][2], p['C_S_S_C_2'][2], i), cal(p['0_1'][3], p['C_S_S_C_2'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', newsul, 'ds', 'CT'), cal(p['0_8'][0], p['C_S_S_C_1'][0], i), cal(p['0_8'][1], p['C_S_S_C_1'][1], i), cal(p['0_8'][2], p['C_S_S_C_1'][2], i), cal(p['0_8'][3], p['C_S_S_C_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', 'CT', newsul, 'ds'), cal(p['0_5'][0], p['X_C_SH_X'][0], i), cal(p['0_5'][1], p['X_C_SH_X'][1], i), cal(p['0_5'][2], p['X_C_SH_X'][2], i), cal(p['0_5'][3], p['X_C_SH_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', newsul, 'ds'), cal(p['0_5'][0], p['X_C_SH_X'][0], i), cal(p['0_5'][1], p['X_C_SH_X'][1], i), cal(p['0_5'][2], p['X_C_SH_X'][2], i), cal(p['0_5'][3], p['X_C_SH_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newsul, 'ds', 'CT'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newsul, 'ds', 'dh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newsul, cal(p['SH'][2], p['SH'][2], i), cal(p['SH'][3], p['SH'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HS'][2], p['0_H'][2], i), cal(p['HS'][3], p['0_H'][3], i))

def stock_add_to_all_MCX(dist, newcar='dc', gonhyd='gh'):
        p = {}
        with open('Param_files/Stock/Stock.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
                p[line.split()[0]] = []
                for point in line.split()[1:]:
                        p[line.split()[0]].append(float(point))
        b.close()
        for i in range(11):
                a = i*10
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['CT'][0], p['CT'][0], i), cal(p['CT'][1], p['CT'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'ds'), cal(p['0_0'][0], p['CT_S'][0], i), cal(p['0_0'][1], p['CT_S'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'dc'), cal(p['0_0'][0], p['CT_CT'][0], i), cal(dist, p['CT_CT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, gonhyd), cal(p['CT_HC'][0], p['HC_sS'][0], i), cal(p['CT_HC'][1], p['HC_sS'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'H1'), cal(p['CT_HC'][0], p['CT_HC'][0], i), cal(p['CT_HC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'S '), cal(p['CT_S'][0], p['CT_S'][0], i), cal(p['CT_S'][1], p['CT_S'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'S ', newcar), cal(p['C_S_C'][0], p['C_S_C'][0], i), cal(p['C_S_C'][1], p['C_S_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', newcar, 'H1'), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(gonhyd, newcar, 'H1'), cal(p['H_C_H'][0], p['C_C_H'][0], i), cal(p['H_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('S ', newcar, 'H1'), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('S ', newcar, gonhyd), cal(p['C_C_H'][0], p['X_CM_CM_0'][0], i), cal(p['C_C_H'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, 'ds', 'CT'), cal(p['X_CM_CM_0'][0], p['C_S_C'][0], i), cal(p['X_CM_CM_0'][1], p['C_S_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('S ', newcar, 'ds'), cal(p['X_CM_CM_0'][0], p['S_C_S'][0], i), cal(p['X_CM_CM_0'][1], p['S_C_S'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('S ', newcar, 'dc'), cal(p['X_CM_CM_0'][0], p['C_C_S'][0], i), cal(p['X_CM_CM_0'][1], p['C_C_S'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', newcar, 'ds'), cal(p['X_CM_CM_0'][0], p['C_C_H'][0], i), cal(p['X_CM_CM_0'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', newcar, 'dc'), cal(p['X_CM_CM_0'][0], p['C_C_H'][0], i), cal(p['X_CM_CM_0'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(gonhyd, newcar, 'ds'), cal(p['X_CM_CM_0'][0], p['X_CM_CM_0'][0], i), cal(p['X_CM_CM_0'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(gonhyd, newcar, 'dc'), cal(p['X_CM_CM_0'][0], p['X_CM_CM_0'][0], i), cal(p['X_CM_CM_0'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, 'ds', 'dh'), cal(p['X_CM_CM_0'][0], p['X_CM_CM_0'][0], i), cal(p['X_CM_CM_0'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'S ', 'CT'), cal(p['X_C_S_X'][0], p['0_5'][0], i), cal(p['X_C_S_X'][1], p['0_5'][1], i), cal(p['X_C_S_X'][2], p['0_5'][2], i), cal(p['X_C_S_X'][3], p['0_5'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', newcar, 'S ', 'CT'), cal(p['X_C_S_X'][0], p['X_C_S_X'][0], i), cal(p['X_C_S_X'][1], p['X_C_S_X'][1], i), cal(p['X_C_S_X'][2], p['X_C_S_X'][2], i), cal(p['X_C_S_X'][3], p['X_C_S_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', newcar, 'dc', 'S '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', newcar, 'dc', 'H1'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', newcar, 'dc', 'H1'), cal(p['0_4'][0], p['H_C_C_H'][0], i), cal(p['0_4'][1], p['H_C_C_H'][1], i), cal(p['0_4'][2], p['H_C_C_H'][2], i), cal(p['0_4'][3], p['H_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newcar, 'dc', 'S ', 'CT'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'S ', newcar, 'ds'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', newcar, 'ds', 'CT'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', newcar, 'ds', 'CT'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newcar, 'ds', 'CT', 'CT'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newcar, 'ds', 'CT', 'H1'), cal(p['0_5'][0], p['X_C_S_X'][0], i), cal(p['0_5'][1], p['X_C_S_X'][1], i), cal(p['0_5'][2], p['X_C_S_X'][2], i), cal(p['0_5'][3], p['X_C_S_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'ds', 'CT'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'dc', 'S '), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'ds', 'dh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'dc', 'gh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', newcar, 'ds', 'dh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('H1', newcar, 'dc', 'gh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', newcar, 'ds', 'dh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('S ', newcar, 'dc', 'gh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['CT'][2], p['CT'][2], i), cal(p['CT'][3], p['CT'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))

def stock_add_to_all_ABX(dist, newcar='dc', gonhyd='dh'):
        p = {}
        with open('Param_files/Stock/Stock.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
                p[line.split()[0]] = []
                for point in line.split()[1:]:
                        p[line.split()[0]].append(float(point))
        b.close()
        for i in range(11):
                a = i*10
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['CT'][0], p['CT'][0], i), cal(p['CT'][1], p['CT'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, gonhyd), cal(p['CT_HC'][0], p['H_sCT'][0], i), cal(p['CT_HC'][1], p['H_sCT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'dc'), cal(p['0_0'][0], p['CT_CT'][0], i), cal(dist, p['CT_CT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'CT'), cal(p['CT_CT'][0], p['CT_CT'][0], i), cal(p['CT_CT'][1], p['CT_CT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'HC'), cal(p['CT_HC'][0], p['CT_HC'][0], i), cal(p['CT_HC'][1], p['CT_HC'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', newcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
		Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', newcar), cal(p['C_C_C'][0], p['C_C_C'][0], i), cal(p['C_C_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', newcar, 'HC'), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', newcar, 'HC'), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('dc', newcar, 'HC'), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', newcar, gonhyd), cal(p['H_C_H'][0], p['X_CM_CM_0'][0], i), cal(p['H_C_H'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', newcar, gonhyd), cal(p['C_C_H'][0], p['X_CM_CM_0'][0], i), cal(p['C_C_H'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, 'dc', 'CT'), cal(p['X_CM_CM_0'][0], p['C_C_C'][0], i), cal(p['X_CM_CM_0'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, 'dc', 'HC'), cal(p['X_CM_CM_0'][0], p['C_C_H'][0], i), cal(p['X_CM_CM_0'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(gonhyd, newcar, 'dc'), cal(p['X_CM_CM_0'][0], p['X_CM_CM_0'][0], i), cal(p['X_CM_CM_0'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'dc', newcar, 'CT'), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'dc', newcar, 'HC'), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('dh', 'dc', newcar, gonhyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'dc', newcar, gonhyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'dc', newcar, gonhyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', newcar, 'dc'), cal(p['0_12'][0], p['C_C_C_C_1'][0], i), cal(p['0_12'][1], p['C_C_C_C_1'][1], i), cal(p['0_12'][2], p['C_C_C_C_1'][2], i), cal(p['0_12'][3], p['C_C_C_C_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', newcar, 'dc'), cal(p['0_11'][0], p['C_C_C_C_2'][0], i), cal(p['0_11'][1], p['C_C_C_C_2'][1], i), cal(p['0_11'][2], p['C_C_C_C_2'][2], i), cal(p['0_11'][3], p['C_C_C_C_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', newcar, 'dc'), cal(p['0_2'][0], p['C_C_C_C_3'][0], i), cal(p['0_2'][1], p['C_C_C_C_3'][1], i), cal(p['0_2'][2], p['C_C_C_C_3'][2], i), cal(p['0_2'][3], p['C_C_C_C_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', newcar, 'dc', 'CT'), cal(p['0_12'][0], p['C_C_C_C_1'][0], i), cal(p['0_12'][1], p['C_C_C_C_1'][1], i), cal(p['0_12'][2], p['C_C_C_C_1'][2], i), cal(p['0_12'][3], p['C_C_C_C_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', newcar, 'dc', 'CT'), cal(p['0_11'][0], p['C_C_C_C_2'][0], i), cal(p['0_11'][1], p['C_C_C_C_2'][1], i), cal(p['0_11'][2], p['C_C_C_C_2'][2], i), cal(p['0_11'][3], p['C_C_C_C_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', newcar, 'dc', 'CT'), cal(p['0_2'][0], p['C_C_C_C_3'][0], i), cal(p['0_2'][1], p['C_C_C_C_3'][1], i), cal(p['0_2'][2], p['C_C_C_C_3'][2], i), cal(p['0_2'][3], p['C_C_C_C_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('dc', newcar, 'CT', 'HC'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'CT', 'HC'), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'CT', 'CT'), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', newcar, 'CT', 'HC'), cal(p['X_C_C_X'][0], p['X_C_C_X'][0], i), cal(p['X_C_C_X'][1], p['X_C_C_X'][1], i), cal(p['X_C_C_X'][2], p['X_C_C_X'][2], i), cal(p['X_C_C_X'][3], p['X_C_C_X'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', newcar, 'CT', 'CT'), cal(p['C_C_C_H'][0], p['C_C_C_H'][0], i), cal(p['C_C_C_H'][1], p['C_C_C_H'][1], i), cal(p['C_C_C_H'][2], p['C_C_C_H'][2], i), cal(p['C_C_C_H'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['CT'][2], p['CT'][2], i), cal(p['CT'][3], p['CT'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))

def stock_add_to_all_ALX(dist, newcar='dc', gonhyd='dh'):
        p = {}
        with open('Param_files/Stock/Stock.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
                p[line.split()[0]] = []
                for point in line.split()[1:]:
                        p[line.split()[0]].append(float(point))
        b.close()
        for i in range(11):
                a = i*10
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['CT'][0], p['CT'][0], i), cal(p['CT'][1], p['CT'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, gonhyd), cal(p['CT_HC'][0], p['H_sCT'][0], i), cal(p['CT_HC'][1], p['H_sCT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'dc'), cal(p['0_0'][0], p['CT_CT'][0], i), cal(dist, p['CT_CT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'CT'), cal(p['CT_CT'][0], p['CT_CT'][0], i), cal(p['CT_CT'][1], p['CT_CT'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'HC'), cal(p['CT_HC'][0], p['CT_HC'][0], i), cal(p['CT_HC'][1], p['CT_HC'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('H1', 'CT', newcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C', 'CT', newcar), cal(p['CT_CT_C'][0], p['CT_CT_C'][0], i), cal(p['CT_CT_C'][1], p['CT_CT_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('N', 'CT', newcar), cal(p['CT_CT_N'][0], p['CT_CT_N'][0], i), cal(p['CT_CT_N'][1], p['CT_CT_N'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', newcar, 'HC'), cal(p['H_C_H'][0], p['H_C_H'][0], i), cal(p['H_C_H'][1], p['H_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', newcar, 'HC'), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('dc', newcar, 'HC'), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', newcar, gonhyd), cal(p['H_C_H'][0], p['X_CM_CM_0'][0], i), cal(p['H_C_H'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', newcar, gonhyd), cal(p['C_C_H'][0], p['X_CM_CM_0'][0], i), cal(p['C_C_H'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, 'dc', 'CT'), cal(p['X_CM_CM_0'][0], p['C_C_C'][0], i), cal(p['X_CM_CM_0'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, 'dc', 'HC'), cal(p['X_CM_CM_0'][0], p['C_C_H'][0], i), cal(p['X_CM_CM_0'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(gonhyd, newcar, 'dc'), cal(p['X_CM_CM_0'][0], p['X_CM_CM_0'][0], i), cal(p['X_CM_CM_0'][1], p['X_CM_CM_0'][1], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'dc', newcar, 'CT'), cal(p['0_1'][0], p['C_C_C_H'][0], i), cal(p['0_1'][1], p['C_C_C_H'][1], i), cal(p['0_1'][2], p['C_C_C_H'][2], i), cal(p['0_1'][3], p['C_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'dc', newcar, 'HC'), cal(p['0_1'][0], p['H_C_C_H'][0], i), cal(p['0_1'][1], p['H_C_C_H'][1], i), cal(p['0_1'][2], p['H_C_C_H'][2], i), cal(p['0_1'][3], p['H_C_C_H'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('dh', 'dc', newcar, gonhyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'dc', newcar, gonhyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'dc', newcar, gonhyd), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', newcar, 'dc', 'CT'), cal(p['0_12'][0], p['C_C_C_C_1'][0], i), cal(p['0_12'][1], p['C_C_C_C_1'][1], i), cal(p['0_12'][2], p['C_C_C_C_1'][2], i), cal(p['0_12'][3], p['C_C_C_C_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', newcar, 'dc', 'CT'), cal(p['0_11'][0], p['C_C_C_C_2'][0], i), cal(p['0_11'][1], p['C_C_C_C_2'][1], i), cal(p['0_11'][2], p['C_C_C_C_2'][2], i), cal(p['0_11'][3], p['C_C_C_C_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', newcar, 'dc', 'CT'), cal(p['0_2'][0], p['C_C_C_C_3'][0], i), cal(p['0_2'][1], p['C_C_C_C_3'][1], i), cal(p['0_2'][2], p['C_C_C_C_3'][2], i), cal(p['0_2'][3], p['C_C_C_C_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('dc', newcar, 'CT', 'N '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('dc', newcar, 'CT', 'C '), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('dc', newcar, 'CT', 'H1'), cal(p['0_4'][0], p['X_C_C_X'][0], i), cal(p['0_4'][1], p['X_C_C_X'][1], i), cal(p['0_4'][2], p['X_C_C_X'][2], i), cal(p['0_4'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'CT', 'N '), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'CT', 'C '), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'CT', 'H1'), cal(p['X_C_C_X'][0], p['0_4'][0], i), cal(p['X_C_C_X'][1], p['0_4'][1], i), cal(p['X_C_C_X'][2], p['0_4'][2], i), cal(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', newcar, 'CT', 'N '), cal(p['X_C_C_X'][0], p['X_C_C_X'][0], i), cal(p['X_C_C_X'][1], p['X_C_C_X'][1], i), cal(p['X_C_C_X'][2], p['X_C_C_X'][2], i), cal(p['X_C_C_X'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', newcar, 'CT', 'C '), cal(p['X_C_C_X'][0], p['X_C_C_X'][0], i), cal(p['X_C_C_X'][1], p['X_C_C_X'][1], i), cal(p['X_C_C_X'][2], p['X_C_C_X'][2], i), cal(p['X_C_C_X'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', newcar, 'CT', 'H1'), cal(p['X_C_C_X'][0], p['X_C_C_X'][0], i), cal(p['X_C_C_X'][1], p['X_C_C_X'][1], i), cal(p['X_C_C_X'][2], p['X_C_C_X'][2], i), cal(p['X_C_C_X'][3], p['X_C_C_X'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['CT'][2], p['CT'][2], i), cal(p['CT'][3], p['CT'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))

def stock_add_to_all_I4X(dist, metcar='mc', methyd='mh', hydhyd='hh', dobcar='dc', dobhyd='dh', gonhyd='gh'):
        p = {}
        with open('Param_files/Stock/Stock.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
                p[line.split()[0]] = []
                for point in line.split()[1:]:
                        p[line.split()[0]].append(float(point))
        b.close()
        for i in range(11):
                a = i*10
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), metcar, lac(p['CT'][0], p['0_C'][0], i), lac(p['CT'][1], p['0_C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), methyd, lac(p['HC'][0], p['0_H'][0], i), lac(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, lac(p['0_H'][0], p['H1'][0], i), lac(p['0_H'][1], p['H1'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', metcar), lac(p['CT_CT'][0], p['CT_mH'][0], i), lac(p['CT_CT'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd), lac(p['HC_sC'][0], p['CT_HC'][0], i), lac(p['HC_sC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(metcar, methyd), lac(p['CT_HC'][0], p['HC_mH'][0], i), lac(p['CT_HC'][1], p['HC_mH'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', metcar, methyd), lac(p['C_C_H'][0], p['Dritt'][0], i), lac(p['C_C_H'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(methyd, metcar, methyd), lac(p['H_C_H'][0], p['Close'][0], i), lac(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C', 'CT', metcar), lac(p['CT_CT_C'][0], p['C_C_C'][0], i), lac(p['CT_CT_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C', 'CT', hydhyd), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('N', 'CT', metcar), lac(p['CT_CT_N'][0], p['C_C_C'][0], i), lac(p['CT_CT_N'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('N', 'CT', hydhyd), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', metcar), lac(p['C_C_C'][0], p['C_C_C'][0], i), lac(p['C_C_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd, 'CT', metcar), lac(p['Close'][0], p['Close'][0], i), lac(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), lac(p['C_C_N_C_1'][0], p['0_3'][0], i), lac(p['C_C_N_C_1'][1], p['0_3'][1], i), lac(p['C_C_N_C_1'][2], p['0_3'][2], i), lac(p['C_C_N_C_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), lac(p['C_C_N_C_2'][0], p['0_8'][0], i), lac(p['C_C_N_C_2'][1], p['0_8'][1], i), lac(p['C_C_N_C_2'][2], p['0_8'][2], i), lac(p['C_C_N_C_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), lac(p['C_C_N_C_3'][0], p['0_2'][0], i), lac(p['C_C_N_C_3'][1], p['0_2'][1], i), lac(p['C_C_N_C_3'][2], p['0_2'][2], i), lac(p['C_C_N_C_3'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), lac(p['C_C_N_C_4'][0], p['0_7'][0], i), lac(p['C_C_N_C_4'][1], p['0_7'][1], i), lac(p['C_C_N_C_4'][2], p['0_7'][2], i), lac(p['C_C_N_C_4'][3], p['0_7'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), lac(p['C_C_C_N_1'][0], p['0_3'][0], i), lac(p['C_C_C_N_1'][1], p['0_3'][1], i), lac(p['C_C_C_N_1'][2], p['0_3'][2], i), lac(p['C_C_C_N_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), lac(p['C_C_C_N_2'][0], p['0_8'][0], i), lac(p['C_C_C_N_2'][1], p['0_8'][1], i), lac(p['C_C_C_N_2'][2], p['0_8'][2], i), lac(p['C_C_C_N_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), lac(p['C_C_C_N_3'][0], p['0_2'][0], i), lac(p['C_C_C_N_3'][1], p['0_2'][1], i), lac(p['C_C_C_N_3'][2], p['0_2'][2], i), lac(p['C_C_C_N_3'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), lac(p['C_C_C_N_4'][0], p['0_7'][0], i), lac(p['C_C_C_N_4'][1], p['0_7'][1], i), lac(p['C_C_C_N_4'][2], p['0_7'][2], i), lac(p['C_C_C_N_4'][3], p['0_7'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', metcar, methyd), lac(p['C_C_C_H'][0], p['0_1'][0], i), lac(p['C_C_C_H'][1], p['0_1'][1], i), lac(p['C_C_C_H'][2], p['0_1'][2], i), lac(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd, 'CT', metcar, methyd), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'CT', metcar, methyd), lac(p['X_C_C_X'][0], p['0_4'][0], i), lac(p['X_C_C_X'][1], p['0_4'][1], i), lac(p['X_C_C_X'][2], p['0_4'][2], i), lac(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'CT', metcar, methyd), lac(p['X_C_C_X'][0], p['0_4'][0], i), lac(p['X_C_C_X'][1], p['0_4'][1], i), lac(p['X_C_C_X'][2], p['0_4'][2], i), lac(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), metcar, lac(p['CT'][2], p['0_C'][2], i), lac(p['CT'][3], p['0_C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), methyd, lac(p['HC'][2], p['0_H'][2], i), lac(p['HC'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, lac(p['0_H'][2], p['H1'][2], i), lac(p['0_H'][3], p['H1'][3], i))

                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), dobcar, cal(p['CT'][0], p['CA'][0], i), cal(p['CT'][1], p['CA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), dobhyd, cal(p['HC'][0], p['HC'][0], i), cal(p['HC'][1], p['HC'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', dobcar), cal(p['CT_CT'][0], p['CT_CA'][0], i), cal(p['CT_CT'][1], p['CT_CA'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(dobcar, 'dc'), cal(p['0_0'][0], p['CM_CM'][0], i), cal(dist, p['CM_CM'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(dobcar, dobhyd), cal(p['CT_HC'][0], p['CA_HA'][0], i), cal(p['CT_HC'][1], p['CA_HA'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(dobcar, gonhyd), cal(p['CT_HC'][0], p['CA_HA2'][0], i), cal(p['CT_HC'][1], p['CA_HA2'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', dobcar), cal(p['C_C_C'][0], p['CT_CT_C'][0], i), cal(p['C_C_C'][1], p['CT_CT_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', dobcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', dobcar, dobhyd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', dobcar, gonhyd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(dobhyd, dobcar, gonhyd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['C_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(gonhyd, dobcar, gonhyd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('dc', dobcar, gonhyd), cal(p['X_CM_CM_0'][0], p['H_CM_CM'][0], i), cal(p['X_CM_CM_0'][1], p['H_CM_CM'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('dc', dobcar, dobhyd), cal(p['X_CM_CM_0'][0], p['H_CM_CM'][0], i), cal(p['X_CM_CM_0'][1], p['H_CM_CM'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(dobcar, 'dc', 'CT'), cal(p['X_CM_CM_0'][0], p['C_CM_CM'][0], i), cal(p['X_CM_CM_0'][1], p['C_CM_CM'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', dobcar, dobhyd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', 'CT', dobcar), lac(p['C_C_C_CM_1'][0], p['0_1'][0], i), lac(p['C_C_C_CM_1'][1], p['0_1'][1], i), lac(p['C_C_C_CM_1'][2], p['0_1'][2], i), lac(p['C_C_C_CM_1'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', 'CT', dobcar), lac(p['C_C_C_CM_2'][0], p['0_11'][0], i), lac(p['C_C_C_CM_2'][1], p['0_11'][1], i), lac(p['C_C_C_CM_2'][2], p['0_11'][2], i), lac(p['C_C_C_CM_2'][3], p['0_11'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', 'CT', dobcar), lac(p['C_C_C_CM_3'][0], p['0_13'][0], i), lac(p['C_C_C_CM_3'][1], p['0_13'][1], i), lac(p['C_C_C_CM_3'][2], p['0_13'][2], i), lac(p['C_C_C_CM_3'][3], p['0_13'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, 'dc'), lac(p['C_C_CM_CM_1'][0], p['0_1'][0], i), lac(p['C_C_CM_CM_1'][1], p['0_1'][1], i), lac(p['C_C_CM_CM_1'][2], p['0_1'][2], i), lac(p['C_C_CM_CM_1'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, 'dc'), lac(p['C_C_CM_CM_2'][0], p['0_8'][0], i), lac(p['C_C_CM_CM_2'][1], p['0_8'][1], i), lac(p['C_C_CM_CM_2'][2], p['0_8'][2], i), lac(p['C_C_CM_CM_2'][3], p['0_8'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, 'dc'), lac(p['C_C_CM_CM_3'][0], p['0_9'][0], i), lac(p['C_C_CM_CM_3'][1], p['0_9'][1], i), lac(p['C_C_CM_CM_3'][2], p['0_9'][2], i), lac(p['C_C_CM_CM_3'][3], p['0_9'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', dobcar, 'dc'), lac(p['H_C_CM_CM_1'][0], p['0_3'][0], i), lac(p['C_C_CM_CM_1'][1], p['0_3'][1], i), lac(p['C_C_CM_CM_1'][2], p['0_3'][2], i), lac(p['C_C_CM_CM_1'][3], p['0_3'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', dobcar, 'dc'), lac(p['H_C_CM_CM_2'][0], p['0_14'][0], i), lac(p['C_C_CM_CM_1'][1], p['0_14'][1], i), lac(p['C_C_CM_CM_1'][2], p['0_14'][2], i), lac(p['C_C_CM_CM_1'][3], p['0_14'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, dobhyd), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, gonhyd), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', dobcar, dobhyd), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', dobcar, gonhyd), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, 'dc', 'CT'), lac(p['C_CM_CM_C_1'][0], p['0_1'][0], i), lac(p['C_CM_CM_C_1'][1], p['0_1'][1], i), lac(p['C_CM_CM_C_1'][2], p['0_1'][2], i), lac(p['C_CM_CM_C_1'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, 'dc', 'CT'), lac(p['C_CM_CM_C_2'][0], p['0_11'][0], i), lac(p['C_CM_CM_C_2'][1], p['0_11'][1], i), lac(p['C_CM_CM_C_2'][2], p['0_11'][2], i), lac(p['C_CM_CM_C_2'][3], p['0_11'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, 'dc', 'CT'), lac(p['C_CM_CM_C_3'][0], p['0_13'][0], i), lac(p['C_CM_CM_C_3'][1], p['0_13'][1], i), lac(p['C_CM_CM_C_3'][2], p['0_13'][2], i), lac(p['C_CM_CM_C_3'][3], p['0_13'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'dc', dobcar, dobhyd), cal(p['0_15'][0], p['0_15'][0], i), cal(p['0_15'][1], p['0_15'][1], i), cal(p['0_15'][2], p['0_15'][2], i), cal(p['0_15'][3], p['0_15'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(dobhyd, dobcar, 'dc', 'dh'), cal(p['0_16'][0], p['0_16'][0], i), cal(p['0_16'][1], p['0_16'][1], i), cal(p['0_16'][2], p['0_16'][2], i), cal(p['0_16'][3], p['0_16'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, 'dc', 'gh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, dobcar, 'dc', 'gh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(dobhyd, dobcar, 'dc', 'gh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), dobcar, cal(p['CT'][2], p['CA'][2], i), cal(p['CT'][3], p['CA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), dobhyd, cal(p['HC'][2], p['HC'][2], i), cal(p['HC'][3], p['HC'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))

def stock_add_to_all_ADX(dist, newcar='dc', newhyd='dh', gonhyd='gh'):
        p = {}
        with open('Param_files/Stock/Stock.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
                p[line.split()[0]] = []
                for point in line.split()[1:]:
                        p[line.split()[0]].append(float(point))
        b.close()
        for i in range(11):
                a = i*10
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['CT'][0], p['CA'][0], i), cal(p['CT'][1], p['CA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, cal(p['HC'][0], p['HA'][0], i), cal(p['HC'][1], p['HA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', newcar), cal(p['CT_CT'][0], p['CT_CA'][0], i), cal(p['CT_CT'][1], p['CT_CA'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, 'dc'), cal(p['0_0'][0], p['CM_CM'][0], i), cal(dist, p['CM_CM'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, newhyd), cal(p['CT_HC'][0], p['CA_HA'][0], i), cal(p['CT_HC'][1], p['CA_HA'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(newcar, gonhyd), cal(p['CT_HC'][0], p['CA_HA2'][0], i), cal(p['CT_HC'][1], p['CA_HA2'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', newcar), cal(p['C_C_C'][0], p['CT_CT_C'][0], i), cal(p['C_C_C'][1], p['CT_CT_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', newcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', newcar, newhyd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', newcar, gonhyd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newhyd, newcar, gonhyd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['C_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(gonhyd, newcar, gonhyd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('dc', newcar, gonhyd), cal(p['X_CM_CM_0'][0], p['H_CM_CM'][0], i), cal(p['X_CM_CM_0'][1], p['H_CM_CM'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('dc', newcar, newhyd), cal(p['X_CM_CM_0'][0], p['H_CM_CM'][0], i), cal(p['X_CM_CM_0'][1], p['H_CM_CM'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(newcar, 'dc', 'CT'), cal(p['X_CM_CM_0'][0], p['C_CM_CM'][0], i), cal(p['X_CM_CM_0'][1], p['C_CM_CM'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', newcar, newhyd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', 'CT', newcar), lac(p['C_C_C_CM_1'][0], p['0_1'][0], i), lac(p['C_C_C_CM_1'][1], p['0_1'][1], i), lac(p['C_C_C_CM_1'][2], p['0_1'][2], i), lac(p['C_C_C_CM_1'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', 'CT', newcar), lac(p['C_C_C_CM_2'][0], p['0_11'][0], i), lac(p['C_C_C_CM_2'][1], p['0_11'][1], i), lac(p['C_C_C_CM_2'][2], p['0_11'][2], i), lac(p['C_C_C_CM_2'][3], p['0_11'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', 'CT', newcar), lac(p['C_C_C_CM_3'][0], p['0_13'][0], i), lac(p['C_C_C_CM_3'][1], p['0_13'][1], i), lac(p['C_C_C_CM_3'][2], p['0_13'][2], i), lac(p['C_C_C_CM_3'][3], p['0_13'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', newcar, 'dc'), lac(p['C_C_CM_CM_1'][0], p['0_1'][0], i), lac(p['C_C_CM_CM_1'][1], p['0_1'][1], i), lac(p['C_C_CM_CM_1'][2], p['0_1'][2], i), lac(p['C_C_CM_CM_1'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', newcar, 'dc'), lac(p['C_C_CM_CM_2'][0], p['0_8'][0], i), lac(p['C_C_CM_CM_2'][1], p['0_8'][1], i), lac(p['C_C_CM_CM_2'][2], p['0_8'][2], i), lac(p['C_C_CM_CM_2'][3], p['0_8'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', newcar, 'dc'), lac(p['C_C_CM_CM_3'][0], p['0_9'][0], i), lac(p['C_C_CM_CM_3'][1], p['0_9'][1], i), lac(p['C_C_CM_CM_3'][2], p['0_9'][2], i), lac(p['C_C_CM_CM_3'][3], p['0_9'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', newcar, 'dc'), lac(p['H_C_CM_CM_1'][0], p['0_3'][0], i), lac(p['C_C_CM_CM_1'][1], p['0_3'][1], i), lac(p['C_C_CM_CM_1'][2], p['0_3'][2], i), lac(p['C_C_CM_CM_1'][3], p['0_3'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', newcar, 'dc'), lac(p['H_C_CM_CM_2'][0], p['0_14'][0], i), lac(p['C_C_CM_CM_1'][1], p['0_14'][1], i), lac(p['C_C_CM_CM_1'][2], p['0_14'][2], i), lac(p['C_C_CM_CM_1'][3], p['0_14'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', newcar, newhyd), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', newcar, gonhyd), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', newcar, newhyd), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', newcar, gonhyd), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', newcar, 'dc', 'CT'), lac(p['C_CM_CM_C_1'][0], p['0_1'][0], i), lac(p['C_CM_CM_C_1'][1], p['0_1'][1], i), lac(p['C_CM_CM_C_1'][2], p['0_1'][2], i), lac(p['C_CM_CM_C_1'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', newcar, 'dc', 'CT'), lac(p['C_CM_CM_C_2'][0], p['0_11'][0], i), lac(p['C_CM_CM_C_2'][1], p['0_11'][1], i), lac(p['C_CM_CM_C_2'][2], p['0_11'][2], i), lac(p['C_CM_CM_C_2'][3], p['0_11'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', newcar, 'dc', 'CT'), lac(p['C_CM_CM_C_3'][0], p['0_13'][0], i), lac(p['C_CM_CM_C_3'][1], p['0_13'][1], i), lac(p['C_CM_CM_C_3'][2], p['0_13'][2], i), lac(p['C_CM_CM_C_3'][3], p['0_13'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'dc', newcar, newhyd), cal(p['0_15'][0], p['0_15'][0], i), cal(p['0_15'][1], p['0_15'][1], i), cal(p['0_15'][2], p['0_15'][2], i), cal(p['0_15'][3], p['0_15'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, newcar, 'dc', 'dh'), cal(p['0_16'][0], p['0_16'][0], i), cal(p['0_16'][1], p['0_16'][1], i), cal(p['0_16'][2], p['0_16'][2], i), cal(p['0_16'][3], p['0_16'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', newcar, 'dc', 'gh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, newcar, 'dc', 'gh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(newhyd, newcar, 'dc', 'gh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newcar, cal(p['CT'][2], p['CA'][2], i), cal(p['CT'][3], p['CA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), newhyd, cal(p['HC'][2], p['HA'][2], i), cal(p['HC'][3], p['HA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))

def stock_add_to_all_I7X(dist, metcar='mc', methyd='mh', hydhyd='hh', dobcar='dc', dobhyd='dh', gonhyd='gh'):
        p = {}
        with open('Param_files/Stock/Stock.param', 'r') as b:
                data = b.readlines()[1:]
        for line in data:
                p[line.split()[0]] = []
                for point in line.split()[1:]:
                        p[line.split()[0]].append(float(point))
        b.close()
        for i in range(11):
                a = i*10
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), metcar, lac(p['CT'][0], p['0_C'][0], i), lac(p['CT'][1], p['0_C'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), methyd, lac(p['HC'][0], p['0_H'][0], i), lac(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, lac(p['0_H'][0], p['H1'][0], i), lac(p['0_H'][1], p['H1'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', metcar), lac(p['CT_CT'][0], p['CT_mH'][0], i), lac(p['CT_CT'][1], p['CT_mH'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', hydhyd), lac(p['HC_sC'][0], p['CT_HC'][0], i), lac(p['HC_sC'][1], p['CT_HC'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(metcar, methyd), lac(p['CT_HC'][0], p['HC_mH'][0], i), lac(p['CT_HC'][1], p['HC_mH'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', metcar, methyd), lac(p['C_C_H'][0], p['Dritt'][0], i), lac(p['C_C_H'][1], p['Dritt'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(methyd, metcar, methyd), lac(p['H_C_H'][0], p['Close'][0], i), lac(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C', 'CT', metcar), lac(p['CT_CT_C'][0], p['C_C_C'][0], i), lac(p['CT_CT_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('C', 'CT', hydhyd), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('N', 'CT', metcar), lac(p['CT_CT_N'][0], p['C_C_C'][0], i), lac(p['CT_CT_N'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('N', 'CT', hydhyd), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', metcar), lac(p['C_C_C'][0], p['C_C_C'][0], i), lac(p['C_C_C'][1], p['C_C_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', hydhyd), lac(p['C_C_H'][0], p['C_C_H'][0], i), lac(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(hydhyd, 'CT', metcar), lac(p['Close'][0], p['Close'][0], i), lac(p['Close'][1], p['Close'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), lac(p['C_C_N_C_1'][0], p['0_3'][0], i), lac(p['C_C_N_C_1'][1], p['0_3'][1], i), lac(p['C_C_N_C_1'][2], p['0_3'][2], i), lac(p['C_C_N_C_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), lac(p['C_C_N_C_2'][0], p['0_8'][0], i), lac(p['C_C_N_C_2'][1], p['0_8'][1], i), lac(p['C_C_N_C_2'][2], p['0_8'][2], i), lac(p['C_C_N_C_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), lac(p['C_C_N_C_3'][0], p['0_2'][0], i), lac(p['C_C_N_C_3'][1], p['0_2'][1], i), lac(p['C_C_N_C_3'][2], p['0_2'][2], i), lac(p['C_C_N_C_3'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'N ', 'CT', metcar), lac(p['C_C_N_C_4'][0], p['0_7'][0], i), lac(p['C_C_N_C_4'][1], p['0_7'][1], i), lac(p['C_C_N_C_4'][2], p['0_7'][2], i), lac(p['C_C_N_C_4'][3], p['0_7'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), lac(p['C_C_C_N_1'][0], p['0_3'][0], i), lac(p['C_C_C_N_1'][1], p['0_3'][1], i), lac(p['C_C_C_N_1'][2], p['0_3'][2], i), lac(p['C_C_C_N_1'][3], p['0_3'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), lac(p['C_C_C_N_2'][0], p['0_8'][0], i), lac(p['C_C_C_N_2'][1], p['0_8'][1], i), lac(p['C_C_C_N_2'][2], p['0_8'][2], i), lac(p['C_C_C_N_2'][3], p['0_8'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), lac(p['C_C_C_N_3'][0], p['0_2'][0], i), lac(p['C_C_C_N_3'][1], p['0_2'][1], i), lac(p['C_C_C_N_3'][2], p['0_2'][2], i), lac(p['C_C_C_N_3'][3], p['0_2'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'C ', 'CT', metcar), lac(p['C_C_C_N_4'][0], p['0_7'][0], i), lac(p['C_C_C_N_4'][1], p['0_7'][1], i), lac(p['C_C_C_N_4'][2], p['0_7'][2], i), lac(p['C_C_C_N_4'][3], p['0_7'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', metcar, methyd), lac(p['C_C_C_H'][0], p['0_1'][0], i), lac(p['C_C_C_H'][1], p['0_1'][1], i), lac(p['C_C_C_H'][2], p['0_1'][2], i), lac(p['C_C_C_H'][3], p['0_1'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(hydhyd, 'CT', metcar, methyd), lac(p['0_Dihe'][0], p['0_Dihe'][0], i), lac(p['0_Dihe'][1], p['0_Dihe'][1], i), lac(p['0_Dihe'][2], p['0_Dihe'][2], i), lac(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('N ', 'CT', metcar, methyd), lac(p['X_C_C_X'][0], p['0_4'][0], i), lac(p['X_C_C_X'][1], p['0_4'][1], i), lac(p['X_C_C_X'][2], p['0_4'][2], i), lac(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('C ', 'CT', metcar, methyd), lac(p['X_C_C_X'][0], p['0_4'][0], i), lac(p['X_C_C_X'][1], p['0_4'][1], i), lac(p['X_C_C_X'][2], p['0_4'][2], i), lac(p['X_C_C_X'][3], p['0_4'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), metcar, lac(p['CT'][2], p['0_C'][2], i), lac(p['CT'][3], p['0_C'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), methyd, lac(p['HC'][2], p['0_H'][2], i), lac(p['HC'][3], p['0_H'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), hydhyd, lac(p['0_H'][2], p['H1'][2], i), lac(p['0_H'][3], p['H1'][3], i))

                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), dobcar, cal(p['CT'][0], p['CA'][0], i), cal(p['CT'][1], p['CA'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), dobhyd, cal(p['HC'][0], p['HC'][0], i), cal(p['HC'][1], p['HC'][1], i))
                Frcmod_creator.MASS_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][0], p['0_H'][0], i), cal(p['HC'][1], p['0_H'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format('CT', dobcar), cal(p['CT_CT'][0], p['CT_CA'][0], i), cal(p['CT_CT'][1], p['CT_CA'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(dobcar, 'dc'), cal(p['0_0'][0], p['CM_CM'][0], i), cal(dist, p['CM_CM'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(dobcar, dobhyd), cal(p['CT_HC'][0], p['CA_HA'][0], i), cal(p['CT_HC'][1], p['CA_HA'][1], i))
                Frcmod_creator.BOND_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}'.format(dobcar, gonhyd), cal(p['CT_HC'][0], p['CA_HA2'][0], i), cal(p['CT_HC'][1], p['CA_HA2'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', 'CT', dobcar), cal(p['C_C_C'][0], p['CT_CT_C'][0], i), cal(p['C_C_C'][1], p['CT_CT_C'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('HC', 'CT', dobcar), cal(p['C_C_H'][0], p['C_C_H'][0], i), cal(p['C_C_H'][1], p['C_C_H'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', dobcar, dobhyd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', dobcar, gonhyd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(dobhyd, dobcar, gonhyd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['C_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(gonhyd, dobcar, gonhyd), cal(p['H_C_H'][0], p['Close'][0], i), cal(p['H_C_H'][1], p['Close'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('dc', dobcar, gonhyd), cal(p['X_CM_CM_0'][0], p['H_CM_CM'][0], i), cal(p['X_CM_CM_0'][1], p['H_CM_CM'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('dc', dobcar, dobhyd), cal(p['X_CM_CM_0'][0], p['H_CM_CM'][0], i), cal(p['X_CM_CM_0'][1], p['H_CM_CM'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format(dobcar, 'dc', 'CT'), cal(p['X_CM_CM_0'][0], p['C_CM_CM'][0], i), cal(p['X_CM_CM_0'][1], p['C_CM_CM'][1], i))
                Frcmod_creator.ANGLE_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}'.format('CT', dobcar, dobhyd), cal(p['C_C_H'][0], p['F_CA_CA_HA'][0], i), cal(p['C_C_H'][1], p['F_CA_CA_HA'][1], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', 'CT', dobcar), lac(p['C_C_C_CM_1'][0], p['0_1'][0], i), lac(p['C_C_C_CM_1'][1], p['0_1'][1], i), lac(p['C_C_C_CM_1'][2], p['0_1'][2], i), lac(p['C_C_C_CM_1'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', 'CT', dobcar), lac(p['C_C_C_CM_2'][0], p['0_11'][0], i), lac(p['C_C_C_CM_2'][1], p['0_11'][1], i), lac(p['C_C_C_CM_2'][2], p['0_11'][2], i), lac(p['C_C_C_CM_2'][3], p['0_11'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', 'CT', dobcar), lac(p['C_C_C_CM_3'][0], p['0_13'][0], i), lac(p['C_C_C_CM_3'][1], p['0_13'][1], i), lac(p['C_C_C_CM_3'][2], p['0_13'][2], i), lac(p['C_C_C_CM_3'][3], p['0_13'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, 'dc'), lac(p['C_C_CM_CM_1'][0], p['0_1'][0], i), lac(p['C_C_CM_CM_1'][1], p['0_1'][1], i), lac(p['C_C_CM_CM_1'][2], p['0_1'][2], i), lac(p['C_C_CM_CM_1'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, 'dc'), lac(p['C_C_CM_CM_2'][0], p['0_8'][0], i), lac(p['C_C_CM_CM_2'][1], p['0_8'][1], i), lac(p['C_C_CM_CM_2'][2], p['0_8'][2], i), lac(p['C_C_CM_CM_2'][3], p['0_8'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, 'dc'), lac(p['C_C_CM_CM_3'][0], p['0_9'][0], i), lac(p['C_C_CM_CM_3'][1], p['0_9'][1], i), lac(p['C_C_CM_CM_3'][2], p['0_9'][2], i), lac(p['C_C_CM_CM_3'][3], p['0_9'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', dobcar, 'dc'), lac(p['H_C_CM_CM_1'][0], p['0_3'][0], i), lac(p['C_C_CM_CM_1'][1], p['0_3'][1], i), lac(p['C_C_CM_CM_1'][2], p['0_3'][2], i), lac(p['C_C_CM_CM_1'][3], p['0_3'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', dobcar, 'dc'), lac(p['H_C_CM_CM_2'][0], p['0_14'][0], i), lac(p['C_C_CM_CM_1'][1], p['0_14'][1], i), lac(p['C_C_CM_CM_1'][2], p['0_14'][2], i), lac(p['C_C_CM_CM_1'][3], p['0_14'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, dobhyd), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'CT', dobcar, gonhyd), cal(p['C_C_C_H'][0], p['0_1'][0], i), cal(p['C_C_C_H'][1], p['0_1'][1], i), cal(p['C_C_C_H'][2], p['0_1'][2], i), cal(p['C_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', dobcar, dobhyd), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('HC', 'CT', dobcar, gonhyd), cal(p['H_C_C_H'][0], p['0_1'][0], i), cal(p['H_C_C_H'][1], p['0_1'][1], i), cal(p['H_C_C_H'][2], p['0_1'][2], i), cal(p['H_C_C_H'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, 'dc', 'CT'), lac(p['C_CM_CM_C_1'][0], p['0_1'][0], i), lac(p['C_CM_CM_C_1'][1], p['0_1'][1], i), lac(p['C_CM_CM_C_1'][2], p['0_1'][2], i), lac(p['C_CM_CM_C_1'][3], p['0_1'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, 'dc', 'CT'), lac(p['C_CM_CM_C_2'][0], p['0_11'][0], i), lac(p['C_CM_CM_C_2'][1], p['0_11'][1], i), lac(p['C_CM_CM_C_2'][2], p['0_11'][2], i), lac(p['C_CM_CM_C_2'][3], p['0_11'][3], i))
		Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, 'dc', 'CT'), lac(p['C_CM_CM_C_3'][0], p['0_13'][0], i), lac(p['C_CM_CM_C_3'][1], p['0_13'][1], i), lac(p['C_CM_CM_C_3'][2], p['0_13'][2], i), lac(p['C_CM_CM_C_3'][3], p['0_13'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', 'dc', dobcar, dobhyd), cal(p['0_15'][0], p['0_15'][0], i), cal(p['0_15'][1], p['0_15'][1], i), cal(p['0_15'][2], p['0_15'][2], i), cal(p['0_15'][3], p['0_15'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(dobhyd, dobcar, 'dc', 'dh'), cal(p['0_16'][0], p['0_16'][0], i), cal(p['0_16'][1], p['0_16'][1], i), cal(p['0_16'][2], p['0_16'][2], i), cal(p['0_16'][3], p['0_16'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format('CT', dobcar, 'dc', 'gh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(gonhyd, dobcar, 'dc', 'gh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.DIHEDRAL_insert('{}_{}.frcmod'.format(a, 100-a), '{}-{}-{}-{}'.format(dobhyd, dobcar, 'dc', 'gh'), cal(p['0_Dihe'][0], p['0_Dihe'][0], i), cal(p['0_Dihe'][1], p['0_Dihe'][1], i), cal(p['0_Dihe'][2], p['0_Dihe'][2], i), cal(p['0_Dihe'][3], p['0_Dihe'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), dobcar, cal(p['CT'][2], p['CA'][2], i), cal(p['CT'][3], p['CA'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), dobhyd, cal(p['HC'][2], p['HC'][2], i), cal(p['HC'][3], p['HC'][3], i))
                Frcmod_creator.NONBON_insert('{}_{}.frcmod'.format(a, 100-a), gonhyd, cal(p['HC'][2], p['0_H'][2], i), cal(p['HC'][3], p['0_H'][3], i))
