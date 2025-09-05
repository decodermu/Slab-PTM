import numpy as np
import pandas as pd
import mdtraj as md
import itertools
import os
import math
import MDAnalysis
from MDAnalysis import transformations
from argparse import ArgumentParser
from simtk import openmm, unit
from simtk.openmm import app
from simtk.openmm import XmlSerializer
import time

# -------------------
# PTM Transfer
# -------------------
def trans_ptm(seq, phos_list=None, ace_list=None, met_list=None):
    fasta = list(seq)  # fasta to list
    # check phos
    if phos_list:
        for position in phos_list:
            if fasta[position] == 'S':
                fasta[position] = 'Sp'
            elif fasta[position] == 'T':
                fasta[position] = 'Tp'
            elif fasta[position] == 'Y':
                fasta[position] = 'Yp'
            else:
                raise ValueError(f"Position {position} is not Ser(S) or Thr(T) or Tyr(Y)!")
    # check ace
    if ace_list:
        for position in ace_list:
            if fasta[position] == 'K':
                fasta[position] = 'acK'
            else:
                raise ValueError(f"Position {position} is not Lys(K)!")
    # check met
    if met_list:
        for position in met_list:
            if fasta[position] == 'R':
                fasta[position] = 'Rm'
            else:
                raise ValueError(f"Position {position} is not Arg(R)!")
    return fasta

# -------------------
# genDCD post simulation analysis
# -------------------
def genDCD(residues, name, prot, temp, n_chains):
    path = '{:s}/{:d}'.format(name, temp)
    top = md.Topology()
    for _ in range(n_chains):
        chain = top.add_chain()
        for resname in prot.fasta:
            residue = top.add_residue(residues.loc[resname, 'three'], chain)
            top.add_atom(residues.loc[resname, 'three'],
                         element=md.element.carbon, residue=residue)
        for i in range(chain.n_atoms-1):
            top.add_bond(chain.atom(i), chain.atom(i+1))

    t = md.load(path+'/{:s}.dcd'.format(name), top=top)
    t.xyz *= 10
    t.unitcell_lengths *= 10
    lz = t.unitcell_lengths[0, 2]
    edges = np.arange(-lz/2., lz/2., 1)
    dz = (edges[1]-edges[0])/2.
    z = edges[:-1]+dz
    h = np.apply_along_axis(lambda a: np.histogram(a, bins=edges)[0], 1, t.xyz[:, :, 2])
    zmid = np.apply_along_axis(lambda a: z[a.argmax()], 1, h)
    indices = np.argmin(np.abs(t.xyz[:, :, 2]-zmid[:, np.newaxis]), axis=1)
    t[0].save_pdb(path+'/top.pdb')
    t.save_dcd(path+'/traj4.dcd')

    u = MDAnalysis.Universe(path+'/top.pdb', path+'/traj4.dcd')
    ag = u.atoms
    with MDAnalysis.Writer(path+'/traj3.dcd', ag.n_atoms) as W:
        for ts, ndx in zip(u.trajectory, indices):
            ts = transformations.unwrap(ag)(ts)
            ts = transformations.center_in_box(
                u.select_atoms('index {:d}'.format(ndx)), center='geometry')(ts)
            ts = transformations.wrap(ag)(ts)
            W.write(ag)

    t = md.load(path+'/traj3.dcd', top=path+'/top.pdb')
    edges = np.arange(0, lz, 1)
    dz = (edges[1]-edges[0])/2.
    z = edges[:-1]+dz
    h = np.apply_along_axis(lambda a: np.histogram(a, bins=edges)[0], 1, t.xyz[:, :, 2])
    h = np.mean(h[:120], axis=0)
    maxoverlap = np.apply_along_axis(lambda a: np.correlate(h, np.histogram(a,
                bins=edges)[0], 'full').argmax()-h.size+dz, 1, t.xyz[:, :, 2])

    u = MDAnalysis.Universe(path+'/top.pdb', path+'/traj3.dcd')
    ag = u.atoms
    with MDAnalysis.Writer(path+'/traj2.dcd', ag.n_atoms) as W:
        for ts, mo in zip(u.trajectory, maxoverlap):
            ts = transformations.unwrap(ag)(ts)
            ts = transformations.translate([0, 0, mo*10])(ts)
            ts = transformations.wrap(ag)(ts)
            W.write(ag)

    t = md.load(path+'/traj2.dcd', top=path+'/top.pdb')
    h = np.apply_along_axis(lambda a: np.histogram(a, bins=edges)[0], 1, t.xyz[:, :, 2])
    zmid = np.apply_along_axis(lambda a: z[a > np.quantile(a, .98)].mean(), 1, h)
    indices = np.argmin(np.abs(t.xyz[:, :, 2]-zmid[:, np.newaxis]), axis=1)

    u = MDAnalysis.Universe(path+'/top.pdb', path+'/traj2.dcd')
    ag = u.atoms
    with MDAnalysis.Writer(path+'/traj1.dcd', ag.n_atoms) as W:
        for ts, ndx in zip(u.trajectory, indices):
            ts = transformations.unwrap(ag)(ts)
            ts = transformations.center_in_box(
                u.select_atoms('index {:d}'.format(ndx)), center='geometry')(ts)
            ts = transformations.wrap(ag)(ts)
            W.write(ag)

    t = md.load(path+'/traj1.dcd', top=path+'/top.pdb')
    h = np.apply_along_axis(lambda a: np.histogram(a, bins=edges)[0], 1, t.xyz[:, :, 2])
    h = np.mean(h[120:], axis=0)
    maxoverlap = np.apply_along_axis(lambda a: np.correlate(h, np.histogram(a,
                bins=edges)[0], 'full').argmax()-h.size+dz, 1, t.xyz[:, :, 2])

    u = MDAnalysis.Universe(path+'/top.pdb', path+'/traj1.dcd')
    ag = u.atoms
    with MDAnalysis.Writer(path+'/traj.dcd', ag.n_atoms) as W:
        for ts, mo in zip(u.trajectory, maxoverlap):
            ts = transformations.unwrap(ag)(ts)
            ts = transformations.translate([0, 0, mo*10])(ts)
            ts = transformations.wrap(ag)(ts)
            W.write(ag)

    t = md.load(path+'/traj.dcd', top=path+'/top.pdb')
    h = np.apply_along_axis(lambda a: np.histogram(a, bins=edges)[0], 1, t.xyz[:, :, 2])
    np.save('{:s}_{:d}.npy'.format(name, temp), h, allow_pickle=False)
    os.remove(path+'/traj1.dcd')
    os.remove(path+'/traj2.dcd')
    os.remove(path+'/traj3.dcd')
    os.remove(path+'/traj4.dcd')
    t.xyz /= 10
    t.unitcell_lengths /= 10
    t[0].save_pdb(path+'/top.pdb')
    t.save_dcd(path+'/traj.dcd')

# -------------------
# Define Forces
# -------------------
def genParamsLJ(df,name,prot):
    fasta = prot.fasta.copy()
    r = df.copy()
    r.loc['X'] = r.loc[fasta[0]]
    r.loc['Z'] = r.loc[fasta[-1]]
    r.loc['X','MW'] += 2
    r.loc['Z','MW'] += 16
    fasta[0] = 'X'
    fasta[-1] = 'Z'
    types = list(np.unique(fasta))
    MWs = [r.loc[a,'MW'] for a in types]
    lj_eps = prot.eps_factor*4.184
    return lj_eps, fasta, types, MWs

def genParamsDH(df,name,prot,temp):
    kT = 8.3145*temp*1e-3
    fasta = prot.fasta.copy()
    r = df.copy()
    # Set the charge on HIS based on the pH of the protein solution
    #cal H
    r.loc['H','q'] = 0.5 / ( 1 + 10**(prot.pH-6) )
    #cal Sp
    pka1 = 2.0
    pka2 = 6.5
    a0 = 1.0 / (1+ math.pow(10,prot.pH-pka1) + math.pow(10,2*prot.pH-pka1-pka2))
    a1 = math.pow(10,prot.pH-pka1) / (1+ math.pow(10,prot.pH-pka1) + math.pow(10,2*prot.pH-pka1-pka2))
    a2 = math.pow(10,2*prot.pH-pka1-pka2) / (1 + math.pow(10,prot.pH-pka1)+ math.pow(10,2*prot.pH-pka1-pka2))
    r.loc['Sp','q'] = 0 * a0 - a1 - 2*a2
    #cal Tp
    pka1 = 2.1
    pka2 = 7.0
    a0 = 1.0 / (1+ math.pow(10,prot.pH-pka1) + math.pow(10,2*prot.pH-pka1-pka2))
    a1 = math.pow(10,prot.pH-pka1) / (1+ math.pow(10,prot.pH-pka1) + math.pow(10,2*prot.pH-pka1-pka2))
    a2 = math.pow(10,2*prot.pH-pka1-pka2) / (1 + math.pow(10,prot.pH-pka1)+ math.pow(10,2*prot.pH-pka1-pka2))
    r.loc['Tp','q'] = 0 * a0 - a1 - 2*a2
    #cal Yp
    pka1 = 2.0
    pka2 = 6.0
    a0 = 1.0 / (1+ math.pow(10,prot.pH-pka1) + math.pow(10,2*prot.pH-pka1-pka2))
    a1 = math.pow(10,prot.pH-pka1) / (1+ math.pow(10,prot.pH-pka1) + math.pow(10,2*prot.pH-pka1-pka2))
    a2 = math.pow(10,2*prot.pH-pka1-pka2) / (1 + math.pow(10,prot.pH-pka1)+ math.pow(10,2*prot.pH-pka1-pka2))
    r.loc['Yp','q'] = 0 * a0 - a1 - 2*a2

    r.loc['X'] = r.loc[fasta[0]]
    r.loc['Z'] = r.loc[fasta[-1]]
    fasta[0] = 'X'
    fasta[-1] = 'Z'
    r.loc['X','q'] = r.loc[prot.fasta[0],'q'] + 1.
    r.loc['Z','q'] = r.loc[prot.fasta[-1],'q'] - 1.
    # Calculate the prefactor for the Yukawa potential
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(temp)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/kT
    yukawa_eps = [r.loc[a].q*np.sqrt(lB*kT) for a in fasta]
    # Calculate the inverse of the Debye length
    yukawa_kappa = np.sqrt(8*np.pi*lB*prot.ionic*6.022/10)
    return yukawa_eps, yukawa_kappa


# -------------------
# Simulation
# -------------------
def simulate(all_parms,name,prot,temp,cutoff):
    residues = all_parms['residues'].set_index('one')

    lj_eps, fasta, types, MWs = genParamsLJ(residues,name,prot)
    _, fasta, _, MWs = genParamsLJ(residues,name,prot)
    yukawa_eps, yukawa_kappa = genParamsDH(residues,name,prot,temp)

    N = len(fasta)

    # set parameters
    L = 15.
    margin = 2
    if N > 350:
        L = 25.
        Lz = 300.
        margin = 8
        Nsteps = int(2e7)
    elif N > 200:
        L = 17.
        Lz = 300.
        margin = 4
        Nsteps = int(6e7)
    else:
        Lz = 10*L

    system = openmm.System()

    # set box vectors
    a = unit.Quantity(np.zeros([3]), unit.nanometers)
    a[0] = L * unit.nanometers
    b = unit.Quantity(np.zeros([3]), unit.nanometers)
    b[1] = L * unit.nanometers
    c = unit.Quantity(np.zeros([3]), unit.nanometers)
    c[2] = Lz * unit.nanometers
    system.setDefaultPeriodicBoxVectors(a, b, c)
    
    # initial config
    xy = np.empty(0)
    xy = np.append(xy,np.random.rand(2)*(L-margin)-(L-margin)/2).reshape((-1,2))
    for x,y in np.random.rand(1000,2)*(L-margin)-(L-margin)/2:
        x1 = x-L if x>0 else x+L
        y1 = y-L if y>0 else y+L
        if np.all(np.linalg.norm(xy-[x,y],axis=1)>.7):
            if np.all(np.linalg.norm(xy-[x1,y],axis=1)>.7):
                if np.all(np.linalg.norm(xy-[x,y1],axis=1)>.7):
                    xy = np.append(xy,[x,y]).reshape((-1,2))
        if xy.shape[0] == 100:
            break

    n_chains = xy.shape[0]

    pdb_file = name+'/{:d}/top.pdb'.format(temp)

    if os.path.isfile(pdb_file):
        pdb = app.pdbfile.PDBFile(pdb_file)
    else:
        top = md.Topology()
        pos = []
        for x,y in xy:
            chain = top.add_chain()
            pos.append([[x,y,Lz/2+(i-N/2.)*.38] for i in range(N)])
            for resname in fasta:
                # rename to norm resi name for mdtraj to know
                if resname=='Sp':
                    resname='SEP'
                elif resname=='Tp':
                    resname='TPO'
                elif resname=='Yp':
                    resname='PTR'
                elif resname=='acK':
                    resname='ALY'
                elif resname=='Rm':
                    resname='DA2'
                residue = top.add_residue(resname, chain)
                top.add_atom(resname, element=md.element.carbon, residue=residue)
            for i in range(chain.n_atoms-1):
                top.add_bond(chain.atom(i),chain.atom(i+1))
        md.Trajectory(np.array(pos).reshape(n_chains*N,3), top, 0, [L,L,Lz], [90,90,90]).save_pdb(pdb_file)
        pdb = app.pdbfile.PDBFile(pdb_file)

    for _ in range(n_chains):
        system.addParticle((residues.loc[prot.fasta[0]].MW+2)*unit.amu)
        for a in prot.fasta[1:-1]:
            system.addParticle(residues.loc[a].MW*unit.amu) 
        system.addParticle((residues.loc[prot.fasta[-1]].MW+16)*unit.amu)

    hb = openmm.openmm.HarmonicBondForce()

    energy_expression = 'select(step(r-2^(1/6)*s),4*eps*l*((s/r)^12-(s/r)^6-shift),4*eps*((s/r)^12-(s/r)^6-l*shift)+eps*(1-l))'
    ah = openmm.openmm.CustomNonbondedForce(energy_expression+'; shift=(s/rc)^12-(s/rc)^6; s=s_unit*sTable(type1,type2); l=l_unit*lTable(type1,type2)')

    ah.addGlobalParameter('eps',lj_eps*unit.kilojoules_per_mole)
    ah.addGlobalParameter('rc',cutoff*unit.nanometer)
    ah.addGlobalParameter('s_unit',1*unit.nanometer)
    ah.addGlobalParameter('l_unit',1*unit.dimensionless)
    ah.addPerParticleParameter('type')

    ntypes = 25
    
    sTable = all_parms['s'].to_numpy().flatten()
    sTable = openmm.openmm.Discrete2DFunction(ntypes, ntypes, sTable)
    ah.addTabulatedFunction("sTable", sTable)
    
    muTable = all_parms['l'].to_numpy().flatten()
    muTable = openmm.openmm.Discrete2DFunction(ntypes, ntypes, muTable)
    ah.addTabulatedFunction("lTable", muTable)
 
    yu = openmm.openmm.CustomNonbondedForce('q*(exp(-kappa*r)/r-shift); q=q1*q2')
    yu.addGlobalParameter('kappa',yukawa_kappa/unit.nanometer)
    yu.addGlobalParameter('shift',np.exp(-yukawa_kappa*4.0)/4.0/unit.nanometer)
    yu.addPerParticleParameter('q')

    for j in range(n_chains):
        begin = j*N
        end = j*N+N
       
        for a,e in zip(prot.fasta,yukawa_eps):
            yu.addParticle([e*unit.nanometer*unit.kilojoules_per_mole])
            ah.addParticle([residues.loc[a].type])

        for i in range(begin,end-1):
            hb.addBond(i, i+1, 0.38*unit.nanometer, 8033.0*unit.kilojoules_per_mole/(unit.nanometer**2))
            yu.addExclusion(i, i+1)
            ah.addExclusion(i, i+1)

    yu.setForceGroup(0)
    ah.setForceGroup(1)
    yu.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    ah.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    hb.setUsesPeriodicBoundaryConditions(True)
    yu.setCutoffDistance(4*unit.nanometer)
    ah.setCutoffDistance(cutoff*unit.nanometer)

    system.addForce(hb)
    system.addForce(yu)
    system.addForce(ah)

    print(ah.usesPeriodicBoundaryConditions())
    print(yu.usesPeriodicBoundaryConditions())
    print(hb.usesPeriodicBoundaryConditions())

    integrator = openmm.openmm.LangevinIntegrator(temp*unit.kelvin,0.01/unit.picosecond,0.01*unit.picosecond)

    print(integrator.getFriction(),integrator.getTemperature())

    platform = openmm.Platform.getPlatformByName('CUDA')

    simulation = app.simulation.Simulation(pdb.topology, system, integrator, platform, dict(CudaPrecision='mixed'))

    check_point = name+'/{:d}/restart.chk'.format(temp)

    if os.path.isfile(check_point):
        print('Reading check point file')
        simulation.loadCheckpoint(check_point)
        simulation.reporters.append(app.dcdreporter.DCDReporter(name+'/{:d}/{:s}.dcd'.format(temp,name),int(1e5),append=True))
    else:
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(app.dcdreporter.DCDReporter(name+'/{:d}/{:s}.dcd'.format(temp,name),int(1e5)))

    simulation.reporters.append(app.statedatareporter.StateDataReporter('{:s}_{:d}.log'.format(name,temp),50000,
             step=True,speed=True,elapsedTime=True,separator='\t'))

    simulation.step(args.time*100000) #ns

    simulation.saveCheckpoint(check_point)

    genDCD(residues,name,prot,temp,n_chains)

# -------------------
# Main Function
# -------------------
if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--name', type=str, required=True, help="Protein name (must exactly match an ID in the FASTA file, without '>').")
    parser.add_argument('--fasta_file', type=str, required=True, help="Input protein FASTA file (can contain multiple sequences).")
    parser.add_argument('--parm_file', type=str, default="calvados_ptm_parms.xlsx", help="Forcefield parameter file.")
    parser.add_argument('--phos', type=str, default="", help="Comma-separated list of 0-based positions for phosphorylation PTMs (S/T/Y).")
    parser.add_argument('--ace', type=str, default="", help="Comma-separated list of 0-based positions for acetylation PTMs (K).")
    parser.add_argument('--met', type=str, default="", help="Comma-separated list of 0-based positions for methylation PTMs (R).")
    parser.add_argument('--time', type=int, default=100, help="The simulation time (ns)")
    parser.add_argument('--cutoff', type=float, default=2.0)
    parser.add_argument('--pH', type=float, default=7.0)
    parser.add_argument('--eps_factor', type=float, default=0.2)
    parser.add_argument('--ionic', type=float, default=0.15)
    parser.add_argument('--temp', type=int, default=300)
    args = parser.parse_args()

    # 读取 FASTA 文件成字典
    fasta_dict = {}
    current_name = None
    seq_lines = []
    with open(args.fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # 保存上一个蛋白
                if current_name and seq_lines:
                    fasta_dict[current_name] = "".join(seq_lines)
                current_name = line[1:].strip()  # 移除 '>'
                seq_lines = []
            else:
                seq_lines.append(line)
        # 保存最后一个蛋白
        if current_name and seq_lines:
            fasta_dict[current_name] = "".join(seq_lines)

    # 检查是否有目标蛋白
    if args.name not in fasta_dict:
        raise ValueError(f"Could not find protein named '{args.name}' in the given fasta file, please chech the --name parameter")

    seq = fasta_dict[args.name]
    print(f"protein {args.name} sequence length: {len(seq)}")

    # 转换 PTM 位点
    phos_list = [int(x) for x in args.phos.split(",") if x.strip().isdigit()]
    ace_list = [int(x) for x in args.ace.split(",") if x.strip().isdigit()]
    met_list = [int(x) for x in args.met.split(",") if x.strip().isdigit()]
    print("phosphorylation sites:", phos_list)
    print("acetylation sites:", ace_list)
    print("methylation:", met_list)

    ptm_seq = trans_ptm(seq, phos_list=phos_list, ace_list=ace_list, met_list=met_list)
    print("PTM protein sequence:", ptm_seq)

    # 加载 calvados 参数表
    parms_file = args.parm_file
    sheets_config = {
        "l": {"header": 0, "index_col": 0},
        "s": {"header": 0, "index_col": 0},
        "residues": {"header": 0, "index_col": None}
    }
    all_parms = {}
    for sheet_name in sheets_config:
        params = sheets_config[sheet_name]
        df = pd.read_excel(parms_file, sheet_name=sheet_name, **params)
        all_parms[sheet_name] = df

    # 构造 prot 数据并运行模拟
    proteins = pd.DataFrame(
        [[args.eps_factor, args.pH, args.ionic, args.cutoff, ptm_seq]],
        columns=["eps_factor", "pH", "ionic", "cutoff", "fasta"],
        index=[args.name]
    )
    simulate(all_parms, args.name, proteins.loc[args.name], args.temp, args.cutoff)
