import os
import json
import argparse
import tempfile
import time
import hashlib
import numpy as np
from collections import Counter
from Bio.PDB import PDBParser, MMCIFParser, Select, PDBIO, PPBuilder
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, MACCSkeys
from rdkit.Chem.rdchem import HybridizationType
try:
    from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
    _HAS_MORGAN_GEN = True
except Exception:
    _HAS_MORGAN_GEN = False


AA_LETTERS = 'ACDEFGHIKLMNPQRSTVWY'
AA_SET = set(AA_LETTERS)


TIMES = {}
def tick():
    return time.perf_counter()
def tock(key, t0):
    TIMES[key] = TIMES.get(key, 0.0) + (time.perf_counter() - t0)


class ProteinOnlySelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] == ' '


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def three_to_one(aa):
    d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
    return d[aa.upper()]


def read_ligand_from_sdf(sdf_path, add_h=True, gen_conf_if_missing=True):
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError('failed to read ligand from sdf')
    if add_h:
        mol = Chem.AddHs(mol, addCoords=True)
    if mol.GetNumConformers() == 0 and gen_conf_if_missing:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    return mol


def extract_protein_and_ligand_from_complex_pdb(pdb_path, out_dir):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', pdb_path)
    io = PDBIO()
    prot_path = os.path.join(out_dir, 'protein_only.pdb')
    io.set_structure(structure)
    io.save(prot_path, ProteinOnlySelect())

    lig_records = []
    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag = residue.id[0].strip()
                resname = residue.get_resname().strip()
                if hetflag != '' and resname not in ['HOH', 'WAT']:
                    lig_records.append((len(list(residue.get_atoms())), model.id, chain.id, residue.id, resname))
    if len(lig_records) == 0:
        raise ValueError('no ligand hetatm found in complex pdb')

    lig_records.sort(key=lambda x: x[0], reverse=True)
    _, model_id, chain_id, resid, resname = lig_records[0]

    lig_atoms = []
    for model in structure:
        if model.id != model_id:
            continue
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                if residue.id == resid:
                    for atom in residue.get_atoms():
                        lig_atoms.append(atom)
                    break

    lig_pdb_path = os.path.join(out_dir, 'ligand_only.pdb')
    with open(lig_pdb_path, 'w') as f:
        f.write('MODEL\n')
        for i, atom in enumerate(lig_atoms, start=1):
            name = atom.get_name().rjust(4)
            rname = resname.rjust(3)
            ch = chain_id
            resseq = resid[1]
            x, y, z = atom.coord
            occ = atom.get_occupancy() if atom.get_occupancy() is not None else 1.0
            b = atom.get_bfactor() if atom.get_bfactor() is not None else 0.0
            elem = atom.element.strip().rjust(2) if atom.element is not None else ' C'
            f.write(f'HETATM{i:5d} {name} {rname} {ch}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}          {elem}\n')
        f.write('ENDMDL\n')

    lig_mol = Chem.MolFromPDBFile(lig_pdb_path, removeHs=False)
    if lig_mol is None:
        lig_mol = Chem.MolFromPDBFile(lig_pdb_path, removeHs=True)
        if lig_mol is None:
            raise ValueError('failed to parse ligand from complex pdb')
        lig_mol = Chem.AddHs(lig_mol, addCoords=True)
    if lig_mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(lig_mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(lig_mol, maxIters=200)

    return prot_path, lig_mol


def extract_protein_and_ligand_from_complex_cif(cif_path, out_dir):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('complex', cif_path)
    io = PDBIO()
    prot_path = os.path.join(out_dir, 'protein_only.pdb')
    io.set_structure(structure)
    io.save(prot_path, ProteinOnlySelect())

    lig_records = []
    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag = residue.id[0].strip()
                resname = residue.get_resname().strip()
                if hetflag != '' and resname not in ['HOH', 'WAT']:
                    lig_records.append((len(list(residue.get_atoms())), model.id, chain.id, residue.id, resname))
    if len(lig_records) == 0:
        raise ValueError('no ligand hetatm found in complex cif')
    lig_records.sort(key=lambda x: x[0], reverse=True)
    _, model_id, chain_id, resid, resname = lig_records[0]

    lig_atoms = []
    for model in structure:
        if model.id != model_id:
            continue
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                if residue.id == resid:
                    for atom in residue.get_atoms():
                        lig_atoms.append(atom)
                    break

    lig_pdb_path = os.path.join(out_dir, 'ligand_only.pdb')
    with open(lig_pdb_path, 'w') as f:
        f.write('MODEL\n')
        for i, atom in enumerate(lig_atoms, start=1):
            name = atom.get_name().rjust(4)
            rname = resname.rjust(3)
            ch = chain_id
            resseq = resid[1]
            x, y, z = atom.coord
            occ = atom.get_occupancy() if atom.get_occupancy() is not None else 1.0
            b = atom.get_bfactor() if atom.get_bfactor() is not None else 0.0
            elem = (atom.element.strip() if atom.element else 'C').rjust(2)
            f.write(f'HETATM{i:5d} {name} {rname} {ch}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}          {elem}\n')
        f.write('ENDMDL\n')

    lig_mol = Chem.MolFromPDBFile(lig_pdb_path, removeHs=False)
    if lig_mol is None:
        lig_mol = Chem.MolFromPDBFile(lig_pdb_path, removeHs=True)
        if lig_mol is None:
            raise ValueError('failed to parse ligand from complex cif')
        lig_mol = Chem.AddHs(lig_mol, addCoords=True)
    if lig_mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(lig_mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(lig_mol, maxIters=200)
    return prot_path, lig_mol


def read_protein_sequence_and_ca(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('prot', pdb_path)
    ppb = PPBuilder()
    seq = ''
    ca_coords = []
    b_factors = []
    residues = []
    for pp in ppb.build_peptides(structure):
        s = pp.get_sequence()
        seq += str(s)
        for res in pp:
            ca = res['CA'] if 'CA' in res else None
            if ca is not None:
                ca_coords.append(ca.get_coord())
                b_factors.append(ca.get_bfactor())
                try:
                    aa = three_to_one(res.get_resname())
                except KeyError:
                    aa = 'X'
                residues.append(aa)
    if len(seq) == 0 and len(residues) > 0:
        seq = ''.join([a if a in AA_SET else 'X' for a in residues])
    return seq, np.array(ca_coords, dtype=np.float32), np.array(b_factors, dtype=np.float32)


def seq_one_hot(seq):
    mat = np.zeros((len(seq), len(AA_LETTERS)), dtype=np.float32)
    idx = {aa:i for i, aa in enumerate(AA_LETTERS)}
    for i, a in enumerate(seq):
        if a in idx:
            mat[i, idx[a]] = 1.0
    return mat


def seq_aac(seq):
    cnt = Counter(a for a in seq if a in AA_SET)
    total = sum(cnt.values()) if sum(cnt.values()) > 0 else 1
    vec = np.array([cnt.get(a, 0) / total for a in AA_LETTERS], dtype=np.float32)
    return vec


def ligand_morgan(mol, n_bits=2048, radius=2):
    if _HAS_MORGAN_GEN:
        gen = GetMorganGenerator(radius=radius, fpSize=n_bits)
        fp = gen.GetFingerprint(mol)
    else:
        fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=np.int8)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def ligand_maccs(mol):
    fp = MACCSkeys.GenMACCSKeys(mol)
    arr = np.zeros((len(fp),), dtype=np.int8)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def atom_feature_vector(atom):
    feats = []
    feats.append(atom.GetAtomicNum())
    feats.append(atom.GetTotalDegree())
    feats.append(atom.GetFormalCharge())
    feats.append(int(atom.GetIsAromatic()))
    hyb = atom.GetHybridization()
    feats.extend([
        int(hyb == HybridizationType.SP),
        int(hyb == HybridizationType.SP2),
        int(hyb == HybridizationType.SP3)
    ])
    feats.append(atom.GetTotalNumHs(includeNeighbors=True))
    chiral_tag = atom.GetChiralTag()
    feats.extend([
        int(chiral_tag.name == 'CHI_TETRAHEDRAL_CW'),
        int(chiral_tag.name == 'CHI_TETRAHEDRAL_CCW')
    ])
    return feats


def bond_feature_vector(bond):
    btype = bond.GetBondType()
    feats = [
        int(btype == Chem.rdchem.BondType.SINGLE),
        int(btype == Chem.rdchem.BondType.DOUBLE),
        int(btype == Chem.rdchem.BondType.TRIPLE),
        int(btype == Chem.rdchem.BondType.AROMATIC),
        int(bond.GetIsConjugated()),
        int(bond.IsInRing())
    ]
    return feats


def build_ligand_graph(mol):
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    n = mol.GetNumAtoms()
    x = np.array([atom_feature_vector(mol.GetAtomWithIdx(i)) for i in range(n)], dtype=np.float32)
    edge_index = []
    edge_attr = []
    for b in mol.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        bf = bond_feature_vector(b)
        edge_index.append([i, j])
        edge_attr.append(bf)
        edge_index.append([j, i])
        edge_attr.append(bf)
    edge_index = np.array(edge_index, dtype=np.int64).T if len(edge_index) > 0 else np.zeros((2, 0), dtype=np.int64)
    edge_attr = np.array(edge_attr, dtype=np.float32) if len(edge_attr) > 0 else np.zeros((0, 6), dtype=np.float32)
    coords = mol.GetConformer().GetPositions().astype(np.float32)
    return x, edge_index, edge_attr, coords


def residue_graph_from_ca(ca_coords, b_factors, seq, radius=8.0):
    n = ca_coords.shape[0]
    if n == 0:
        return np.zeros((0, len(AA_LETTERS) + 1), dtype=np.float32), np.zeros((2, 0), dtype=np.int64), ca_coords
    aa_onehot = seq_one_hot(seq[:n])
    bf = b_factors[:n].reshape(-1, 1).astype(np.float32) if b_factors.size >= n else np.zeros((n, 1), dtype=np.float32)
    x = np.concatenate([aa_onehot, bf], axis=1)
    diff = ca_coords[:, None, :] - ca_coords[None, :, :]
    D = np.sqrt((diff ** 2).sum(-1)).astype(np.float32)
    edges = np.argwhere((D < radius) & (D > 0))
    if edges.size == 0:
        edge_index = np.zeros((2, 0), dtype=np.int64)
    else:
        edge_index = edges.T.astype(np.int64)
    return x, edge_index, ca_coords.astype(np.float32)


def simple_pocket_grid(prot_pdb_path, lig_mol, grid_size=24.0, spacing=1.0):
    conf = lig_mol.GetConformer()
    lig_xyz = conf.GetPositions()
    center = lig_xyz.mean(axis=0)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('prot', prot_pdb_path)
    atoms = []
    for atom in structure.get_atoms():
        el = atom.element.strip() if atom.element is not None else 'C'
        atoms.append((el, atom.coord))
    if len(atoms) == 0:
        return np.zeros((0, 0, 0, 0), dtype=np.float32), center.astype(np.float32)

    half = grid_size / 2.0
    xs = np.arange(center[0] - half, center[0] + half, spacing)
    ys = np.arange(center[1] - half, center[1] + half, spacing)
    zs = np.arange(center[2] - half, center[2] + half, spacing)
    nx, ny, nz = len(xs), len(ys), len(zs)

    channels = ['C', 'N', 'O', 'S', 'Hal', 'H', 'Metal']
    ch_idx = {c:i for i, c in enumerate(channels)}
    grid = np.zeros((len(channels), nx, ny, nz), dtype=np.float32)

    def elem_channel(e):
        e = e.upper()
        if e in ['F', 'CL', 'BR', 'I']:
            return 'Hal'
        if e in ['NA', 'K', 'CA', 'MG', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN']:
            return 'Metal'
        if e in ['C', 'N', 'O', 'S', 'H']:
            return e
        return 'C'

    for e, coord in atoms:
        ch = elem_channel(e)
        i = int(np.floor((coord[0] - xs[0]) / spacing))
        j = int(np.floor((coord[1] - ys[0]) / spacing))
        k = int(np.floor((coord[2] - zs[0]) / spacing))
        if 0 <= i < nx and 0 <= j < ny and 0 <= k < nz:
            grid[ch_idx[ch], i, j, k] += 1.0

    return grid, center.astype(np.float32)


def build_complex_heterograph(lig_coords, prot_coords, cutoff=6.0):
    if lig_coords.size == 0 or prot_coords.size == 0:
        return np.zeros((2, 0), dtype=np.int64), np.zeros((0, 1), dtype=np.float32)
    L = lig_coords.shape[0]
    P = prot_coords.shape[0]
    diff = lig_coords[:, None, :] - prot_coords[None, :, :]
    D = np.sqrt((diff ** 2).sum(-1)).astype(np.float32)
    i_idx, j_idx = np.where(D < cutoff)
    if i_idx.size == 0:
        return np.zeros((2, 0), dtype=np.int64), np.zeros((0, 1), dtype=np.float32)
    edge_index = np.stack([i_idx, j_idx], axis=0).astype(np.int64)
    edge_attr = D[i_idx, j_idx].reshape(-1, 1)
    return edge_index, edge_attr


def save_npz(path, **arrays):
    ensure_dir(os.path.dirname(path))
    np.savez_compressed(path, **arrays)


def write_fasta(path, header, sequence):
    with open(path, 'w') as f:
        f.write(f'>{header}\n')
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + '\n')


def stable_hash(text):
    return hashlib.sha1(text.encode('utf-8')).hexdigest()[:16]


class GPUTimer:
    def __init__(self, device):
        self.device = device
        self.enabled = False
        try:
            import torch
            self.enabled = (device == 'cuda' and torch.cuda.is_available())
            self.torch = torch
        except Exception:
            self.enabled = False
            self.torch = None
        self.times = {}
    def start(self, key):
        if self.enabled:
            self.torch.cuda.synchronize()
        self.times[key] = [-1.0, -1.0]
        self.times[key][0] = time.perf_counter()
    def stop(self, key):
        if self.enabled:
            self.torch.cuda.synchronize()
        self.times[key][1] = time.perf_counter()
    def as_seconds(self, prefix=''):
        out = {}
        for k, (s, e) in self.times.items():
            if s >= 0 and e >= 0:
                out[(prefix + k) if prefix else k] = round(e - s, 6)
        return out


def plm_embed_esm2(seq, model_name, cache_dir, device, fp16, batch_size):
    ensure_dir(cache_dir)
    h = stable_hash(model_name + '|' + seq)
    cache_path = os.path.join(cache_dir, f'esm2_{h}.npy')
    if os.path.exists(cache_path):
        emb = np.load(cache_path)
        return emb, {'plm_cached': True, 'plm_backend': 'cache'}

    # try esm backend first (ESM3-style callable, then legacy load_model_and_alphabet)
    try:
        import torch
        from esm import pretrained
        dev = device if device != 'auto' else ('cuda' if torch.cuda.is_available() else 'cpu')
        timer = GPUTimer(dev)

        timer.start('load_model')
        model = None
        alphabet = None
        fn = getattr(pretrained, model_name, None)
        if callable(fn):
            model, alphabet = fn()
        elif hasattr(pretrained, 'load_model_and_alphabet'):
            model, alphabet = pretrained.load_model_and_alphabet(model_name)
        else:
            raise RuntimeError('esm.pretrained does not expose the requested model')
        model.eval().to(dev)
        if fp16 and dev == 'cuda':
            model.half()
        batch_converter = alphabet.get_batch_converter()
        timer.stop('load_model')

        data = [('0', seq)]
        timer.start('tokenize')
        _, _, toks = batch_converter(data)
        timer.stop('tokenize')

        timer.start('to_device')
        toks = toks.to(dev)
        timer.stop('to_device')

        with torch.no_grad():
            timer.start('forward')
            out = model(toks, repr_layers=[model.num_layers], return_contacts=False)
            timer.stop('forward')
            timer.start('pool')
            reps = out['representations'][model.num_layers]  # [1, L, D]
            mask = (toks != alphabet.padding_idx)            # [1, L]
            pooled = (reps * mask.unsqueeze(-1)).sum(dim=1) / mask.sum(dim=1, keepdim=True).clamp(min=1)
            emb = pooled.float().cpu().numpy()
            timer.stop('pool')

        np.save(cache_path, emb)
        t = timer.as_seconds(prefix='plm_')
        t.update({'plm_cached': False, 'plm_backend': 'esm', 'plm_device': dev, 'plm_fp16': bool(fp16 and dev == 'cuda'),
                  'plm_batch_size': int(batch_size), 'plm_model': model_name})
        return emb, t

    except Exception as esm_err:
        # fallback to huggingface
        try:
            import torch
            from transformers import EsmTokenizer, EsmModel
            dev = device if device != 'auto' else ('cuda' if torch.cuda.is_available() else 'cpu')
            timer = GPUTimer(dev)

            hf_name = model_name if model_name.startswith('facebook/') else f'facebook/{model_name}'
            timer.start('load_model')
            tok = EsmTokenizer.from_pretrained(hf_name, use_fast=False)
            model = EsmModel.from_pretrained(hf_name)
            model.eval().to(dev)
            if fp16 and dev == 'cuda':
                model.half()
            timer.stop('load_model')

            timer.start('tokenize')
            toks = tok(seq, return_tensors='pt', add_special_tokens=True)
            timer.stop('tokenize')

            timer.start('to_device')
            toks = {k: v.to(dev) for k, v in toks.items()}
            timer.stop('to_device')

            with torch.no_grad():
                timer.start('forward')
                out = model(**toks, output_hidden_states=False)
                timer.stop('forward')
                timer.start('pool')
                reps = out.last_hidden_state
                mask = toks['attention_mask']
                pooled = (reps * mask.unsqueeze(-1)).sum(dim=1) / mask.sum(dim=1, keepdim=True).clamp(min=1)
                emb = pooled.float().cpu().numpy()
                timer.stop('pool')

            np.save(cache_path, emb)
            t = timer.as_seconds(prefix='plm_')
            t.update({'plm_cached': False, 'plm_backend': 'hf', 'plm_device': dev, 'plm_fp16': bool(fp16 and dev == 'cuda'),
                      'plm_batch_size': int(batch_size), 'plm_model': hf_name})
            return emb, t

        except Exception as hf_err:
            raise RuntimeError('failed to load ESM2 with esm or transformers') from hf_err


def main():
    parser = argparse.ArgumentParser(description='Produce block-ready inputs with timings for affinity models')
    parser.add_argument('--ligand_sdf', type=str, default=None)
    parser.add_argument('--protein_pdb', type=str, default=None)
    parser.add_argument('--complex_pdb', type=str, default=None)
    parser.add_argument('--outdir', type=str, required=True)
    parser.add_argument('--grid_size', type=float, default=24.0)
    parser.add_argument('--grid_spacing', type=float, default=1.0)
    parser.add_argument('--prot_radius', type=float, default=8.0)
    parser.add_argument('--cxn_cutoff', type=float, default=6.0)
    parser.add_argument('--protein_plm', type=str, default=None)
    parser.add_argument('--plm_device', type=str, default='auto')
    parser.add_argument('--plm_fp16', action='store_true')
    parser.add_argument('--plm_batch_size', type=int, default=4)
    parser.add_argument('--plm_cache_dir', type=str, default=None)
    args = parser.parse_args()

    ensure_dir(args.outdir)
    lig_dir = os.path.join(args.outdir, 'ligand')
    prot_dir = os.path.join(args.outdir, 'protein')
    joint_dir = os.path.join(args.outdir, 'joint')
    ensure_dir(lig_dir)
    ensure_dir(prot_dir)
    ensure_dir(joint_dir)
    plm_cache_dir = args.plm_cache_dir or os.path.join(prot_dir, 'plm_cache')


    t0 = tick()
    if args.complex_pdb:
        with tempfile.TemporaryDirectory() as tmpd:
            protein_pdb_path, ligand_mol = extract_protein_and_ligand_from_complex_pdb(args.complex_pdb, tmpd)
            dest_prot = os.path.join(args.outdir, 'protein_from_complex.pdb')
            if os.path.exists(dest_prot):
                os.remove(dest_prot)
            os.replace(protein_pdb_path, dest_prot)
            protein_pdb_path = dest_prot
    else:
        if not args.ligand_sdf or not args.protein_pdb:
            raise ValueError('provide either --complex_pdb or both --ligand_sdf and --protein_pdb')
        protein_pdb_path = args.protein_pdb
        ligand_mol = read_ligand_from_sdf(args.ligand_sdf)
    tock('load_inputs', t0)

    t0 = tick()
    lig_smiles = Chem.MolToSmiles(Chem.RemoveHs(ligand_mol)) if ligand_mol is not None else ''
    tock('ligand_smiles', t0)

    t0 = tick()
    morgan = ligand_morgan(ligand_mol, n_bits=2048, radius=2)
    maccs = ligand_maccs(ligand_mol)
    fp = np.concatenate([morgan.astype(np.int8), maccs.astype(np.int8)], axis=0)
    fp_path = os.path.join(lig_dir, 'ligand_fp.npy')
    np.save(fp_path, fp)
    tock('ligand_fp', t0)

    t0 = tick()
    x_lig, ei_lig, ea_lig, pos_lig = build_ligand_graph(ligand_mol)
    lig_graph_path = os.path.join(lig_dir, 'ligand_graph.npz')
    save_npz(lig_graph_path, x=x_lig, edge_index=ei_lig, edge_attr=ea_lig, pos=pos_lig)
    tock('ligand_graph', t0)

    t0 = tick()
    seq, ca_coords, b_factors = read_protein_sequence_and_ca(protein_pdb_path)
    onehot_path = os.path.join(prot_dir, 'protein_seq_onehot.npy')
    aac_path = os.path.join(prot_dir, 'protein_aac.npy')
    fasta_path = os.path.join(prot_dir, 'protein_seq.fasta')
    np.save(onehot_path, seq_one_hot(seq))
    np.save(aac_path, seq_aac(seq))
    write_fasta(fasta_path, 'protein', seq)
    tock('protein_sequence', t0)

    t0 = tick()
    x_res, ei_res, pos_res = residue_graph_from_ca(ca_coords, b_factors, seq, radius=args.prot_radius)
    prot_graph_path = os.path.join(prot_dir, 'protein_graph.npz')
    save_npz(prot_graph_path, x=x_res, edge_index=ei_res, pos=pos_res)
    tock('protein_graph', t0)

    t0 = tick()
    grid, center = simple_pocket_grid(protein_pdb_path, ligand_mol, grid_size=args.grid_size, spacing=args.grid_spacing)
    voxel_path = os.path.join(joint_dir, 'pocket_voxel.npy')
    center_path = os.path.join(joint_dir, 'pocket_center.npy')
    np.save(voxel_path, grid)
    np.save(center_path, center)
    tock('pocket_voxel', t0)

    t0 = tick()
    ei_cxn, ea_cxn = build_complex_heterograph(pos_lig, pos_res, cutoff=args.cxn_cutoff)
    cxn_graph_path = os.path.join(joint_dir, 'complex_graph.npz')
    save_npz(cxn_graph_path, lig_pos=pos_lig, prot_pos=pos_res, edge_index=ei_cxn, edge_attr=ea_cxn)
    tock('complex_graph', t0)
    
    plm_info = None
    if args.protein_plm:
        t0 = tick()
        emb, extra = plm_embed_esm2(seq=seq,
                                    model_name=args.protein_plm,
                                    cache_dir=plm_cache_dir,
                                    device=args.plm_device,
                                    fp16=args.plm_fp16,
                                    batch_size=args.plm_batch_size)
        plm_path = os.path.join(prot_dir, f'plm_{args.protein_plm.replace("/", "_")}.npy')
        np.save(plm_path, emb)
        tock('protein_plm_total', t0)
        plm_info = {
            'file': os.path.relpath(plm_path, args.outdir),
            'shape': list(emb.shape),
            'timings': extra
        }
    
    meta = {
        'inputs': {
            'mode': 'complex' if args.complex_pdb else 'separate',
            'protein_pdb': protein_pdb_path,
            'ligand_source': args.complex_pdb if args.complex_pdb else args.ligand_sdf,
            'ligand_smiles': lig_smiles
        },
        'blocks': {
            'ligand_fp_mlp': {'file': os.path.relpath(fp_path, args.outdir), 'shape': list(fp.shape)},
            'ligand_gnn': {
                'file': os.path.relpath(lig_graph_path, args.outdir),
                'x': list(x_lig.shape),
                'edge_index': [2, int(ei_lig.shape[1])],
                'edge_attr': list(ea_lig.shape),
                'pos': list(pos_lig.shape)
            },
            'protein_seq_plm': {'fasta': os.path.relpath(fasta_path, args.outdir)},
            'protein_seq_light': {
                'onehot': os.path.relpath(onehot_path, args.outdir),
                'aac': os.path.relpath(aac_path, args.outdir)
            },
            'protein_gnn': {
                'file': os.path.relpath(prot_graph_path, args.outdir),
                'x': list(x_res.shape),
                'edge_index': [2, int(ei_res.shape[1])],
                'pos': list(pos_res.shape)
            },
            'voxel_3dcnn': {
                'file': os.path.relpath(voxel_path, args.outdir),
                'center': os.path.relpath(center_path, args.outdir),
                'grid_channels': ['C', 'N', 'O', 'S', 'Hal', 'H', 'Metal']
            },
            'hetero_complex_gnn': {
                'file': os.path.relpath(cxn_graph_path, args.outdir),
                'edge_index': [2, int(ei_cxn.shape[1])],
                'edge_attr': list(ea_cxn.shape)
            }
        },
        'timings_seconds': {k: round(v, 6) for k, v in TIMES.items()}
    }
    if plm_info is not None:
        meta['blocks']['protein_plm'] = plm_info

    with open(os.path.join(args.outdir, 'meta.json'), 'w') as f:
        json.dump(meta, f, indent=2)

if __name__ == '__main__':
    main()
