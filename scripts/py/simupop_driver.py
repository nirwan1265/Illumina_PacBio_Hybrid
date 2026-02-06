#!/usr/bin/env python3
"""
SimuPOP driver for R (reticulate) integration.
Creates a population, applies a mating scheme, and exports genotypes.
"""

from __future__ import annotations
from typing import Dict, Any, List
import random

try:
    import simuPOP as sim
except Exception as e:  # pragma: no cover
    raise RuntimeError("simuPOP not available. Activate the 'simitall' env with simupop installed.") from e


def _read_fasta_multi(path: str):
    ids = []
    seqs = []
    with open(path) as fh:
        cur_id = None
        cur_seq = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    ids.append(cur_id)
                    seqs.append("".join(cur_seq).upper())
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            ids.append(cur_id)
            seqs.append("".join(cur_seq).upper())
    return ids, seqs

def _panel_to_loci(seqs, max_loci=None):
    if not seqs:
        return [], [], []
    L = len(seqs[0])
    for s in seqs:
        if len(s) != L:
            raise ValueError("All haplotypes must be same length")
    var_positions = []
    alleles_per_pos = []
    for i in range(L):
        alleles = sorted(set(s[i] for s in seqs))
        if len(alleles) > 1:
            var_positions.append(i + 1)
            alleles_per_pos.append(alleles)
    if max_loci is not None and max_loci > 0 and len(var_positions) > max_loci:
        # downsample evenly
        idxs = [int(i * len(var_positions) / max_loci) for i in range(max_loci)]
        var_positions = [var_positions[i] for i in idxs]
        alleles_per_pos = [alleles_per_pos[i] for i in idxs]
    return var_positions, alleles_per_pos

def _hap_to_geno(hap_seq, var_positions, allele_maps):
    out = []
    for pos, amap in zip(var_positions, allele_maps):
        allele = hap_seq[pos - 1]
        out.append(amap.get(allele, 0))
    return out

def _make_mating_scheme(cfg: Dict[str, Any]) -> sim.MatingScheme:
    scheme = (cfg.get("scheme") or "RandomMating").strip()
    size = cfg.get("offspring")

    # map a subset of common schemes
    mapping = {
        "RandomMating": sim.RandomMating,
        "RandomSelection": sim.RandomSelection,
        "MonogamousMating": sim.MonogamousMating,
        "PolygamousMating": sim.PolygamousMating,
        "SelfMating": sim.SelfMating,
        "HermaphroditicMating": sim.HermaphroditicMating,
        "ControlledRandomMating": sim.ControlledRandomMating,
        "CloneMating": sim.CloneMating,
        "HaplodiploidMating": sim.HaplodiploidMating,
    }
    if scheme not in mapping:
        raise ValueError(f"Unknown mating scheme: {scheme}")

    if size is not None:
        return mapping[scheme](subPopSize=size)
    return mapping[scheme]()


def _ensure_info_fields(pop: sim.Population, fields: List[str]) -> None:
    if not fields:
        return
    existing = set(pop.infoFields())
    to_add = [f for f in fields if f not in existing]
    if to_add:
        pop.addInfoFields(to_add, init=0)


def _assign_ids(pop: sim.Population, field: str = "ind_id") -> None:
    _ensure_info_fields(pop, [field])
    for i, ind in enumerate(pop.individuals()):
        ind.setInfo(i + 1, field)


def _export_vcf(pop: sim.Population, out_prefix: str) -> None:
    vcf_path = f"{out_prefix}.vcf"
    loci = pop.totNumLoci()
    ploidy = pop.ploidy()
    samples = []
    for ind in pop.individuals():
        sid = ind.info("ind_id") if "ind_id" in pop.infoFields() else "NA"
        samples.append(str(sid))
    with open(vcf_path, "w") as vf:
        vf.write("##fileformat=VCFv4.2
")
        vf.write(f"##contig=<ID=chr1,length={loci}>
")
        vf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
")
        vf.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	" + "	".join(samples) + "
")
        # use allele indices, with REF=0 ALT=1 for biallelic
        for l in range(loci):
            row = ["chr1", str(l+1), f"var{l+1}", "0", "1", ".", "PASS", ".", "GT"]
            gts = []
            for ind in pop.individuals():
                alleles = [str(ind.allele(l, ploidy=p)) for p in range(ploidy)]
                if ploidy == 1:
                    gt = alleles[0]
                else:
                    gt = "/".join(alleles[:2])
                gts.append(gt)
            row.extend(gts)
            vf.write("	".join(row) + "
")

def _export_genotypes(pop: sim.Population, out_prefix: str) -> None:
    geno_path = f"{out_prefix}.genotypes.tsv"
    meta_path = f"{out_prefix}.meta.tsv"

    loci = pop.totNumLoci()
    ploidy = pop.ploidy()

    # header
    locus_headers = [f"locus_{i+1}" for i in range(loci)]

    with open(geno_path, "w") as gf:
        gf.write("sample\t" + "\t".join(locus_headers) + "\n")
        for ind in pop.individuals():
            sid = ind.info("ind_id") if "ind_id" in pop.infoFields() else "NA"
            row = []
            for l in range(loci):
                # allele indices for each ploidy
                alleles = [str(ind.allele(l, ploidy=p)) for p in range(ploidy)]
                row.append("|".join(alleles))
            gf.write(str(sid) + "\t" + "\t".join(row) + "\n")

    with open(meta_path, "w") as mf:
        mf.write("sample\n")
        for ind in pop.individuals():
            sid = ind.info("ind_id") if "ind_id" in pop.infoFields() else "NA"
            mf.write(str(sid) + "\n")


def run_simupop(config: Dict[str, Any], out_prefix: str) -> Dict[str, str]:
    pop_cfg = config.get("population", {})
    mating_cfg = config.get("mating", {})
    init_cfg = config.get("init", {})
    preset = config.get("preset")
    gens = int(config.get("generations", 1))

    size = pop_cfg.get("size", 100)
    ploidy = pop_cfg.get("ploidy", 2)
    loci = pop_cfg.get("loci", [100])
    chrom_names = pop_cfg.get("chromNames", [])
    loci_pos = pop_cfg.get("lociPos", [])
    loci_names = pop_cfg.get("lociNames", [])
    allele_names = pop_cfg.get("alleleNames", [])
    info_fields = pop_cfg.get("infoFields", ["ind_id"])

    init_from_fasta = init_cfg.get("from_fasta")
    max_loci = init_cfg.get("max_loci")
    sample_haplotypes = init_cfg.get("sample_haplotypes", True)
    track_hap_ids = init_cfg.get("track_hap_ids", False)

    if init_from_fasta:
        hap_ids, hap_seqs = _read_fasta_multi(init_from_fasta)
        var_positions, alleles_per_pos = _panel_to_loci(hap_seqs, max_loci=max_loci)
        if not var_positions:
            raise ValueError("No variable sites found in haplotype panel")
        loci = [len(var_positions)]
        loci_pos = var_positions
        allele_names = alleles_per_pos
        chrom_names = ["chr1"]

    pop = sim.Population(
        size=size,
        ploidy=ploidy,
        loci=loci,
        chromNames=chrom_names,
        lociPos=loci_pos,
        lociNames=loci_names,
        alleleNames=allele_names,
        infoFields=info_fields,
    )

    _assign_ids(pop, "ind_id")

    # preset mating schemes
    if preset == "F1":
        mating_cfg = {"scheme": "RandomMating", "offspring": pop.popSize()}
        gens = 1
    elif preset == "F2":
        mating_cfg = {"scheme": "SelfMating", "offspring": pop.popSize()}
        gens = 2
    elif preset == "MAGIC":
        mating_cfg = {"scheme": "RandomMating", "offspring": pop.popSize()}
        gens = 5
    elif preset == "NAM":
        mating_cfg = {"scheme": "RandomMating", "offspring": pop.popSize()}
        gens = 3

    if init_from_fasta:
        if track_hap_ids:
            _ensure_info_fields(pop, ["hap1", "hap2"])
        allele_maps = [ {a:i for i,a in enumerate(alleles)} for alleles in allele_names ]
        for ind in pop.individuals():
            if sample_haplotypes:
                hap_idxs = [random.randrange(len(hap_seqs)) for _ in range(ploidy)]
            else:
                hap_idxs = [0 for _ in range(ploidy)]
            geno = []
            for hi in hap_idxs:
                geno.extend(_hap_to_geno(hap_seqs[hi], loci_pos, allele_maps))
            ind.setGenotype(geno)
            if track_hap_ids and ploidy >= 2:
                ind.setInfo(hap_ids[hap_idxs[0]], "hap1")
                ind.setInfo(hap_ids[hap_idxs[1]], "hap2")

    mating_scheme = _make_mating_scheme(mating_cfg)

    sim.Simulator(pop).evolve(matingScheme=mating_scheme, gen=gens)

    _export_genotypes(pop, out_prefix)

    return {
        "genotypes": f"{out_prefix}.genotypes.tsv",
        "meta": f"{out_prefix}.meta.tsv",
    }

