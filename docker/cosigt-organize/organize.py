#!/usr/bin/env python3

"""
organize.py - Create resources directory structure and config.yaml for COSIGT pipeline
Handles both native execution and Docker/Singularity container paths
"""

import os
import sys
import argparse
import yaml
from pathlib import Path

def strip_container_prefix(path, workdir="/work"):
    """
    Convert container path back to host path if running in container.

    If path starts with /work and workdir exists, strip it.
    Otherwise return path unchanged.

    Args:
        path: Absolute or relative path (may include container prefix)
        workdir: Container working directory prefix to strip

    Returns:
        Path usable from host system
    """
    path_str = str(path)

    # If running inside container and path starts with /work
    if path_str.startswith(workdir + "/"):
        # Get the part after /work/
        stripped = path_str[len(workdir)+1:]

        # Convert to absolute path relative to current directory
        # (which is /work inside container, but actual dir on host)
        return os.path.abspath(stripped)

    return path_str

def create_symlink(source, target):
    """
    Create symlink with container path handling.

    Args:
        source: Source file (may have container prefix)
        target: Where to create symlink
    """
    # Strip container prefix from source path
    source_clean = strip_container_prefix(source)

    # Create parent directory if needed
    os.makedirs(os.path.dirname(target), exist_ok=True)

    # Create symlink
    if os.path.lexists(target):
        os.remove(target)

    os.symlink(source_clean, target)

def main():
    parser = argparse.ArgumentParser(
        description="Organize COSIGT pipeline resources and create config"
    )
    parser.add_argument("-g", "--graphs", required=True, 
                       help="TSV file mapping chromosomes to graph files")
    parser.add_argument("-r", "--reads", required=True,
                       help="TSV file mapping samples to alignment files")
    parser.add_argument("-b", "--bed", required=True,
                       help="BED file with regions of interest")
    parser.add_argument("-f", "--fasta", required=True,
                       help="Reference genome FASTA file")
    parser.add_argument("-o", "--outdir", required=True,
                       help="Output directory for pipeline results")
    parser.add_argument("-p", "--prefix", default="grch38#1#",
                       help="PanSN prefix for path names (default: grch38#1#)")

    args = parser.parse_args()

    # Create directory structure
    resources_dir = Path("resources")
    config_dir = Path("config")

    print("=" * 60)
    print("Reading inputs...")
    print("=" * 60)

    # Read graph mapping
    graphs = {}
    with open(args.graphs) as f:
        for line in f:
            if line.strip():
                chrom, graph_path = line.strip().split("\t")
                # Verify file exists
                clean_path = strip_container_prefix(graph_path)
                if not os.path.exists(clean_path):
                    print(f"WARNING: Graph file does not exist: {clean_path}")
                graphs[chrom] = graph_path

    print(f"Loaded table with graphs: {args.graphs}!")

    # Read BED file to get regions
    regions = {}
    with open(args.bed) as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                parts = line.strip().split("\t")
                chrom = parts[0]
                start = parts[1]
                end = parts[2]
                region_id = f"{chrom}_{start}_{end}"

                if chrom not in regions:
                    regions[chrom] = []
                regions[chrom].append({
                    "id": region_id,
                    "line": line.strip()
                })

    print(f"Loaded bed file: {args.bed}!")

    # Read alignment mapping
    samples = {}
    with open(args.reads) as f:
        for line in f:
            if line.strip():
                sample, aln_path = line.strip().split("\t")
                # Verify file exists
                clean_path = strip_container_prefix(aln_path)
                if not os.path.exists(clean_path):
                    print(f"WARNING: Alignment file does not exist: {clean_path}")
                samples[sample] = aln_path

    print(f"Loaded alignments from: {args.reads}!")

    # Verify reference exists
    clean_fasta = strip_container_prefix(args.fasta)
    if not os.path.exists(clean_fasta):
        print(f"WARNING: Reference file does not exist: {clean_fasta}")

    print(f"Loaded reference: {args.fasta}!")

    print("=" * 60)
    print("Writing resources...")
    print("=" * 60)

    # Create assemblies symlinks (organized by chromosome)
    assemblies_dir = resources_dir / "assemblies"
    for chrom, graph_path in graphs.items():
        chrom_dir = assemblies_dir / chrom
        chrom_dir.mkdir(parents=True, exist_ok=True)

        graph_name = os.path.basename(graph_path)
        target = chrom_dir / graph_name
        create_symlink(graph_path, target)

    print(f"Added graphs to: {assemblies_dir}/!")

    # Create alignments symlinks
    alignments_dir = resources_dir / "alignments"
    alignments_dir.mkdir(parents=True, exist_ok=True)

    for sample, aln_path in samples.items():
        aln_name = os.path.basename(aln_path)
        target = alignments_dir / aln_name
        create_symlink(aln_path, target)

        # Also symlink index files if they exist
        for ext in [".bai", ".crai", ".csi"]:
            idx_path = aln_path + ext
            clean_idx = strip_container_prefix(idx_path)
            if os.path.exists(clean_idx):
                idx_target = alignments_dir / (aln_name + ext)
                create_symlink(idx_path, idx_target)

    print(f"Added alignments to: {alignments_dir}/!")

    # Create regions BED files (organized by chromosome)
    regions_dir = resources_dir / "regions"
    for chrom, chrom_regions in regions.items():
        chrom_dir = regions_dir / chrom
        chrom_dir.mkdir(parents=True, exist_ok=True)

        for region in chrom_regions:
            bed_file = chrom_dir / f"{region['id']}.bed"
            with open(bed_file, 'w') as f:
                f.write(region['line'] + "\n")

    print(f"Added regions to: {regions_dir}/!")

    # Create reference symlinks
    reference_dir = resources_dir / "reference"
    reference_dir.mkdir(parents=True, exist_ok=True)

    ref_name = os.path.basename(args.fasta)
    ref_target = reference_dir / ref_name
    create_symlink(args.fasta, ref_target)

    # Also symlink FAI index if exists
    fai_path = args.fasta + ".fai"
    clean_fai = strip_container_prefix(fai_path)
    if os.path.exists(clean_fai):
        fai_target = reference_dir / (ref_name + ".fai")
        create_symlink(fai_path, fai_target)

    print(f"Added reference to: {reference_dir}/!")

    print("=" * 60)
    print("Writing configuration...")
    print("=" * 60)

    # Create config.yaml
    config_dir.mkdir(parents=True, exist_ok=True)

    # Strip /work prefix from paths in config
    config = {
        "output": strip_container_prefix(args.outdir),
        "pansn_prefix": args.prefix,
        "reference": strip_container_prefix(str(ref_target)),
        "regions": sorted([r["id"] for regions_list in regions.values() for r in regions_list]),
        "samples": sorted(samples.keys())
    }

    config_file = config_dir / "config.yaml"
    with open(config_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False)

    print(f"Wrote config to: {config_file}!")

    print("=" * 60)
    print("✓ Setup complete!")
    print("=" * 60)

if __name__ == "__main__":
    main()

