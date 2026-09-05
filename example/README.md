# Pre-computed Example Reference Model Bundle

## Overview

This directory contains the baseline two-layer reference model derived from **3,494 high-coverage genomes** (1000G + HGDP) across **603,906 quality-filtered SNPs** on **GRCh38**.

The bundle is designed for deterministic two-layer ancestry coordinate assignment in clinical and translational workflows, where each sample is projected into a global ancestry space first and then refined in a superpopulation-specific sublayer.

## Directory Layout & Asset Footprint

Directory structure:

- `models/layer1/global/`
- `models/layer2/<SUPERPOP>/`
- `provenance.json`

Asset categories:

- Tracked in Git (Lightweight): `*.eigenval`, `*.eigenvec`, `*_samples.tsv`, `*.json`
- Distributed via Release Tarball (`nf_PopPCA_refgen_v1.0.0_models.tar.gz`): `*.normalized_weights.tsv`, `*.eigenvec.allele`, `*.afreq`

## Model Unpacking & Staging Instructions

Heavy matrices are distributed via `nf_PopPCA_refgen_v1.0.0_models.tar.gz` on GitHub Releases.

Run all commands below from the repository root: `repository root`.

Option A (Automated):

```bash
./bin/download_models.sh --outdir example/models
```

Option B (Manual):

```bash
tar -xzvf nf_PopPCA_refgen_v1.0.0_models.tar.gz -C example/
```

Verification:

```bash
ls -lh example/models/layer1/global/*.normalized_weights.tsv example/models/layer2/*/*.normalized_weights.tsv
```

## Downstream Usage Example (Two-Layer Projection Workflow)

### Workflow Diagram

```mermaid
graph TD
    A[Mapped BAM/CRAM or Unphased VCF] --> B[Targeted Genotyping (603k Loci via bcftools)]
    B --> C[Layer 1 Global Projection (PC1-PC10)]
    C --> D[Superpopulation Classifier]
    D --> E[Layer 2 Sublayer Routing (AFR/AMR/EAS/EUR/SAS)]
    E --> F[Downstream Applications (Phasing Reference Selection / ACMG Variant Filtering)]
```

### Targeted Input Options & Fast-Path Mode

Inputs can be full `VCF/BCF` files or targeted `BAM/CRAM` pileups restricted to the 603,906 loci. For fast pre-call ancestry determination, run `bcftools mpileup` against the panel BED file `example/models/layer1/global/603k_loci.bed`.

```bash
bcftools mpileup \
  -Ou \
  -f GRCh38.fa \
  -R example/models/layer1/global/603k_loci.bed \
  sample.bam \
| bcftools call -mv -Ob -o patient_603k_targets.bcf
```

This targeted path typically finishes in about 5 minutes on a 30x WGS sample before whole-genome calling or phasing.

### Layer 1 Projection Node

1. Harmonize targeted calls and build a normalized genotype matrix `X_{ij}` in the same variant order used by the global model.
2. Compute global coordinates with scoring `PC = X · W`.
3. Assign continental superpopulation class: `AFR`, `AMR`, `EAS`, `EUR`, or `SAS`.

Exact PLINK 2 scoring command:

```bash
plink2 --bcf patient_603k_targets.bcf --score example/models/layer1/global/global_pca.normalized_weights.tsv 1 2 3 header --out patient_global_pcs
```

### Layer 2 Routing

Route the sample to `example/models/layer2/<SUPERPOP>/` using the Layer 1 class label, then project with the corresponding sublayer matrix for fine-grained continental refinement.

### Upstream Value

- Dynamic phasing reference selection for Eagle2/SHAPEIT4 based on inferred superpopulation.
- Population-aware ACMG variant frequency selection using matched gnomAD populations.

### Implementation Snippets

Python `TwoLayerProjector` class:

```python
from pathlib import Path


class TwoLayerProjector:
    """Run two-layer ancestry projection from targeted genotype evidence."""

    def __init__(self, model_root: Path):
        """Set layer-1 global and layer-2 superpopulation model roots."""
        self.model_root = model_root
        self.global_dir = model_root / "layer1" / "global"
        self.layer2_root = model_root / "layer2"

    def targeted_pileup(self, sample_bam: Path, reference_fasta: Path) -> Path:
        """Generate a targeted BCF at the 603k loci from BAM/CRAM alignments."""
        return run_targeted_mpileup(
            sample_bam=sample_bam,
            reference_fasta=reference_fasta,
            loci_bed=self.global_dir / "603k_loci.bed",
            output_bcf=Path("patient_603k_targets.bcf"),
        )

    def project_global_pca(self, patient_bcf: Path):
        """Project a sample into Layer 1 global PCA coordinates (PC1-PC10)."""
        return project_with_weights(
            bcf=patient_bcf,
            allele_model=self.global_dir / "global_pca.eigenvec.allele",
            weights=self.global_dir / "global_pca.normalized_weights.tsv",
            samples=self.global_dir / "global_samples.tsv",
        )

    def classify_superpop(self, global_scores) -> str:
        """Assign AFR/AMR/EAS/EUR/SAS from Layer 1 projected coordinates."""
        return classify_superpopulation(global_scores)

    def project_sublayer_pca(self, patient_bcf: Path, superpop: str):
        """Project into superpopulation-specific Layer 2 PCA coordinates."""
        model_dir = self.layer2_root / superpop
        return project_with_weights(
            bcf=patient_bcf,
            allele_model=model_dir / f"{superpop}_pca.eigenvec.allele",
            weights=model_dir / f"{superpop}_pca.normalized_weights.tsv",
            samples=model_dir / f"{superpop}_samples.tsv",
        )

    def run_two_layer_projection(self, sample_bam: Path, reference_fasta: Path):
        """Execute targeted pileup, global projection, superpop assignment, and sublayer projection."""
        patient_bcf = self.targeted_pileup(sample_bam, reference_fasta)
        global_scores = self.project_global_pca(patient_bcf)
        superpop = self.classify_superpop(global_scores)
        sublayer_scores = self.project_sublayer_pca(patient_bcf, superpop)
        return global_scores, superpop, sublayer_scores
```

Nextflow DSL2 process outline:

```nextflow
nextflow.enable.dsl = 2

process TARGETED_PILEUP {
    input:
    path sample_bam
    path reference_fasta

    output:
    path "patient_603k_targets.bcf"

    script:
    """
    bcftools mpileup -Ou -f ${reference_fasta} -R example/models/layer1/global/603k_loci.bed ${sample_bam} \
      | bcftools call -mv -Ob -o patient_603k_targets.bcf
    """
}

process PROJECT_GLOBAL_PCA {
    input:
    path patient_bcf

    output:
    path "patient_global_pcs.sscore"

    script:
    """
    plink2 --bcf ${patient_bcf} --score example/models/layer1/global/global_pca.normalized_weights.tsv 1 2 3 header --out patient_global_pcs
    """
}

process CLASSIFY_SUPERPOP {
    input:
    path layer1_scores

    output:
    val superpop_label

    script:
    """
    python classify_superpop.py --scores ${layer1_scores} > superpop_label.txt
    cat superpop_label.txt
    """
}

process PROJECT_SUBLAYER_PCA {
    input:
    path patient_bcf
    val superpop_label

    output:
    path "patient_sublayer_pcs.tsv"

    script:
    """
    model_dir=example/models/layer2/${superpop_label}
    python project_sublayer.py \
      --bcf ${patient_bcf} \
      --weights ${model_dir}/${superpop_label}_pca.normalized_weights.tsv \
      --alleles ${model_dir}/${superpop_label}_pca.eigenvec.allele \
      --samples ${model_dir}/${superpop_label}_samples.tsv \
      --out patient_sublayer_pcs.tsv
    """
}

workflow TWO_LAYER_PROJECTION {
    take:
    ch_bam
    ch_ref

    main:
    TARGETED_PILEUP(ch_bam, ch_ref)
    PROJECT_GLOBAL_PCA(TARGETED_PILEUP.out)
    CLASSIFY_SUPERPOP(PROJECT_GLOBAL_PCA.out)
    PROJECT_SUBLAYER_PCA(TARGETED_PILEUP.out, CLASSIFY_SUPERPOP.out)

    emit:
    PROJECT_GLOBAL_PCA.out
    CLASSIFY_SUPERPOP.out
    PROJECT_SUBLAYER_PCA.out
}
```

Practical note: keep harmonization, normalization, and variant-ordering identical between reference training and downstream projection to prevent coordinate drift.
