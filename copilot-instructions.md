# GitHub Copilot Rules & Constraints for `nf_PopPCA_refgen`

## Directory & Context Rules
1. **Target Directory Constraint:** ALL generated files, scripts, configs, and modules MUST be created strictly inside the root directory `nf_PopPCA_refgen/`. Never create nested wrapper directories like `nf_PopPCA_refgen/nf_PopPCA_refgen/`.
2. **DSL2 Syntax Standard:** Write modern Nextflow DSL2 syntax (`nextflow.enable.dsl=2`). Do NOT use legacy DSL1 code constructs.
3. **Module Isolation:** Every process MUST live in its own file under `modules/local/<PROCESS_NAME>.nf` and be explicitly imported into `main.nf`. Do NOT write monolithic single-file scripts.

## Code Style & Formatting Guidelines
1. **Naming Conventions:**
   - Processes: UPPERCASE with underscores (e.g., `NORMALIZE_INPUTS`).
   - Channels: Prefix with `ch_` in lowercase (e.g., `ch_active_datasets`).
   - Parameters: Lowercase with underscores (e.g., `params.ref_manifest`).
2. **Publish Directory:** Always specify `publishDir "${params.outdir}/models/...", mode: 'copy', overwrite: true` inside processes emitting final reference assets.
3. **Container Definition:** Define the default container (`wes-onco-core:1.0.0` or `biocontainers/plink2:2.00a5.10--h9a82719_0`) in `nextflow.config`—do not hardcode container strings repeatedly across modules.

## Biological & Technical Guardrails
1. **Microarray Palindrome Filter:** Always include `plink2 --snps-only just-acgt` and explicit filtering for A/T, T/A, C/G, and G/C palindromic SNPs in `NORMALIZE_INPUTS.nf`.
2. **Intersection Guardrail:** In `INTERSECT_VARIANTS.nf`, assert that the shared variant count is >= 100,000 before proceeding to PCA.
3. **Sample Threshold:** Do not attempt SVD/PCA on Layer 2 sub-population groups with fewer than N = 30 samples.