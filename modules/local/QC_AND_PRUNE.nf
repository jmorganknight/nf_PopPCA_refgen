process QC_AND_PRUNE {
    tag { "qc:${meta.id}" }

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam), path(variant_table), path(sample_table), path(normalized_summary)
    path shared_variants
    path high_ld_regions

    output:
    tuple val(meta),
        path("${meta.id}.qc.pgen"),
        path("${meta.id}.qc.pvar"),
        path("${meta.id}.qc.psam"),
        path("${meta.id}.qc.samples.tsv"),
        path("${meta.id}.qc.summary.json"),
        emit: qc_panels

    stub:
    def sampleCount = Math.max((meta.sample_count ?: 4) as Integer, 4)
    def variantCount = Math.max(params.min_shared_variants as Integer, 10)
    """
    cp ${pvar} ${meta.id}.qc.pvar
    cp ${psam} ${meta.id}.qc.psam
    cp ${sample_table} ${meta.id}.qc.samples.tsv
    : > ${meta.id}.qc.pgen
    cat <<EOF > ${meta.id}.qc.summary.json
{
  "dataset_id": "${meta.id}",
  "superpop": "${meta.superpop}",
  "sample_count": ${sampleCount},
  "variant_count": ${variantCount},
  "maf_filter": ${params.maf},
  "hwe_filter": ${params.hwe},
  "geno_filter": ${params.geno},
  "king_cutoff": ${params.king_cutoff},
  "high_ld_mask": "${high_ld_regions.getName()}",
  "ld_pruning": "${params.ld_window} ${params.ld_step} ${params.ld_r2}"
}
EOF
    """

    script:
    """
    input_prefix="${pgen.baseName}"

    plink2 \
      --pfile "\${input_prefix}" \
      --extract ${shared_variants} \
      --maf ${params.maf} \
      --hwe ${params.hwe} midp \
      --geno ${params.geno} \
      --snps-only just-acgt \
      --make-pgen \
      --out ${meta.id}.qc_step1

    plink2 \
      --pfile ${meta.id}.qc_step1 \
      --make-king triangle bin \
      --out ${meta.id}.king

    plink2 \
      --pfile ${meta.id}.qc_step1 \
      --king-cutoff ${meta.id}.king ${params.king_cutoff} \
      --make-pgen \
      --out ${meta.id}.unrelated

    plink2 \
      --pfile ${meta.id}.unrelated \
      --exclude range ${high_ld_regions} \
      --indep-pairwise ${params.ld_window} ${params.ld_step} ${params.ld_r2} \
      --out ${meta.id}.pruned

    plink2 \
      --pfile ${meta.id}.unrelated \
      --extract ${meta.id}.pruned.prune.in \
      --make-pgen \
      --out ${meta.id}.qc

    awk 'NR==FNR && FNR>1 { keep[\$1 FS \$2]=1; next } FNR==1 { print; next } ((\$1 FS \$2) in keep)' ${meta.id}.qc.psam ${sample_table} > ${meta.id}.qc.samples.tsv

    sample_count=\$(( \$(wc -l < ${meta.id}.qc.psam) - 1 ))
    variant_count=\$(( \$(wc -l < ${meta.id}.qc.pvar) - 1 ))

    cat <<EOF > ${meta.id}.qc.summary.json
{
  "dataset_id": "${meta.id}",
  "superpop": "${meta.superpop}",
  "sample_count": \${sample_count},
  "variant_count": \${variant_count},
  "maf_filter": ${params.maf},
  "hwe_filter": ${params.hwe},
  "geno_filter": ${params.geno},
  "king_cutoff": ${params.king_cutoff},
  "high_ld_mask": "${high_ld_regions.getName()}",
  "ld_pruning": "${params.ld_window} ${params.ld_step} ${params.ld_r2}"
}
EOF
    """
}