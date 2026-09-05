process GENERATE_GLOBAL_PCA {
    tag 'global-pca'
    publishDir "${params.outdir}/models/layer1/global", mode: 'copy', overwrite: true

    input:
    tuple path(pgen_files), path(pvar_files), path(psam_files), path(sample_tables), path(qc_summaries)

    output:
    path 'global_pca.eigenvec', emit: eigenvec
    path 'global_pca.eigenval', emit: eigenval
    path 'global_pca.eigenvec.allele', emit: allele_weights
    path 'global_pca.normalized_weights.tsv', emit: normalized_weights
    path 'global_pca.afreq', emit: allele_freqs
    path 'global_samples.tsv', emit: sample_manifest
    path 'global_summary.json', emit: summary

    stub:
    """
    printf '2.0\n1.8\n1.6\n1.4\n1.2\n1.0\n0.9\n0.8\n0.7\n0.6\n' > global_pca.eigenval
    printf '#FID\tIID\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\nGLOB\tS1\t0.1\t0.1\t0.1\t0.1\t0.1\t0.1\t0.1\t0.1\t0.1\t0.1\n' > global_pca.eigenvec
    printf 'CHROM\tPOS\tID\tREF\tALT\tA1\tAX\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tPC7\tPC8\tPC9\tPC10\nchr1\t1000\tchr1:1000:A:G\tA\tG\tG\tA\t0.1\t0.1\t0.1\t0.1\t0.1\t0.1\t0.1\t0.1\t0.1\t0.1\n' > global_pca.eigenvec.allele
    cp global_pca.eigenvec.allele global_pca.normalized_weights.tsv
    printf 'ID\tALT_FREQS\nchr1:1000:A:G\t0.4\n' > global_pca.afreq
    printf 'FID\tIID\tSUPERPOP\tDATASET\tASSAY\tBUILD\n' > global_samples.tsv
    for f in *.samples.tsv; do awk 'FNR>1' "\${f}" >> global_samples.tsv; done
    sha256sum global_pca.eigenvec global_pca.eigenval global_pca.eigenvec.allele global_pca.normalized_weights.tsv global_pca.afreq global_samples.tsv > checksums.txt
    eigenvec_sha=\$(sha256sum global_pca.eigenvec | cut -d ' ' -f 1)
    eigenval_sha=\$(sha256sum global_pca.eigenval | cut -d ' ' -f 1)
    allele_sha=\$(sha256sum global_pca.eigenvec.allele | cut -d ' ' -f 1)
    weights_sha=\$(sha256sum global_pca.normalized_weights.tsv | cut -d ' ' -f 1)
    afreq_sha=\$(sha256sum global_pca.afreq | cut -d ' ' -f 1)
    samples_sha=\$(sha256sum global_samples.tsv | cut -d ' ' -f 1)

    cat <<EOF > global_summary.json
{
  "layer": "global",
  "sample_count": 1,
  "variant_count": ${params.min_shared_variants},
  "build": "${params.target_build}",
  "files": {
    "layer1/global/global_pca.eigenvec": "\${eigenvec_sha}",
    "layer1/global/global_pca.eigenval": "\${eigenval_sha}",
    "layer1/global/global_pca.eigenvec.allele": "\${allele_sha}",
    "layer1/global/global_pca.normalized_weights.tsv": "\${weights_sha}",
    "layer1/global/global_pca.afreq": "\${afreq_sha}",
    "layer1/global/global_samples.tsv": "\${samples_sha}"
  }
}
EOF
    """

    script:
    """
    prefixes=()
    for pgen in *.qc.pgen; do
      prefixes+=("\${pgen%.pgen}")
    done

    [[ "\${#prefixes[@]}" -gt 0 ]] || { echo 'No QCed panels were staged for global PCA.' >&2; exit 1; }

    if [[ "\${#prefixes[@]}" -eq 1 ]]; then
      cp "\${prefixes[0]}.pgen" global_merged.pgen
      cp "\${prefixes[0]}.pvar" global_merged.pvar
      cp "\${prefixes[0]}.psam" global_merged.psam
    else
      # Build a strict common-variant set across all inputs before merge.
      _first=1
      for _pfx in "\${prefixes[@]}"; do
        awk '!/^#/ { print \$3 }' "\${_pfx}.pvar" | sort > "_pqc_vars_\${_pfx}.txt"
        if [[ "\${_first}" -eq 1 ]]; then
          cp "_pqc_vars_\${_pfx}.txt" _common_vars.txt
          _first=0
        else
          comm -12 _common_vars.txt "_pqc_vars_\${_pfx}.txt" > _common_new.txt
          mv _common_new.txt _common_vars.txt
        fi
      done

      [[ -s _common_vars.txt ]] || { echo 'Common variant set is empty after QC/pruning.' >&2; exit 1; }

      # Convert each common-variant panel to a PLINK1 binary set so we can
      # merge with --bmerge-list (plink2 --pmerge-list is not fully implemented
      # for this non-concatenating PGEN merge pattern in a5.12).
      for _pfx in "\${prefixes[@]}"; do
        plink2 \
          --pfile "\${_pfx}" \
          --extract _common_vars.txt \
          --make-bed \
          --out "\${_pfx}.bedset"
      done

      printf '%s.bedset\n' "\${prefixes[@]:1}" > _bmerge_list.txt
      plink \
        --bfile "\${prefixes[0]}.bedset" \
        --merge-list _bmerge_list.txt \
        --make-bed \
        --out global_merged_bed

      plink2 \
        --bfile global_merged_bed \
        --make-pgen \
        --out global_merged
    fi

    plink2 --pfile global_merged --freq --out global_pca
    plink2 --pfile global_merged --pca ${params.pca_components} allele-wts --out global_pca

    awk 'NR==FNR {eig[NR]=\$1; next} NR==1 {printf "%s", \$0; for (i=1; i<=length(eig); i++) printf "\\tNORM_PC%d", i; printf "\\n"; next} {printf "%s", \$0; for (i=1; i<=length(eig); i++) {col=(NF-length(eig))+i; printf "\\t%.10f", (\$(col)+0)/sqrt(eig[i])} printf "\\n"}' global_pca.eigenval global_pca.eigenvec.allele > global_pca.normalized_weights.tsv

    printf 'FID\tIID\tSUPERPOP\tDATASET\tASSAY\tBUILD\n' > global_samples.tsv
    for f in *.samples.tsv; do
      awk 'FNR>1' "\${f}" >> global_samples.tsv
    done

    sample_count=\$(( \$(wc -l < global_samples.tsv) - 1 ))
    variant_count=\$(( \$(wc -l < global_merged.pvar) - 1 ))
    sha256sum global_pca.eigenvec global_pca.eigenval global_pca.eigenvec.allele global_pca.normalized_weights.tsv global_pca.afreq global_samples.tsv > checksums.txt
    eigenvec_sha=\$(sha256sum global_pca.eigenvec | cut -d ' ' -f 1)
    eigenval_sha=\$(sha256sum global_pca.eigenval | cut -d ' ' -f 1)
    allele_sha=\$(sha256sum global_pca.eigenvec.allele | cut -d ' ' -f 1)
    weights_sha=\$(sha256sum global_pca.normalized_weights.tsv | cut -d ' ' -f 1)
    afreq_sha=\$(sha256sum global_pca.afreq | cut -d ' ' -f 1)
    samples_sha=\$(sha256sum global_samples.tsv | cut -d ' ' -f 1)

    cat <<EOF > global_summary.json
{
  "layer": "global",
  "sample_count": \${sample_count},
  "variant_count": \${variant_count},
  "build": "${params.target_build}",
  "files": {
    "layer1/global/global_pca.eigenvec": "\${eigenvec_sha}",
    "layer1/global/global_pca.eigenval": "\${eigenval_sha}",
    "layer1/global/global_pca.eigenvec.allele": "\${allele_sha}",
    "layer1/global/global_pca.normalized_weights.tsv": "\${weights_sha}",
    "layer1/global/global_pca.afreq": "\${afreq_sha}",
    "layer1/global/global_samples.tsv": "\${samples_sha}"
  }
}
EOF
    """
}