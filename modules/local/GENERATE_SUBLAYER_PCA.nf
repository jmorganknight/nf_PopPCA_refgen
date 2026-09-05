process GENERATE_SUBLAYER_PCA {
  tag 'sub-pca'
  publishDir "${params.outdir}/models/layer2", mode: 'copy', overwrite: true, saveAs: { name ->
    name.startsWith('sub/') ? name.substring(4) : name
  }

    input:
  tuple path(pgen_files), path(pvar_files), path(psam_files), path(sample_tables), path(qc_summaries)

    output:
  path 'sub/*/*_pca.eigenvec', emit: eigenvec
  path 'sub/*/*_pca.eigenval', emit: eigenval
  path 'sub/*/*_pca.eigenvec.allele', emit: allele_weights
  path 'sub/*/*_pca.normalized_weights.tsv', emit: normalized_weights
  path 'sub/*/*_pca.afreq', emit: allele_freqs
  path 'sub/*/*_samples.tsv', emit: sample_manifest
  path 'sub/sub_summary.json', emit: summary

    stub:
    """
    mkdir -p sub/AFR
    cat <<'EOF' > sub/AFR/AFR_pca.eigenval
2.0
1.8
1.6
1.4
1.2
1.0
0.9
0.8
0.7
0.6
EOF
    cat <<'EOF' > sub/AFR/AFR_pca.eigenvec
#FID	IID	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10
AFR	S1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1
EOF
    cat <<'EOF' > sub/AFR/AFR_pca.eigenvec.allele
CHROM	POS	ID	REF	ALT	A1	AX	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10
chr1	1000	chr1:1000:A:G	A	G	G	A	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1	0.1
EOF
    cp sub/AFR/AFR_pca.eigenvec.allele sub/AFR/AFR_pca.normalized_weights.tsv
    cat <<'EOF' > sub/AFR/AFR_pca.afreq
ID	ALT_FREQS
chr1:1000:A:G	0.4
EOF
    cat <<'EOF' > sub/AFR/AFR_samples.tsv
FID	IID	SUPERPOP	DATASET	ASSAY	BUILD
AFR	S1	AFR	stub	WGS	GRCh38
EOF
    eigenvec_sha=\$(sha256sum sub/AFR/AFR_pca.eigenvec | cut -d ' ' -f 1)
    eigenval_sha=\$(sha256sum sub/AFR/AFR_pca.eigenval | cut -d ' ' -f 1)
    allele_sha=\$(sha256sum sub/AFR/AFR_pca.eigenvec.allele | cut -d ' ' -f 1)
    weights_sha=\$(sha256sum sub/AFR/AFR_pca.normalized_weights.tsv | cut -d ' ' -f 1)
    afreq_sha=\$(sha256sum sub/AFR/AFR_pca.afreq | cut -d ' ' -f 1)
    samples_sha=\$(sha256sum sub/AFR/AFR_samples.tsv | cut -d ' ' -f 1)
    cat <<EOF > sub/sub_summary.json
{
  "layer": "sub",
  "build": "${params.target_build}",
  "models": [
    {
      "layer": "sub",
      "superpop": "AFR",
      "sample_count": 1,
      "variant_count": ${params.min_shared_variants},
      "files": {
        "layer2/AFR/AFR_pca.eigenvec": "\${eigenvec_sha}",
        "layer2/AFR/AFR_pca.eigenval": "\${eigenval_sha}",
        "layer2/AFR/AFR_pca.eigenvec.allele": "\${allele_sha}",
        "layer2/AFR/AFR_pca.normalized_weights.tsv": "\${weights_sha}",
        "layer2/AFR/AFR_pca.afreq": "\${afreq_sha}",
        "layer2/AFR/AFR_samples.tsv": "\${samples_sha}"
      }
    }
  ]
}
EOF
    """

    script:
    """
    mkdir -p sub

    prefixes=()
    for pgen in *.qc.pgen; do
      prefixes+=("\${pgen%.pgen}")
    done

    [[ "\${#prefixes[@]}" -gt 0 ]] || { echo 'No QCed panels were staged for sublayer PCA.' >&2; exit 1; }

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

    printf 'IID\tFID\tSUPERPOP\tDATASET\tASSAY\tBUILD\n' > global_samples.tsv
    for f in *.samples.tsv; do
      awk 'BEGIN{FS="[[:space:]]+"; OFS="\t"} FNR>1 { gsub(/\r/,""); if (\$1!="") print \$1,\$2,\$3,\$4,\$5,\$6 }' "\${f}" >> global_samples.tsv
    done

    awk 'NR>1 && \$3!="" && \$3!="GLOBAL" && \$3!="SUPERPOP" { print \$3 }' global_samples.tsv | sort -u > _superpops.txt
    [[ -s _superpops.txt ]] || { echo 'No valid superpopulation labels found.' >&2; exit 1; }

    cat <<EOF > sub/sub_summary.json
{
  "layer": "sub",
  "build": "${params.target_build}",
  "models": [
EOF

    model_count=0
    while IFS= read -r pop; do
      [[ -n "\${pop}" ]] || continue
      mkdir -p "sub/\${pop}"

      awk -v p="\${pop}" 'NR>1 && \$3==p { print \$1 }' global_samples.tsv > "sub/\${pop}/\${pop}_keep.txt"
      if [[ ! -s "sub/\${pop}/\${pop}_keep.txt" ]]; then
        echo "No samples found for superpopulation \${pop}; skipping." >&2
        continue
      fi

      plink2 --pfile global_merged --keep "sub/\${pop}/\${pop}_keep.txt" --make-pgen --out "sub_\${pop}_merged"

      sample_count=\$(( \$(wc -l < sub_\${pop}_merged.psam) - 1 ))
      if [[ "\${sample_count}" -lt 30 ]]; then
        echo "Skipping \${pop}: sample_count=\${sample_count} < 30." >&2
        continue
      fi

      plink2 --pfile "sub_\${pop}_merged" --freq --out "sub/\${pop}/\${pop}_pca"
      plink2 --pfile "sub_\${pop}_merged" --pca ${params.pca_components} allele-wts --out "sub/\${pop}/\${pop}_pca"

      awk 'NR==FNR {eig[NR]=\$1; next} NR==1 {printf "%s", \$0; for (i=1; i<=length(eig); i++) printf "\\tNORM_PC%d", i; printf "\\n"; next} {printf "%s", \$0; for (i=1; i<=length(eig); i++) {col=(NF-length(eig))+i; printf "\\t%.10f", (\$(col)+0)/sqrt(eig[i])} printf "\\n"}' "sub/\${pop}/\${pop}_pca.eigenval" "sub/\${pop}/\${pop}_pca.eigenvec.allele" > "sub/\${pop}/\${pop}_pca.normalized_weights.tsv"

      awk -v p="\${pop}" 'NR==1 { print; next } \$3==p { print }' global_samples.tsv > "sub/\${pop}/\${pop}_samples.tsv"
      variant_count=\$(( \$(wc -l < sub_\${pop}_merged.pvar) - 1 ))

      eigenvec_sha=\$(sha256sum "sub/\${pop}/\${pop}_pca.eigenvec" | cut -d ' ' -f 1)
      eigenval_sha=\$(sha256sum "sub/\${pop}/\${pop}_pca.eigenval" | cut -d ' ' -f 1)
      allele_sha=\$(sha256sum "sub/\${pop}/\${pop}_pca.eigenvec.allele" | cut -d ' ' -f 1)
      weights_sha=\$(sha256sum "sub/\${pop}/\${pop}_pca.normalized_weights.tsv" | cut -d ' ' -f 1)
      afreq_sha=\$(sha256sum "sub/\${pop}/\${pop}_pca.afreq" | cut -d ' ' -f 1)
      samples_sha=\$(sha256sum "sub/\${pop}/\${pop}_samples.tsv" | cut -d ' ' -f 1)

      if [[ "\${model_count}" -gt 0 ]]; then
        printf ',\n' >> sub/sub_summary.json
      fi

      cat <<EOF >> sub/sub_summary.json
    {
      "layer": "sub",
      "superpop": "\${pop}",
      "sample_count": \${sample_count},
      "variant_count": \${variant_count},
      "files": {
        "layer2/\${pop}/\${pop}_pca.eigenvec": "\${eigenvec_sha}",
        "layer2/\${pop}/\${pop}_pca.eigenval": "\${eigenval_sha}",
        "layer2/\${pop}/\${pop}_pca.eigenvec.allele": "\${allele_sha}",
        "layer2/\${pop}/\${pop}_pca.normalized_weights.tsv": "\${weights_sha}",
        "layer2/\${pop}/\${pop}_pca.afreq": "\${afreq_sha}",
        "layer2/\${pop}/\${pop}_samples.tsv": "\${samples_sha}"
      }
    }
EOF

      model_count=\$((model_count + 1))
    done < _superpops.txt

    [[ "\${model_count}" -gt 0 ]] || { echo 'No sublayer models were generated after filtering.' >&2; exit 1; }

    cat <<EOF >> sub/sub_summary.json

  ]
}
EOF
    """
}