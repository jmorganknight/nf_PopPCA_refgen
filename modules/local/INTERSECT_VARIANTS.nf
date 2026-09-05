process INTERSECT_VARIANTS {
    tag 'variant-intersection'

    input:
    path variant_tables

    output:
    path 'shared_variants.txt', emit: shared_variants
    path 'shared_variants.summary.json', emit: summary

    stub:
    """
    seq 1 ${params.min_shared_variants} | awk '{ print "chr1:" 999+\$1 ":A:G" }' > shared_variants.txt
    cat <<EOF > shared_variants.summary.json
{
  "layer": "intersection",
  "shared_variant_count": ${params.min_shared_variants},
  "minimum_required": ${params.min_shared_variants}
}
EOF
    """

    script:
    """
    dataset_count=\$(find . -maxdepth 1 -name '*.variants.tsv' | wc -l)
    [[ "\${dataset_count}" -gt 0 ]] || { echo 'No variant tables were staged for intersection.' >&2; exit 1; }

    tmp_dir=intersection_tmp
    mkdir -p "\${tmp_dir}"
    for f in *.variants.tsv; do
      sort -u "\${f}" > "\${tmp_dir}/\${f}.sorted"
    done

    cat "\${tmp_dir}"/*.sorted | sort | uniq -c | awk -v n="\${dataset_count}" '\$1 == n { print \$2 }' > shared_variants.txt

    shared_count=\$(wc -l < shared_variants.txt)
    if [[ "\${shared_count}" -lt ${params.min_shared_variants} ]]; then
      echo "Shared variant core contains only \${shared_count} SNPs; expected at least ${params.min_shared_variants}." >&2
      exit 1
    fi

    cat <<EOF > shared_variants.summary.json
{
  "layer": "intersection",
  "shared_variant_count": \${shared_count},
  "minimum_required": ${params.min_shared_variants}
}
EOF
    """
}