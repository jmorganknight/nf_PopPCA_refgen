process NORMALIZE_INPUTS {
    tag { "normalize:${meta.id}" }

    input:
    tuple val(meta), path(input_files)

    output:
    tuple val(meta),
        path("${meta.id}.norm.pgen"),
        path("${meta.id}.norm.pvar"),
        path("${meta.id}.norm.psam"),
        path("${meta.id}.norm.variants.tsv"),
        path("${meta.id}.norm.samples.tsv"),
        path("${meta.id}.norm.summary.json"),
        emit: normalized_panels

    stub:
    def sampleCount = Math.max((meta.sample_count ?: 4) as Integer, 4)
    """
    cat <<'EOF' > ${meta.id}.norm.pvar
#CHROM  POS ID  REF ALT
chr1    1000    chr1:1000:A:G   A   G
chr1    1010    chr1:1010:C:T   C   T
chr1    1020    chr1:1020:G:A   G   A
chr1    1030    chr1:1030:T:C   T   C
chr1    1040    chr1:1040:G:T   G   T
chr1    1050    chr1:1050:C:A   C   A
chr1    1060    chr1:1060:A:C   A   C
chr1    1070    chr1:1070:G:C   G   C
chr1    1080    chr1:1080:T:G   T   G
chr1    1090    chr1:1090:A:G   A   G
EOF
    : > ${meta.id}.norm.pgen
    {
      printf '#FID\tIID\tSEX\tPHENO\n'
      for i in \$(seq 1 ${sampleCount}); do
        printf '%s\t%s_%02d\t0\t-9\n' ${meta.id} ${meta.id} "\${i}"
      done
    } > ${meta.id}.norm.psam
    awk 'NR>1 { print \$3 }' ${meta.id}.norm.pvar > ${meta.id}.norm.variants.tsv
    awk 'BEGIN{OFS="\t"; print "IID","FID","SUPERPOP","DATASET","ASSAY","BUILD"} NR>1 { print \$1,"NA","${meta.superpop}","${meta.id}","${meta.assay}","${meta.build}" }' ${meta.id}.norm.psam > ${meta.id}.norm.samples.tsv
    cat <<EOF > ${meta.id}.norm.summary.json
{
  "dataset_id": "${meta.id}",
  "superpop": "${meta.superpop}",
  "assay": "${meta.assay}",
  "build": "${meta.build}",
  "normalized_sample_count": ${sampleCount},
  "normalized_variant_count": 10,
  "palindromic_filter_applied": true
}
EOF
    """

    script:
    def refFromFaArg = params.ref_fasta ? "--ref-from-fa ${params.ref_fasta}" : ''
    def liftoverRef = params.liftover_ref_fasta ?: params.ref_fasta
    def liftoverRefArg = liftoverRef ? "${liftoverRef}" : ''
    def liftoverChainArg = params.liftover_chain ? "${params.liftover_chain}" : ''
    def liftoverMapDirArg = params.liftover_map_dir ? "${params.liftover_map_dir}" : ''
    """
    dataset_build="${meta.build}"
    target_build="${params.target_build}"
    source_prefix=""

    case "${meta.input_format}" in
      vcf)
        shopt -s nullglob
        vcf_candidates=( *.vcf *.vcf.gz )
        [[ \${#vcf_candidates[@]} -gt 0 ]] || { echo "VCF input not found for ${meta.id}" >&2; exit 1; }

        if [[ \${#vcf_candidates[@]} -eq 1 ]]; then
          plink2 \
            --vcf "\${vcf_candidates[0]}" \
            --max-alleles 2 \
            --snps-only just-acgt \
            --rm-dup force-first \
            --make-pgen \
            --out ${meta.id}.raw
        else
          # plink2 --pmerge-list does not support vcf mode in a5.12+;
          # convert each per-chromosome VCF to pgen independently, then merge
          _pidx=0
          for _vcf in "\${vcf_candidates[@]}"; do
            plink2 \
              --vcf "\${_vcf}" \
              --max-alleles 2 \
              --snps-only just-acgt \
              --rm-dup force-first \
              --make-pgen \
              --out "_vcf_part_\${_pidx}"
            echo "_vcf_part_\${_pidx}" >> _vcf_pgen_list.txt
            _pidx=\$(( _pidx + 1 ))
          done
          plink2 \
            --pmerge-list _vcf_pgen_list.txt pfile \
            --make-pgen \
            --out ${meta.id}.raw
        fi
        source_prefix="${meta.id}.raw"
        ;;
      plink1)
        source_bed=\$(find . -maxdepth 1 -name '*.bed' | head -n 1)
        [[ -n "\${source_bed}" ]] || { echo "PLINK1 BED input not found for ${meta.id}" >&2; exit 1; }
        source_prefix=\${source_bed%.bed}
        plink2 \
          --bfile "\${source_prefix}" \
          --snps-only just-acgt \
          --rm-dup force-first \
          --make-pgen \
          --out ${meta.id}.raw
        source_prefix="${meta.id}.raw"
        ;;
      plink2)
        source_pgen=\$(find . -maxdepth 1 -name '*.pgen' | head -n 1)
        [[ -n "\${source_pgen}" ]] || { echo "PLINK2 PGEN input not found for ${meta.id}" >&2; exit 1; }
        source_prefix=\${source_pgen%.pgen}
        plink2 \
          --pfile "\${source_prefix}" \
          --snps-only just-acgt \
          --rm-dup force-first \
          --make-pgen \
          --out ${meta.id}.raw
        source_prefix="${meta.id}.raw"
        ;;
      *)
        echo "Unsupported input format '${meta.input_format}' for ${meta.id}" >&2
        exit 1
        ;;
    esac

    harmonize_input_prefix="\${source_prefix}"
    conversion_method="native"

    if [[ "\${dataset_build}" != "\${target_build}" ]]; then
      if [[ "${params.liftover_enabled}" != "true" ]]; then
        echo "Dataset ${meta.id} build \${dataset_build} differs from \${target_build}, and liftover is disabled." >&2
        exit 1
      fi

      if [[ -n "${liftoverChainArg}" && -n "${liftoverRefArg}" ]]; then
        command -v ${params.crossmap_cmd} >/dev/null 2>&1 || { echo "CrossMap command not found: ${params.crossmap_cmd}" >&2; exit 1; }

        plink2 \
          --pfile "\${source_prefix}" \
          --export vcf bgz id-paste=iid \
          --out ${meta.id}.buildconv

        ${params.crossmap_cmd} vcf "${liftoverChainArg}" ${meta.id}.buildconv.vcf.gz "${liftoverRefArg}" ${meta.id}.buildconv.lifted.vcf

        lifted_vcf=${meta.id}.buildconv.lifted.vcf
        [[ -f "\${lifted_vcf}" ]] || { echo "CrossMap did not produce lifted VCF for ${meta.id}" >&2; exit 1; }

        plink2 \
          --vcf "\${lifted_vcf}" \
          --max-alleles 2 \
          --snps-only just-acgt \
          --rm-dup force-first \
          --make-pgen \
          --out ${meta.id}.buildconv

        harmonize_input_prefix="${meta.id}.buildconv"
        conversion_method="crossmap_vcf"
      elif [[ -n "${liftoverMapDirArg}" && -f "${liftoverMapDirArg}/${meta.id}.map.tsv" ]]; then
        plink2 \
          --pfile "\${source_prefix}" \
          --update-map "${liftoverMapDirArg}/${meta.id}.map.tsv" 2 1 \
          --make-pgen \
          --out ${meta.id}.buildconv

        harmonize_input_prefix="${meta.id}.buildconv"
        conversion_method="plink_update_map"
      else
        echo "Dataset ${meta.id} is \${dataset_build} but target is \${target_build}; provide --liftover_chain and --liftover_ref_fasta (or --liftover_map_dir with ${meta.id}.map.tsv)." >&2
        exit 1
      fi
    fi

    plink2 \
      --pfile "\${harmonize_input_prefix}" \
      ${refFromFaArg} \
      --set-all-var-ids @:#:\\\$r:\\\$a \
      --new-id-max-allele-len 250 \
      --make-pgen \
      --out ${meta.id}.harmonized

    awk 'BEGIN{OFS="\t"} !/^#/ && ((\$4=="A" && \$5=="T") || (\$4=="T" && \$5=="A") || (\$4=="C" && \$5=="G") || (\$4=="G" && \$5=="C")) { print \$3 }' ${meta.id}.harmonized.pvar > ${meta.id}.palindromic.exclude

    plink2 \
      --pfile ${meta.id}.harmonized \
      --exclude ${meta.id}.palindromic.exclude \
      --snps-only just-acgt \
      --make-pgen \
      --out ${meta.id}.norm

    awk '!/^#/ { print \$3 }' ${meta.id}.norm.pvar > ${meta.id}.norm.variants.tsv

    metadata_map_file="${meta.id}.superpop.map.tsv"
    : > "\${metadata_map_file}"

    if [[ -n "${meta.metadata_path}" ]]; then
      metadata_file="\$(basename "${meta.metadata_path}")"
      if [[ ! -f "\${metadata_file}" ]]; then
        echo "Metadata file not staged for ${meta.id}: ${meta.metadata_path}" >&2
        exit 1
      fi

      awk -v sample_col="${meta.metadata_sample_id_col}" -v super_col="${meta.metadata_superpop_col}" -v pop_col="${meta.metadata_pop_col}" '
        function clean(s) {
          gsub(/\r/, "", s)
          gsub(/^[[:space:]]+|[[:space:]]+\$/, "", s)
          return s
        }
        function header_key(h) {
          h=clean(h)
          gsub(/^#/, "", h)
          h=tolower(h)
          gsub(/[[:space:]]+/, "", h)
          return h
        }
        function normalize_header(h) {
          return header_key(h)
        }
        function canon(raw, val) {
          val=toupper(raw)
          gsub(/[^A-Z0-9]+/, "_", val)
          gsub(/^_+|_+\$/, "", val)
          gsub(/_+/, "_", val)

          if (val=="AFR" || val=="EUR" || val=="EAS" || val=="SAS" || val=="AMR") return val

          if (val ~ /^(YRI|LWK|GWD|MSL|ESN|ASW|ACB)\$/) return "AFR"
          if (val ~ /^(CEU|TSI|FIN|GBR|IBS)\$/) return "EUR"
          if (val ~ /^(CHB|JPT|CHS|CDX|KHV)\$/) return "EAS"
          if (val ~ /^(GIH|PJL|BEB|STU|ITU)\$/) return "SAS"
          if (val ~ /^(MXL|PUR|CLM|PEL)\$/) return "AMR"

          if (val ~ /(AFR|AFRICA|SUBSAHARAN)/) return "AFR"
          if (val ~ /(EUR|EUROPE|WEST_EURASIA|NFE|FIN|ASJ|MIDEAST|MIDDLE_EAST|MENA)/) return "EUR"
          if (val ~ /(EAS|EAST_ASIA|CHN|HAN|JPN|KOR|MONG|SEA)/) return "EAS"
          if (val ~ /(SAS|SOUTH_ASIA|CENTRAL_SOUTH_ASIA|CSA|BALOCHI|BRAHUI|BURUSHO|KALASH|PATHAN|PUNJABI|GUJARATI|SINDHI|BENGALI|MAKRANI)/) return "SAS"
          if (val ~ /(AMR|AMERICA|LATIN|MXL|PEL|PUR|CLM|ACB|ASW)/) return "AMR"

          return ""
        }
        BEGIN {
          FS="\t"
          OFS="\t"
        }
        NR==1 {
          for (i=1; i<=NF; i++) {
            key=normalize_header(\$i)
            col[key]=i
          }

          sid_key=normalize_header(sample_col)
          sid_idx=col[sid_key]

          super_key=normalize_header(super_col)
          pop_key=normalize_header(pop_col)
          pop_idx=0
          if (super_key != "") pop_idx=col[super_key]
          if (!pop_idx && pop_key != "") pop_idx=col[pop_key]

          if (!sid_idx) {
            print "Metadata join failed: sample_id_col \"" sample_col "\" not found in " FILENAME > "/dev/stderr"
            exit 2
          }
          if (!pop_idx) {
            print "Metadata join failed: neither superpop_col \"" super_col "\" nor pop_col \"" pop_col "\" found in " FILENAME > "/dev/stderr"
            exit 2
          }
          next
        }
        {
          sid=clean(\$sid_idx)
          mapped=canon(clean(\$pop_idx))
          if (sid != "" && mapped != "") {
            print sid "\t" mapped
          }
        }
      ' "\${metadata_file}" | sort -u > "\${metadata_map_file}"
    fi

    awk -v dataset_superpop="${meta.superpop}" '
      function clean(s) {
        gsub(/\r/, "", s)
        gsub(/^[[:space:]]+|[[:space:]]+\$/, "", s)
        return s
      }
      BEGIN { FS="[[:space:]]+"; OFS="\t"; print "IID","FID","SUPERPOP","DATASET","ASSAY","BUILD" }
      FNR==NR {
        if (NF >= 2) {
          map_id=clean(\$1)
          map_superpop=clean(\$2)
          if (map_id != "" && map_superpop != "") sample_superpop[map_id]=map_superpop
        }
        next
      }
      NR>1 {
        iid=clean(\$1)
        sp=""
        if (iid != "") sp=sample_superpop[iid]
        if (sp == "") sp=dataset_superpop
        print iid, "NA", sp, "${meta.id}", "${meta.assay}", "${meta.build}"
      }
    ' "\${metadata_map_file}" ${meta.id}.norm.psam > ${meta.id}.norm.samples.tsv

    sample_count=\$(( \$(wc -l < ${meta.id}.norm.psam) - 1 ))
    variant_count=\$(wc -l < ${meta.id}.norm.variants.tsv)

    cat <<EOF > ${meta.id}.norm.summary.json
{
  "dataset_id": "${meta.id}",
  "superpop": "${meta.superpop}",
  "assay": "${meta.assay}",
  "source_build": "${meta.build}",
  "target_build": "${params.target_build}",
  "build_conversion": "\${conversion_method}",
  "normalized_sample_count": \${sample_count},
  "normalized_variant_count": \${variant_count},
  "palindromic_filter_applied": true,
  "reference_alignment": "${params.ref_fasta ? 'ref-from-fa' : 'input-build-verified'}"
}
EOF
    """
}