process EMIT_PROVENANCE {
    tag 'provenance'
    publishDir "${params.outdir}", mode: params.publish_dir_mode, overwrite: true

    input:
    path ref_manifest
    path staged_inputs, stageAs: 'staged_inputs/*'
    path params_yaml
    path high_ld_regions_file
    path layer1_summary
    path layer2_summary
    path layer1_global_weights
    val git_commit
    val nextflow_version
    val run_profile
    val execution_timestamp
    val execution_command
    val container_ref

    output:
    path 'provenance.json', emit: provenance

    stub:
    """
    cat <<EOF > provenance.json
{
  "pipeline": {
    "name": "nf_PopPCA_refgen",
    "git_commit": "${git_commit}",
    "nextflow_version": "${nextflow_version}",
    "execution_timestamp": "${execution_timestamp}",
    "profile": "${run_profile}"
  },
  "environment": {
    "container": "${container_ref}",
    "plink2_version": "stub",
    "plink1_version": "stub",
    "awk_version": "stub",
        "python3_version": "stub",
        "execution_command": "${execution_command}",
        "kernel": "stub",
        "cpu_model": "stub",
        "memory_gb": "stub"
  },
    "config_lineage": {
        "params_yaml": {
            "path": "assets/params.yaml",
            "sha256": "stub"
        },
        "high_ld_regions": {
            "path": "${high_ld_regions_file}",
            "sha256": "stub"
        }
    },
    "matrix_qc_summary": {
        "global": {
            "sample_count": null,
            "variant_count": null
        },
        "sublayers": []
    },
  "input_lineage": [],
  "output_manifests": {
        "layer1_global": "models/layer1/global/global_summary.json",
        "layer2_sublayer": "models/layer2/sub_summary.json",
        "global_weights_sha256": "stub"
  }
}
EOF
    """

    script:
    """
    plink2_version=\$(plink2 --version 2>/dev/null | head -n 1 | tr -d '\\r' || true)
    plink1_version=\$(plink --version 2>/dev/null | head -n 1 | tr -d '\\r' || true)
    awk_version=\$(awk --version 2>/dev/null | head -n 1 | tr -d '\\r' || awk -W version 2>/dev/null | head -n 1 | tr -d '\\r' || true)
    python3_version=\$(python3 --version 2>&1 | tr -d '\\r' || true)

    export REF_MANIFEST="${ref_manifest}"
    export PARAMS_YAML="${params_yaml}"
    export HIGH_LD_REGIONS_FILE="${high_ld_regions_file}"
    export LAYER1_SUMMARY="${layer1_summary}"
    export LAYER2_SUMMARY="${layer2_summary}"
    export LAYER1_GLOBAL_WEIGHTS="${layer1_global_weights}"
    export GIT_COMMIT="${git_commit}"
    export NEXTFLOW_VERSION="${nextflow_version}"
    export RUN_PROFILE="${run_profile}"
    export EXECUTION_TIMESTAMP="${execution_timestamp}"
    export EXECUTION_COMMAND="${execution_command}"
    export CONTAINER_REF="${container_ref}"

    export PLINK2_VERSION="\$plink2_version"
    export PLINK1_VERSION="\$plink1_version"
    export AWK_VERSION="\$awk_version"
    export PYTHON3_VERSION="\$python3_version"

    cat << 'EOF' > generate_provenance.py
import hashlib
import json
import os

ref_manifest = os.environ.get("REF_MANIFEST", "")
layer1_summary = os.environ.get("LAYER1_SUMMARY", "")
layer2_summary = os.environ.get("LAYER2_SUMMARY", "")
params_yaml = os.environ.get("PARAMS_YAML", "")
high_ld_regions_file = os.environ.get("HIGH_LD_REGIONS_FILE", "")
layer1_global_weights = os.environ.get("LAYER1_GLOBAL_WEIGHTS", "")

git_commit = os.environ.get("GIT_COMMIT", "") or "unknown"
nextflow_version = os.environ.get("NEXTFLOW_VERSION", "") or "unknown"
run_profile = os.environ.get("RUN_PROFILE", "") or "unknown"
execution_timestamp = os.environ.get("EXECUTION_TIMESTAMP", "") or "unknown"
execution_command = os.environ.get("EXECUTION_COMMAND", "") or "unknown"
container_ref = os.environ.get("CONTAINER_REF", "") or "unknown"

plink2_version = os.environ.get("PLINK2_VERSION", "")
plink1_version = os.environ.get("PLINK1_VERSION", "")
awk_version = os.environ.get("AWK_VERSION", "")
python3_version = os.environ.get("PYTHON3_VERSION", "")

def sha256_file(path: str):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def first_non_empty(lines):
    for item in lines:
        if item and item.strip():
            return item.strip()
    return "unknown"

def get_uname_sr() -> str:
    try:
        return first_non_empty([os.popen("uname -sr 2>/dev/null").read()])
    except Exception:
        return "unknown"

def get_cpu_model() -> str:
    try:
        lscpu_line = first_non_empty(
            [line.split(":", 1)[1].strip() for line in os.popen("lscpu 2>/dev/null").read().splitlines() if line.lower().startswith("model name:")]
        )
        if lscpu_line != "unknown":
            return lscpu_line
    except Exception:
        pass
    try:
        with open("/proc/cpuinfo", "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.lower().startswith("model name") and ":" in line:
                    return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return "unknown"

def get_memory_gb() -> dict:
    metrics = {"total_gb": None, "used_gb": None, "free_gb": None}
    try:
        output = os.popen("free -g 2>/dev/null").read().splitlines()
        for line in output:
            if line.lower().startswith("mem:"):
                parts = line.split()
                if len(parts) >= 4:
                    metrics["total_gb"] = int(parts[1])
                    metrics["used_gb"] = int(parts[2])
                    metrics["free_gb"] = int(parts[3])
                    return metrics
    except Exception:
        pass

    try:
        meminfo = {}
        with open("/proc/meminfo", "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if ":" not in line:
                    continue
                key, val = line.split(":", 1)
                meminfo[key.strip()] = val.strip()

        total_kb = int(meminfo.get("MemTotal", "0 kB").split()[0])
        free_kb = int(meminfo.get("MemAvailable", meminfo.get("MemFree", "0 kB")).split()[0])
        used_kb = max(total_kb - free_kb, 0)
        if total_kb > 0:
            metrics["total_gb"] = total_kb // (1024 * 1024)
            metrics["used_gb"] = used_kb // (1024 * 1024)
            metrics["free_gb"] = free_kb // (1024 * 1024)
    except Exception:
        pass

    return metrics

def load_json_or_empty(path: str):
    if not path or not os.path.exists(path) or not os.path.isfile(path):
        return {}
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            return json.load(fh)
    except Exception:
        return {}

def resolve_staged_or_source(path_value: str):
    filename = os.path.basename(path_value)
    staged_path = os.path.join("staged_inputs", filename)
    if os.path.exists(filename) and os.path.isfile(filename):
        return filename
    if os.path.exists(staged_path) and os.path.isfile(staged_path):
        return staged_path
    if os.path.exists(path_value) and os.path.isfile(path_value):
        return path_value
    return None

# Minimal, robust parser for the manifest shape used by this pipeline.
datasets = []
current = None
in_input_paths = False
paths_indent = -1
in_metadata = False
metadata_indent = -1

with open(ref_manifest, "r", encoding="utf-8", errors="replace") as fh:
    for raw in fh:
        line = raw.rstrip('\\r\\n')
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        indent = len(line) - len(line.lstrip(" "))

        if stripped.startswith("- id:"):
            dataset_id = stripped.split(":", 1)[1].strip().strip('"').strip("'")
            if not dataset_id:
                continue
            if current:
                datasets.append(current)
            current = {
                "dataset_id": dataset_id,
                "paths": []
            }
            in_input_paths = False
            paths_indent = -1
            in_metadata = False
            metadata_indent = -1
            continue

        if current is None:
            continue

        if stripped == "input:":
            in_input_paths = False
            paths_indent = -1
            continue

        if stripped == "paths:":
            in_input_paths = True
            paths_indent = indent
            continue

        if in_input_paths:
            if indent <= paths_indent:
                # End of the paths block at dedent.
                in_input_paths = False
                paths_indent = -1
            elif stripped.startswith("- "):
                parsed_item = stripped[2:].strip()
                if parsed_item.startswith("path:"):
                    parsed_path = parsed_item.split(":", 1)[1].strip().strip('"').strip("'")
                else:
                    parsed_path = parsed_item.strip().strip('"').strip("'")
                if parsed_path:
                    current["paths"].append(parsed_path)
                continue
            else:
                # Continuation lines inside a path object (e.g., expected_sha256).
                continue

        if stripped == "metadata:":
            in_metadata = True
            metadata_indent = indent
            continue

        if in_metadata:
            if indent <= metadata_indent:
                in_metadata = False
            else:
                if stripped.startswith("path:"):
                    parsed_meta_path = stripped.split(":", 1)[1].strip()
                    if parsed_meta_path and parsed_meta_path != "":
                        parsed_meta_path = parsed_meta_path.strip('"').strip("'")
                        current["paths"].append(parsed_meta_path)
                    continue

if current:
    datasets.append(current)

input_lineage = []
seen = set()
for ds in datasets:
    dsid = ds["dataset_id"]
    for p in ds.get("paths", []):
        key = (dsid, p)
        if key in seen:
            continue
        seen.add(key)

        target = resolve_staged_or_source(p)

        if target:
            try:
                digest = sha256_file(target)
                input_lineage.append({
                    "dataset_id": dsid,
                    "file_path": p,
                    "sha256": digest
                })
            except Exception:
                pass

layer1_payload = load_json_or_empty(layer1_summary)
layer2_payload = load_json_or_empty(layer2_summary)

matrix_qc_summary = {
    "global": {
        "sample_count": layer1_payload.get("sample_count"),
        "variant_count": layer1_payload.get("variant_count")
    },
    "sublayers": []
}

for model in (layer2_payload.get("models") or []):
    if not isinstance(model, dict):
        continue
    matrix_qc_summary["sublayers"].append({
        "superpop": model.get("superpop"),
        "sample_count": model.get("sample_count"),
        "variant_count": model.get("variant_count")
    })

global_weights_sha256 = None
if layer1_global_weights and os.path.exists(layer1_global_weights) and os.path.isfile(layer1_global_weights):
    try:
        global_weights_sha256 = sha256_file(layer1_global_weights)
    except Exception:
        global_weights_sha256 = None

config_lineage = {
    "params_yaml": {
        "path": "assets/params.yaml",
        "sha256": None
    },
    "high_ld_regions": {
        "path": high_ld_regions_file,
        "sha256": None
    }
}

if params_yaml and os.path.exists(params_yaml) and os.path.isfile(params_yaml):
    try:
        config_lineage["params_yaml"]["sha256"] = sha256_file(params_yaml)
    except Exception:
        config_lineage["params_yaml"]["sha256"] = None

if high_ld_regions_file and os.path.exists(high_ld_regions_file) and os.path.isfile(high_ld_regions_file):
    try:
        config_lineage["high_ld_regions"]["sha256"] = sha256_file(high_ld_regions_file)
    except Exception:
        config_lineage["high_ld_regions"]["sha256"] = None

environment = {
    "container": container_ref,
    "plink2_version": plink2_version,
    "plink1_version": plink1_version,
    "awk_version": awk_version,
    "python3_version": python3_version,
    "execution_command": execution_command,
    "kernel": get_uname_sr(),
    "cpu_model": get_cpu_model(),
    "memory_gb": get_memory_gb()
}

provenance = {
    "pipeline": {
        "name": "nf_PopPCA_refgen",
        "git_commit": git_commit,
        "nextflow_version": nextflow_version,
        "execution_timestamp": execution_timestamp,
        "profile": run_profile
    },
    "environment": environment,
    "config_lineage": config_lineage,
    "matrix_qc_summary": matrix_qc_summary,
    "input_lineage": input_lineage,
    "output_manifests": {
        "layer1_global": "models/layer1/global/global_summary.json",
        "layer2_sublayer": "models/layer2/sub_summary.json",
        "global_weights_sha256": global_weights_sha256
    }
}

with open("provenance.json", "w", encoding="utf-8") as out:
    json.dump(provenance, out, indent=2)
    out.write("\\n")
EOF

    python3 generate_provenance.py
    """
}
