import json
import random
from collections import defaultdict

random.seed(42)

# =========================
# Load modules
# =========================

print("Loading module_index...")

with open("modules_index.json", encoding="utf-8") as f:
    module_data = json.load(f)

modules = list(set(module_data["modules"]))

print("Modules:", len(modules))


# =========================
# Load file index
# =========================

print("Loading file_index...")

with open("file_index.json", encoding="utf-8") as f:
    file_data = json.load(f)

files = file_data["files"]

print("Files:", len(files))


# =========================
# Load module dependencies
# =========================

print("Loading module dependencies...")

with open("module_dependency.json", encoding="utf-8") as f:
    dep_data = json.load(f)

deps = dep_data["dependencies"]

dep_pairs = set()
outgoing = defaultdict(set)
incoming = defaultdict(set)

for d in deps:
    src = d["source_module"]
    tgt = d["target_module"]

    if src != tgt:
        dep_pairs.add((src, tgt))
        outgoing[src].add(tgt)
        incoming[tgt].add(src)

print("Dependency pairs:", len(dep_pairs))


# =========================
# Component Identification
# =========================

print("Generating component identification...")

component_tasks = []

for f in files:

    true_module = f["module"]
    file_path = f["path"]

    candidates = {true_module}

    while len(candidates) < 4:
        candidates.add(random.choice(modules))

    candidates = list(candidates)
    random.shuffle(candidates)

    component_tasks.append({
        "file": file_path,
        "candidates": candidates,
        "label": true_module
    })

component_tasks = component_tasks[:4000]


# =========================
# Dependency Prediction（论文级🔥）
# =========================

print("Generating dependency prediction (hard negatives)...")

dependency_tasks = []

dep_list = list(dep_pairs)

# -------- 正样本 --------
pos_samples = random.sample(dep_list, min(2000, len(dep_list)))

for src, tgt in pos_samples:
    dependency_tasks.append({
        "module_A": src,
        "module_B": tgt,
        "label": 1
    })


# -------- Hard Negative --------
neg_samples = set()

modules_list = list(modules)

def sample_hard_negative():
    """
    Hard negative strategy:
    1. same source, different target
    2. same target, different source
    3. two-hop neighbors（结构近邻）
    """

    strategy = random.choice([1, 2, 3])

    # -------- Strategy 1 --------
    if strategy == 1:
        src = random.choice(list(outgoing.keys()))
        candidates = list(set(modules_list) - outgoing[src] - {src})
        if candidates:
            tgt = random.choice(candidates)
            return (src, tgt)

    # -------- Strategy 2 --------
    if strategy == 2:
        tgt = random.choice(list(incoming.keys()))
        candidates = list(set(modules_list) - incoming[tgt] - {tgt})
        if candidates:
            src = random.choice(candidates)
            return (src, tgt)

    # -------- Strategy 3 --------
    if strategy == 3:
        src = random.choice(modules_list)
        tgt = random.choice(modules_list)

        if src != tgt and (src, tgt) not in dep_pairs:
            return (src, tgt)

    return None


while len(neg_samples) < len(pos_samples):

    pair = sample_hard_negative()

    if pair and pair not in dep_pairs:
        neg_samples.add(pair)


for m1, m2 in neg_samples:
    dependency_tasks.append({
        "module_A": m1,
        "module_B": m2,
        "label": 0
    })


# =========================
# Architecture Reconstruction
# =========================

print("Generating architecture reconstruction...")

dep_list = list(dep_pairs)

sample_size = min(2000, len(dep_list))

reconstruction_tasks = random.sample(dep_list, sample_size)

reconstruction_tasks = [
    {
        "source_module": src,
        "target_module": tgt
    }
    for src, tgt in reconstruction_tasks
]


# =========================
# Save JSONL
# =========================

def save_jsonl(data, path):
    with open(path, "w", encoding="utf-8") as f:
        for item in data:
            f.write(json.dumps(item) + "\n")


save_jsonl(component_tasks, "component_identification.jsonl")
save_jsonl(dependency_tasks, "dependency_prediction.jsonl")
save_jsonl(reconstruction_tasks, "architecture_reconstruction.jsonl")


# =========================
# Stats
# =========================

print("\nTasks generated successfully.")
print("Component:", len(component_tasks))
print("Dependency:", len(dependency_tasks))
print("Reconstruction:", len(reconstruction_tasks))