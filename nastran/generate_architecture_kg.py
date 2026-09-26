import json

print("Loading file_index.json...")

with open("file_index.json", encoding="utf-8") as f:
    file_data = json.load(f)

files = file_data["files"]

print("Files:", len(files))


print("Loading module_dependency.json...")

with open("module_dependency.json", encoding="utf-8") as f:
    mod_data = json.load(f)

module_deps = mod_data["dependencies"]

print("Module dependencies:", len(module_deps))


kg_edges = []


# -------------------------
# 1️⃣ file -> module
# -------------------------
for f in files:

    kg_edges.append({
        "head": f["path"],
        "head_type": "file",          # ⭐新增
        "relation": "belongs_to",
        "tail": f["module"],
        "tail_type": "module",        # ⭐新增
        "system": "nastran-95"        # ⭐新增
    })


# -------------------------
# 2️⃣ module -> module
# -------------------------
for dep in module_deps:

    kg_edges.append({
        "head": dep["source_module"],
        "head_type": "module",
        "relation": "depends_on",
        "tail": dep["target_module"],
        "tail_type": "module",
        "weight": dep["weight"],
        "system": "nastran-95"
    })


print("Total KG edges:", len(kg_edges))


with open("architecture_kg.jsonl", "w", encoding="utf-8") as f:

    for edge in kg_edges:
        f.write(json.dumps(edge) + "\n")


print("✅ architecture_kg.jsonl generated successfully.")