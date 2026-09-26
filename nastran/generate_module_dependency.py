import json

print("Loading file_index.json...")

with open("file_index.json", encoding="utf-8") as f:
    file_data = json.load(f)

files = file_data["files"]

# 👉 file path -> module
file_to_module = {}

for f in files:
    file_to_module[f["path"]] = f["module"]

print("Indexed files:", len(file_to_module))


print("Loading file_dependency.json...")

with open("file_dependency.json", encoding="utf-8") as f:
    dep_data = json.load(f)

deps = dep_data["dependencies"]

print("File dependencies:", len(deps))


# 👉 module dependency counter
module_counter = {}

missing_target = 0  # 统计异常

for dep in deps:

    src_file = dep["source"]
    tgt_file = dep["target"]   # ✅ NASTRAN是文件路径

    src_module = file_to_module.get(src_file)
    tgt_module = file_to_module.get(tgt_file)  # ⭐关键修改

    if not src_module or not tgt_module:
        missing_target += 1
        continue

    # 👉 忽略同模块依赖（可选）
    if src_module == tgt_module:
        continue

    key = (src_module, tgt_module)

    if key not in module_counter:
        module_counter[key] = 0

    module_counter[key] += 1


# 👉 转换为列表
module_dependencies = []

for (src, tgt), count in module_counter.items():

    module_dependencies.append({
        "source_module": src,
        "target_module": tgt,
        "weight": count
    })


output = {
    "system": "nastran-95",
    "total": len(module_dependencies),
    "missing_target": missing_target,
    "dependencies": module_dependencies
}

with open("module_dependency.json", "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2)

print("Module dependencies:", len(module_dependencies))
print("Missing target mappings:", missing_target)
print("✅ module_dependency.json generated.")