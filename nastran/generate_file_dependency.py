import json
import os

print("Loading file_index.json...")

with open("file_index.json", encoding="utf-8") as f:
    file_data = json.load(f)

files = file_data["files"]

# 👉 所有文件路径集合（用于验证）
all_paths = set(f["path"] for f in files)

print("Indexed files:", len(all_paths))


print("Loading dependency_raw.json...")

with open("dependency_raw.json", encoding="utf-8") as f:
    dep_data = json.load(f)

raw_deps = dep_data["dependencies"]


print("Loading subroutine_index.json...")

with open("subroutine_index.json", encoding="utf-8") as f:
    subroutine_index = json.load(f)

# 👉 统一 key 为大写（防止大小写问题）
subroutine_index = {k.upper(): v for k, v in subroutine_index.items()}

print("Indexed subroutines:", len(subroutine_index))


file_dependencies = []

unresolved = 0   # 找不到的 subroutine
self_calls = 0   # 自调用（可选统计）

for dep in raw_deps:

    src = dep["source"]
    target_symbol = dep["target"].upper()

    # 👉 查找 subroutine 定义
    target_path = subroutine_index.get(target_symbol)

    if target_path is None:
        unresolved += 1
        continue

    # 👉 可选：过滤自己调用自己（减少噪声）
    if src == target_path:
        self_calls += 1
        continue

    # 👉 确认目标文件存在
    if target_path not in all_paths:
        # 有些 subroutine 可能来自未纳入模块的文件
        unresolved += 1
        continue

    file_dependencies.append({
        "source": src,
        "target": target_path,
        "type": "call"
    })


print("Resolved dependencies:", len(file_dependencies))
print("Unresolved symbols:", unresolved)
print("Self calls skipped:", self_calls)


output = {
    "system": "nastran-95",
    "total": len(file_dependencies),
    "unresolved": unresolved,
    "dependencies": file_dependencies
}

with open("file_dependency.json", "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2)

print("✅ file_dependency.json generated.")