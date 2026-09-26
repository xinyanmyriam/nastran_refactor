import os
import json
import re

REPO_PATH = "NASTRAN-95"

# 👉 Fortran 关键字（忽略大小写）
call_pattern = re.compile(r'\bCALL\s+([A-Z0-9_]+)', re.IGNORECASE)
subroutine_pattern = re.compile(r'\bSUBROUTINE\s+([A-Z0-9_]+)', re.IGNORECASE)

with open("file_index.json", encoding="utf-8") as f:
    data = json.load(f)

files = data["files"]

dependencies = []
subroutine_defs = {}  # 👉 symbol -> file（后续Step6要用）

for file_info in files:

    path = file_info["path"]

    # 👉 只分析 Fortran 文件
    if not path.lower().endswith((".f", ".for", ".f90")):
        continue

    full_path = os.path.join(REPO_PATH, path)

    if not os.path.exists(full_path):
        continue

    try:
        with open(full_path, "r", encoding="utf-8", errors="ignore") as f:

            for line in f:

                # 👉 去掉注释（Fortran老格式：C 或 * 开头）
                stripped = line.strip()
                if stripped.startswith("C") or stripped.startswith("*") or stripped.startswith("!"):
                    continue

                # =========================
                # 1️⃣ 解析 CALL（依赖）
                # =========================
                call_matches = call_pattern.findall(line)
                for target in call_matches:
                    dependencies.append({
                        "source": path,
                        "target": target.upper(),
                        "type": "call"
                    })

                # =========================
                # 2️⃣ 解析 SUBROUTINE（定义）
                # =========================
                sub_matches = subroutine_pattern.findall(line)
                for name in sub_matches:
                    subroutine_defs[name.upper()] = path

    except Exception as e:
        print(f"⚠️ Error reading {path}: {e}")
        continue


# 👉 输出 dependency_raw
output = {
    "system": "nastran-95",
    "total_dependencies": len(dependencies),
    "dependencies": dependencies
}

with open("dependency_raw.json", "w", encoding="utf-8") as f:
    json.dump(output, f, indent=2)

# 👉 额外输出（关键🔥）：symbol -> file 映射
with open("subroutine_index.json", "w", encoding="utf-8") as f:
    json.dump(subroutine_defs, f, indent=2)

print("✅ Total dependencies:", len(dependencies))
print("✅ Total subroutines:", len(subroutine_defs))