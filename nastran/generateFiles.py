import os
import json

# 👉 修改1：仓库路径
REPO_PATH = "NASTRAN-95"

# 👉 修改2：忽略目录（NASTRAN 特有）
IGNORE_DIRS = {
    ".git",
    "bin",        # 可执行文件
    "demoout",    # 输出示例
    "inp",        # 输入样例
    "__pycache__"
}

# 👉 读取模块定义
with open("modules_index.json") as f:
    module_data = json.load(f)

# 👉 兼容两种格式
if isinstance(module_data, dict):
    modules = module_data.get("modules", [])
else:
    modules = module_data  # 如果是 list

file_index = []

for module in modules:

    # 👉 支持 {"module": "mds", "path": "mds"} 结构
    if isinstance(module, dict):
        module_name = module["module"]
        module_path = module.get("path", module_name)
    else:
        module_name = module
        module_path = module

    abs_module_path = os.path.join(REPO_PATH, module_path)

    if not os.path.exists(abs_module_path):
        print(f"⚠️ Skip (not found): {abs_module_path}")
        continue

    for root, dirs, files in os.walk(abs_module_path):

        # 👉 过滤目录
        dirs[:] = [d for d in dirs if d not in IGNORE_DIRS]

        for file_name in files:

            # 👉 只保留 Fortran 源码（关键优化）
            if not file_name.lower().endswith((".f", ".for", ".f90")):
                continue

            full_path = os.path.join(root, file_name)
            rel_path = os.path.relpath(full_path, REPO_PATH)

            file_index.append({
                "path": rel_path.replace("\\", "/"),
                "module": module_name
            })

# 👉 输出
output = {
    "system": "nastran-95",
    "total_files": len(file_index),
    "files": file_index
}

with open("file_index.json", "w") as f:
    json.dump(output, f, indent=2)

print("✅ Total files:", len(file_index))